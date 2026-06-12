#!/usr/bin/env python3
import argparse
import csv
import glob
import hashlib
import json
import os
import sys
from pathlib import Path


PRODUCTION_SAMPLE_COLUMNS = [
    "samplesheet_schema_version",
    "run_id",
    "minknow_run_folder",
    "barcode_id",
    "sample_id",
    "specimen_id",
    "sample_type",
    "control_type",
    "expected_host_scientific_name",
    "expected_host_taxid",
    "expected_marker_result",
    "species",
    "sibling_species",
    "feeding_status",
    "collection_date",
    "collection_time",
    "collection_location",
    "collection_region",
    "bioclimatic_zone",
    "latitude",
    "longitude",
    "collection_context",
    "collection_method",
    "specimen_sex",
    "extraction_batch",
    "pcr_batch",
    "library_batch",
    "barcode_kit",
    "flowcell",
    "basecalling_model",
]

REQUIRED_PRIMER_COLUMNS = ["Gene", "Forward_Primer", "Reverse_Primer", "Size"]
REQUIRED_TAXONOMY_COLUMNS = ["seqid", "taxid", "scientific_name", "rank", "source_accession", "sequence_md5"]
REQUIRED_TARGET_COLUMNS = [
    "priority",
    "group",
    "scientific_name",
    "common_name",
    "ncbi_taxid",
    "accession",
    "markers_needed",
    "ghana_relevance",
    "status",
]


def str_to_bool(value):
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def read_csv(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames or [], list(reader)


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def fasta_records(path):
    records = {}
    name = None
    parts = []
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seq = "".join(parts).upper()
                    records[name] = {
                        "length": len(seq),
                        "sequence_md5": hashlib.md5(seq.encode()).hexdigest(),
                    }
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if name is not None:
            seq = "".join(parts).upper()
            records[name] = {
                "length": len(seq),
                "sequence_md5": hashlib.md5(seq.encode()).hexdigest(),
            }
    return records


def path_exists(path):
    return bool(path) and Path(path).expanduser().exists()


def blast_prefix_exists(prefix):
    if not prefix:
        return False
    prefix = str(Path(prefix).expanduser())
    return bool(glob.glob(prefix + ".*"))


def write_tsv(path, rows, fieldnames):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def add_check(checks, name, status, severity, message, evidence=""):
    checks.append(
        {
            "check": name,
            "status": status,
            "severity": severity,
            "message": message,
            "evidence": evidence,
        }
    )


def main():
    parser = argparse.ArgumentParser(description="Strict production preflight for the HÆMA blood-meal pipeline")
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--primers", required=True)
    parser.add_argument("--raw-data-dir", required=True)
    parser.add_argument("--reference-fasta", required=True)
    parser.add_argument("--curated-reference-metadata", default="")
    parser.add_argument("--reference-targets", default="")
    parser.add_argument("--blast-db", default="")
    parser.add_argument("--fallback-blast-db", default="")
    parser.add_argument("--taxdump-dir", default="")
    parser.add_argument("--barcode-kit", default="")
    parser.add_argument("--basecalling-model", default="")
    parser.add_argument("--flowcell", default="")
    parser.add_argument("--demux-strategy", default="")
    parser.add_argument("--production-mode", default="false")
    parser.add_argument("--enable-medaka", default="false")
    parser.add_argument("--medaka-model", default="")
    parser.add_argument("--enable-r-outputs", default="false")
    parser.add_argument("--strict-bioconductor", default="false")
    parser.add_argument("--r-container", default="")
    parser.add_argument("--python-container", default="")
    parser.add_argument("--medaka-container", default="")
    parser.add_argument("--output-dir", default=".")
    args = parser.parse_args()

    production_mode = str_to_bool(args.production_mode)
    enable_medaka = str_to_bool(args.enable_medaka)
    enable_r_outputs = str_to_bool(args.enable_r_outputs)
    strict_bioc = str_to_bool(args.strict_bioconductor)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    checks = []

    for name, path in [
        ("samplesheet_exists", args.samplesheet),
        ("primers_exists", args.primers),
        ("raw_data_dir_exists", args.raw_data_dir),
        ("reference_fasta_exists", args.reference_fasta),
    ]:
        exists = path_exists(path)
        add_check(checks, name, "pass" if exists else "fail", "error", f"{path}", "")

    sample_columns, sample_rows = ([], [])
    if path_exists(args.samplesheet):
        sample_columns, sample_rows = read_csv(args.samplesheet)
        missing = [col for col in PRODUCTION_SAMPLE_COLUMNS if col not in sample_columns]
        status = "pass" if not missing else ("fail" if production_mode else "warn")
        severity = "error" if production_mode and missing else "warning"
        add_check(
            checks,
            "samplesheet_production_columns",
            status,
            severity,
            "Production samplesheet column check",
            ",".join(missing),
        )

        positive_missing = []
        negative_missing = []
        barcode_mismatch = []
        basecalling_mismatch = []
        coordinate_missing = []
        for line_no, row in enumerate(sample_rows, start=2):
            sample_type = (row.get("sample_type") or "").strip().lower()
            control_type = (row.get("control_type") or "").strip().lower()
            if sample_type != "sample" and "control_type" in sample_columns and not control_type:
                negative_missing.append(str(line_no))
            if "positive" in sample_type:
                if not (row.get("expected_host_scientific_name") or "").strip() or not (
                    row.get("expected_host_taxid") or ""
                ).strip():
                    positive_missing.append(str(line_no))
            if sample_type == "sample":
                if "latitude" in sample_columns and "longitude" in sample_columns:
                    if not (row.get("latitude") or "").strip() or not (row.get("longitude") or "").strip():
                        coordinate_missing.append(str(line_no))
            if args.barcode_kit and "barcode_kit" in sample_columns:
                row_kit = (row.get("barcode_kit") or "").strip()
                if row_kit and row_kit != args.barcode_kit:
                    barcode_mismatch.append(f"{line_no}:{row_kit}")
            if args.basecalling_model and "basecalling_model" in sample_columns:
                row_model = (row.get("basecalling_model") or "").strip()
                if row_model and row_model != args.basecalling_model:
                    basecalling_mismatch.append(f"{line_no}:{row_model}")

        add_check(
            checks,
            "positive_control_expectations",
            "pass" if not positive_missing else ("fail" if production_mode else "warn"),
            "error" if production_mode and positive_missing else "warning",
            "Positive controls require expected_host_scientific_name and expected_host_taxid",
            ",".join(positive_missing),
        )
        add_check(
            checks,
            "control_type_populated",
            "pass" if not negative_missing else ("fail" if production_mode else "warn"),
            "error" if production_mode and negative_missing else "warning",
            "Controls require a populated control_type value",
            ",".join(negative_missing),
        )
        add_check(
            checks,
            "sample_coordinates_populated",
            "pass" if not coordinate_missing else ("fail" if production_mode else "warn"),
            "error" if production_mode and coordinate_missing else "warning",
            "Sample rows require latitude and longitude in production metadata",
            ",".join(coordinate_missing),
        )
        add_check(
            checks,
            "barcode_kit_consistency",
            "pass" if not barcode_mismatch else ("fail" if production_mode else "warn"),
            "error" if production_mode and barcode_mismatch else "warning",
            f"Expected barcode kit: {args.barcode_kit}",
            ",".join(barcode_mismatch[:20]),
        )
        add_check(
            checks,
            "basecalling_model_consistency",
            "pass" if not basecalling_mismatch else ("fail" if production_mode else "warn"),
            "error" if production_mode and basecalling_mismatch else "warning",
            f"Expected basecalling model: {args.basecalling_model}",
            ",".join(basecalling_mismatch[:20]),
        )

    if path_exists(args.primers):
        primer_columns, primer_rows = read_csv(args.primers)
        missing = [col for col in REQUIRED_PRIMER_COLUMNS if col not in primer_columns]
        bad_sizes = []
        for line_no, row in enumerate(primer_rows, start=2):
            try:
                int(float(str(row.get("Size", "")).strip()))
            except ValueError:
                bad_sizes.append(str(line_no))
        add_check(
            checks,
            "primer_csv_schema",
            "pass" if not missing and not bad_sizes else "fail",
            "error",
            "Primer file must define Gene, Forward_Primer, Reverse_Primer, and numeric Size",
            f"missing={','.join(missing)}; bad_size_lines={','.join(bad_sizes)}",
        )

    add_check(
        checks,
        "minknow_demux_strategy",
        "pass" if args.demux_strategy == "pre_demultiplexed_minknow_trusted_folder" else (
            "fail" if production_mode else "warn"
        ),
        "error" if production_mode else "warning",
        "Production runs trust the containing MinKNOW barcode folder",
        args.demux_strategy,
    )
    add_check(
        checks,
        "barcode_kit_configured",
        "pass" if args.barcode_kit else "fail",
        "error",
        "A barcode kit must be recorded",
        args.barcode_kit,
    )
    add_check(
        checks,
        "basecalling_model_configured",
        "pass" if args.basecalling_model else "fail",
        "error",
        "A basecalling model must be recorded",
        args.basecalling_model,
    )

    blast_ok = blast_prefix_exists(args.blast_db) or blast_prefix_exists(args.fallback_blast_db)
    add_check(
        checks,
        "local_blast_database_exists",
        "pass" if blast_ok else ("fail" if production_mode else "warn"),
        "error" if production_mode else "warning",
        "At least one local BLAST database prefix must exist",
        f"blast_db={args.blast_db}; fallback_blast_db={args.fallback_blast_db}",
    )

    fasta = {}
    if path_exists(args.reference_fasta):
        fasta = fasta_records(args.reference_fasta)

    if args.curated_reference_metadata:
        if path_exists(args.curated_reference_metadata):
            tax_columns, tax_rows = read_tsv(args.curated_reference_metadata)
            missing = [col for col in REQUIRED_TAXONOMY_COLUMNS if col not in tax_columns]
            seqids = {row.get("seqid", "") for row in tax_rows}
            missing_seqids = sorted(set(fasta) - seqids)
            md5_mismatch = []
            for row in tax_rows:
                seqid = row.get("seqid", "")
                expected_md5 = (row.get("sequence_md5") or "").strip()
                if seqid in fasta and expected_md5 and expected_md5 != fasta[seqid]["sequence_md5"]:
                    md5_mismatch.append(seqid)
            status = "pass" if not missing and not missing_seqids and not md5_mismatch else (
                "fail" if production_mode else "warn"
            )
            add_check(
                checks,
                "curated_taxonomy_sidecar_complete",
                status,
                "error" if production_mode and status == "fail" else "warning",
                "Taxonomy sidecar must cover each reference FASTA record with taxid and checksum",
                f"missing_cols={','.join(missing)}; missing_seqids={','.join(missing_seqids[:20])}; md5_mismatch={','.join(md5_mismatch[:20])}",
            )
        else:
            add_check(
                checks,
                "curated_taxonomy_sidecar_complete",
                "fail" if production_mode else "warn",
                "error" if production_mode else "warning",
                "Curated reference metadata sidecar not found",
                args.curated_reference_metadata,
            )
    else:
        add_check(
            checks,
            "curated_taxonomy_sidecar_complete",
            "fail" if production_mode else "warn",
            "error" if production_mode else "warning",
            "Curated reference metadata sidecar is not configured",
            "",
        )

    if args.reference_targets:
        if path_exists(args.reference_targets):
            target_columns, target_rows = read_tsv(args.reference_targets)
            missing = [col for col in REQUIRED_TARGET_COLUMNS if col not in target_columns]
            incomplete_ready = []
            needs_accession = 0
            ready = 0
            for line_no, row in enumerate(target_rows, start=2):
                status = (row.get("status") or "").strip().lower()
                if status == "ready":
                    ready += 1
                    if not (row.get("accession") or "").strip() or not (row.get("ncbi_taxid") or "").strip():
                        incomplete_ready.append(str(line_no))
                elif status == "needs_accession":
                    needs_accession += 1
            add_check(
                checks,
                "ghana_reference_target_manifest",
                "pass" if not missing and not incomplete_ready else ("fail" if production_mode else "warn"),
                "error" if production_mode and (missing or incomplete_ready) else "warning",
                "Reference target manifest may include placeholders but ready rows need accession and taxid",
                f"missing_cols={','.join(missing)}; ready={ready}; needs_accession={needs_accession}; incomplete_ready_lines={','.join(incomplete_ready)}",
            )
        else:
            add_check(
                checks,
                "ghana_reference_target_manifest",
                "fail" if production_mode else "warn",
                "error" if production_mode else "warning",
                "Reference target manifest not found",
                args.reference_targets,
            )

    taxdump_nodes = Path(args.taxdump_dir).expanduser() / "nodes.dmp" if args.taxdump_dir else None
    taxdump_names = Path(args.taxdump_dir).expanduser() / "names.dmp" if args.taxdump_dir else None
    taxdump_ok = bool(args.taxdump_dir and taxdump_nodes.exists() and taxdump_names.exists())
    add_check(
        checks,
        "taxdump_nodes_names_exist",
        "pass" if taxdump_ok else ("fail" if production_mode else "warn"),
        "error" if production_mode else "warning",
        "Full production LCA requires NCBI taxdump nodes.dmp and names.dmp",
        args.taxdump_dir,
    )

    add_check(
        checks,
        "r_bioconductor_production_config",
        "pass" if (enable_r_outputs and strict_bioc and args.r_container) else ("fail" if production_mode else "warn"),
        "error" if production_mode else "warning",
        "Production mode requires formal R outputs with strict Bioconductor packages and an R container",
        f"enable_r_outputs={enable_r_outputs}; strict_bioconductor={strict_bioc}; r_container={args.r_container}",
    )

    add_check(
        checks,
        "medaka_production_config",
        "pass" if (enable_medaka and args.medaka_model and args.medaka_container) else (
            "fail" if production_mode else "warn"
        ),
        "error" if production_mode else "warning",
        "Production mode requires Medaka polishing with a configured model and container",
        f"enable_medaka={enable_medaka}; medaka_model={args.medaka_model}; medaka_container={args.medaka_container}",
    )

    errors = [
        check
        for check in checks
        if check["status"] == "fail" and (production_mode or check["severity"] == "error")
    ]
    warnings = [check for check in checks if check["status"] == "warn"]
    report = {
        "production_mode": production_mode,
        "pass": not errors,
        "checks": checks,
        "errors": errors,
        "warnings": warnings,
        "summary": {
            "samplesheet_rows": len(sample_rows),
            "reference_fasta_records": len(fasta),
            "barcode_kit": args.barcode_kit,
            "flowcell": args.flowcell,
            "basecalling_model": args.basecalling_model,
            "demux_strategy": args.demux_strategy,
            "medaka_model": args.medaka_model,
            "taxdump_dir": args.taxdump_dir,
            "blast_db": args.blast_db,
            "fallback_blast_db": args.fallback_blast_db,
            "python_container": args.python_container,
            "r_container": args.r_container,
            "medaka_container": args.medaka_container,
        },
    }

    with (outdir / "production_preflight_report.json").open("w") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    write_tsv(outdir / "production_preflight_summary.tsv", checks, ["check", "status", "severity", "message", "evidence"])

    for check in warnings[:30]:
        print(f"WARNING: {check['check']}: {check['message']} [{check['evidence']}]", file=sys.stderr)
    if errors:
        for check in errors:
            print(f"ERROR: {check['check']}: {check['message']} [{check['evidence']}]", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
