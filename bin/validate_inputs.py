#!/usr/bin/env python3
import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


REQUIRED_SAMPLE_COLUMNS = [
    "run_id",
    "barcode_id",
    "sample_id",
    "sample_type",
]

RECOMMENDED_METADATA_COLUMNS = [
    "specimen_id",
    "species",
    "sibling_species",
    "feeding_status",
    "collection_date",
    "collection_time",
    "collection_location",
    "bioclimatic_zone",
    "collection_region",
    "collection_cordinates",
    "collection_context",
    "collection_method",
    "specimen_sex",
]

PRODUCTION_METADATA_COLUMNS = [
    "samplesheet_schema_version",
    "minknow_run_folder",
    "control_type",
    "expected_host_scientific_name",
    "expected_host_taxid",
    "expected_marker_result",
    "latitude",
    "longitude",
    "extraction_batch",
    "pcr_batch",
    "library_batch",
    "barcode_kit",
    "flowcell",
    "basecalling_model",
]

REQUIRED_PRIMER_COLUMNS = ["Gene", "Forward_Primer", "Reverse_Primer", "Size"]


def str_to_bool(value):
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def sanitize(value):
    value = str(value or "").strip()
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("_") or "missing"


def read_csv(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        return reader.fieldnames or [], rows


def write_tsv(path, rows, fieldnames):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_marker_windows(spec):
    windows = {}
    if not spec:
        return windows
    for item in str(spec).split(","):
        item = item.strip()
        if not item:
            continue
        parts = item.split(":")
        if len(parts) != 3:
            raise ValueError(f"Invalid marker window '{item}'. Expected marker:min:max")
        marker, min_len, max_len = parts
        windows[marker.strip()] = (int(min_len), int(max_len))
    return windows


def has_extension(path, extensions):
    name = path.name.lower()
    return any(name.endswith(ext.lower()) for ext in extensions)


def discover_fastqs(raw_data_dir, extensions):
    raw = Path(raw_data_dir).expanduser().resolve()
    if not raw.exists():
        return [], [f"Raw data directory does not exist: {raw}"]
    fastqs = [p for p in raw.rglob("*") if p.is_file() and has_extension(p, extensions)]
    return sorted(fastqs), []


def infer_run_barcode(path, raw_data_dir):
    raw = Path(raw_data_dir).expanduser().resolve()
    try:
        rel = path.resolve().relative_to(raw)
    except ValueError:
        return None, None
    run_id = rel.parts[0] if rel.parts else None
    barcode_id = next((part for part in rel.parts if part.startswith("barcode")), None)
    return run_id, barcode_id


def fastqs_for_sample(fastqs_by_run_barcode, run_id, barcode_id):
    return sorted(fastqs_by_run_barcode.get((run_id, barcode_id), []))


def validate_coordinate(value):
    value = (value or "").strip()
    if not value:
        return False
    if "," in value:
        parts = [p.strip() for p in value.split(",", 1)]
    elif ";" in value:
        parts = [p.strip() for p in value.split(";", 1)]
    else:
        # Existing project metadata uses values like 8.785-1.535; this is
        # ambiguous because the hyphen can be a separator or a negative sign.
        return False
    try:
        lat, lon = float(parts[0]), float(parts[1])
    except ValueError:
        return False
    return -90 <= lat <= 90 and -180 <= lon <= 180


def main():
    parser = argparse.ArgumentParser(description="Validate HÆMA blood-meal pipeline inputs")
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--primers", required=True)
    parser.add_argument("--raw-data-dir", required=True)
    parser.add_argument("--fastq-extensions", default=".fastq.gz,.fq.gz,.fastq,.fq")
    parser.add_argument("--marker-windows", default="")
    parser.add_argument("--marker-window-padding", type=int, default=75)
    parser.add_argument("--require-fastqs", default="true")
    parser.add_argument("--strict-metadata", default="false")
    parser.add_argument("--strict-barcode-filenames", default="false")
    parser.add_argument("--production-mode", default="false")
    parser.add_argument("--barcode-kit", default="")
    parser.add_argument("--basecalling-model", default="")
    parser.add_argument("--output-dir", default=".")
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    require_fastqs = str_to_bool(args.require_fastqs)
    strict_metadata = str_to_bool(args.strict_metadata)
    strict_barcode_filenames = str_to_bool(args.strict_barcode_filenames)
    production_mode = str_to_bool(args.production_mode)
    extensions = [ext.strip() for ext in args.fastq_extensions.split(",") if ext.strip()]

    report = {
        "errors": [],
        "warnings": [],
        "summary": {},
        "missing_metadata_by_column": {},
        "fastq_inventory": {},
    }

    sample_columns, sample_rows = read_csv(args.samplesheet)
    primer_columns, primer_rows = read_csv(args.primers)

    missing_required = [col for col in REQUIRED_SAMPLE_COLUMNS if col not in sample_columns]
    if missing_required:
        report["errors"].append(f"Samplesheet missing required columns: {', '.join(missing_required)}")

    missing_recommended_cols = [col for col in RECOMMENDED_METADATA_COLUMNS if col not in sample_columns]
    if missing_recommended_cols:
        report["warnings"].append(
            "Samplesheet missing recommended MIEM/MIMARKS-style columns: "
            + ", ".join(missing_recommended_cols)
        )

    missing_production_cols = [col for col in PRODUCTION_METADATA_COLUMNS if col not in sample_columns]
    if missing_production_cols:
        message = (
            "Samplesheet missing production metadata columns: "
            + ", ".join(missing_production_cols)
        )
        if production_mode:
            report["errors"].append(message)
        else:
            report["warnings"].append(message)

    missing_primer_cols = [col for col in REQUIRED_PRIMER_COLUMNS if col not in primer_columns]
    if missing_primer_cols:
        report["errors"].append(f"Primer file missing required columns: {', '.join(missing_primer_cols)}")

    marker_windows = parse_marker_windows(args.marker_windows)
    marker_rows = []
    if not missing_primer_cols:
        for row in primer_rows:
            marker = (row.get("Gene") or "").strip()
            if not marker:
                report["errors"].append("Primer row has an empty Gene/marker name")
                continue
            try:
                expected_size = int(float(str(row.get("Size", "")).strip()))
            except ValueError:
                report["errors"].append(f"Primer row for {marker} has invalid Size: {row.get('Size')}")
                expected_size = 0
            if marker in marker_windows:
                min_len, max_len = marker_windows[marker]
            else:
                min_len = max(1, expected_size - args.marker_window_padding)
                max_len = expected_size + args.marker_window_padding
                report["warnings"].append(
                    f"No configured length window for {marker}; using Size ±{args.marker_window_padding} bp "
                    f"from the primer CSV ({min_len}-{max_len})."
                )
            marker_rows.append(
                {
                    "marker": marker,
                    "forward_primer": (row.get("Forward_Primer") or "").strip().replace(" ", ""),
                    "reverse_primer": (row.get("Reverse_Primer") or "").strip().replace(" ", ""),
                    "expected_size": expected_size,
                    "min_len": min_len,
                    "max_len": max_len,
                }
            )

    fastqs, fastq_errors = discover_fastqs(args.raw_data_dir, extensions)
    report["errors"].extend(fastq_errors)
    fastqs_by_run_barcode = defaultdict(list)
    barcode_filename_mismatches = []
    barcode_re = re.compile(r"barcode\d+")
    for fastq in fastqs:
        run_id, barcode_id = infer_run_barcode(fastq, args.raw_data_dir)
        if run_id and barcode_id:
            fastqs_by_run_barcode[(run_id, barcode_id)].append(str(fastq.resolve()))
            match = barcode_re.search(fastq.name)
            file_barcode = match.group(0) if match else ""
            if file_barcode and file_barcode != barcode_id:
                barcode_filename_mismatches.append(
                    {
                        "run_id": run_id,
                        "folder_barcode": barcode_id,
                        "file_barcode": file_barcode,
                        "fastq_path": str(fastq.resolve()),
                    }
                )

    row_keys = set()
    sample_ids = Counter()
    validated_rows = []
    sample_fastq_rows = []
    control_rows = []
    sample_fieldnames = list(sample_columns)
    for col in ["sample_uid"]:
        if col not in sample_fieldnames:
            sample_fieldnames.append(col)

    fastq_fieldnames = sample_fieldnames + ["n_fastq_files", "fastq_paths"]

    for line_no, row in enumerate(sample_rows, start=2):
        for col in REQUIRED_SAMPLE_COLUMNS:
            if col in sample_columns and not (row.get(col) or "").strip():
                report["errors"].append(f"Line {line_no}: required column '{col}' is empty")

        run_id = (row.get("run_id") or "").strip()
        barcode_id = (row.get("barcode_id") or "").strip()
        sample_id = (row.get("sample_id") or "").strip()
        sample_uid = sanitize(f"{run_id}__{barcode_id}__{sample_id}")
        row["sample_uid"] = sample_uid
        row_keys.add((run_id, barcode_id))
        sample_ids[sample_uid] += 1

        for col in RECOMMENDED_METADATA_COLUMNS:
            if col in sample_columns and (row.get("sample_type") or "").strip() == "sample":
                if not (row.get(col) or "").strip():
                    report["missing_metadata_by_column"][col] = report["missing_metadata_by_column"].get(col, 0) + 1

        sample_type = (row.get("sample_type") or "").strip()
        sample_type_lower = sample_type.lower()
        if production_mode:
            if sample_type_lower == "sample":
                for col in ["latitude", "longitude", "extraction_batch", "pcr_batch", "library_batch", "barcode_kit", "flowcell", "basecalling_model"]:
                    if col in sample_columns and not (row.get(col) or "").strip():
                        report["errors"].append(f"Line {line_no}: production sample row missing '{col}'")
            elif "control_type" in sample_columns and not (row.get("control_type") or "").strip():
                report["errors"].append(f"Line {line_no}: production control row missing 'control_type'")

            if "positive" in sample_type_lower:
                for col in ["expected_host_scientific_name", "expected_host_taxid", "expected_marker_result"]:
                    if col in sample_columns and not (row.get(col) or "").strip():
                        report["errors"].append(f"Line {line_no}: positive control row missing '{col}'")

            if args.barcode_kit and "barcode_kit" in sample_columns:
                row_kit = (row.get("barcode_kit") or "").strip()
                if row_kit and row_kit != args.barcode_kit:
                    report["errors"].append(
                        f"Line {line_no}: barcode_kit '{row_kit}' does not match --barcode-kit '{args.barcode_kit}'"
                    )

            if args.basecalling_model and "basecalling_model" in sample_columns:
                row_model = (row.get("basecalling_model") or "").strip()
                if row_model and row_model != args.basecalling_model:
                    report["errors"].append(
                        f"Line {line_no}: basecalling_model '{row_model}' does not match --basecalling-model '{args.basecalling_model}'"
                    )

            if "latitude" in sample_columns and "longitude" in sample_columns and sample_type_lower == "sample":
                try:
                    lat = float((row.get("latitude") or "").strip())
                    lon = float((row.get("longitude") or "").strip())
                    if not (-90 <= lat <= 90 and -180 <= lon <= 180):
                        raise ValueError
                except ValueError:
                    report["errors"].append(f"Line {line_no}: latitude/longitude are not valid decimal degrees")

        coord_col = "collection_cordinates"
        if coord_col in sample_columns and (row.get("sample_type") or "").strip() == "sample":
            coord = (row.get(coord_col) or "").strip()
            if coord and not validate_coordinate(coord):
                report["warnings"].append(
                    f"Line {line_no}: collection_cordinates '{coord}' is not parseable as 'lat,lon'."
                )

        found_fastqs = fastqs_for_sample(fastqs_by_run_barcode, run_id, barcode_id)
        if require_fastqs and not found_fastqs:
            report["errors"].append(
                f"Line {line_no}: no FASTQ files found for run_id={run_id}, barcode_id={barcode_id}"
            )
        row_for_fastq = dict(row)
        row_for_fastq["n_fastq_files"] = len(found_fastqs)
        row_for_fastq["fastq_paths"] = "|".join(found_fastqs)
        validated_rows.append(dict(row))
        sample_fastq_rows.append(row_for_fastq)

        if (row.get("sample_type") or "").strip() != "sample":
            control_rows.append(
                {
                    "sample_uid": sample_uid,
                    "run_id": run_id,
                    "barcode_id": barcode_id,
                    "sample_id": sample_id,
                    "sample_type": row.get("sample_type", ""),
                }
            )

    duplicated_uids = [uid for uid, count in sample_ids.items() if count > 1]
    if duplicated_uids:
        report["errors"].append("Duplicate sample_uid values after sanitisation: " + ", ".join(duplicated_uids))

    inventory_by_run = Counter()
    inventory_by_run_barcode = Counter()
    for (run_id, barcode_id), paths in fastqs_by_run_barcode.items():
        inventory_by_run[run_id] += len(paths)
        inventory_by_run_barcode[f"{run_id}/{barcode_id}"] = len(paths)

    extra_run_barcodes = sorted(set(fastqs_by_run_barcode) - row_keys)
    if extra_run_barcodes:
        preview = ", ".join(f"{run}/{bc}" for run, bc in extra_run_barcodes[:20])
        report["warnings"].append(
            f"FASTQ directories found without samplesheet rows ({len(extra_run_barcodes)}): {preview}"
        )

    if strict_metadata and report["missing_metadata_by_column"]:
        missing = ", ".join(f"{k}={v}" for k, v in sorted(report["missing_metadata_by_column"].items()))
        report["errors"].append(f"Strict metadata mode failed; missing recommended metadata values: {missing}")

    if barcode_filename_mismatches:
        preview = ", ".join(
            f"{row['run_id']}/{row['folder_barcode']} file={row['file_barcode']}"
            for row in barcode_filename_mismatches[:12]
        )
        report["warnings"].append(
            f"FASTQ filename barcode disagrees with containing folder for "
            f"{len(barcode_filename_mismatches)} files: {preview}"
        )
        if strict_barcode_filenames:
            report["errors"].append(
                "Strict barcode filename mode failed; see barcode_filename_mismatches.tsv."
            )

    report["summary"] = {
        "samplesheet_rows": len(sample_rows),
        "sample_rows": sum(1 for row in sample_rows if row.get("sample_type") == "sample"),
        "control_rows": sum(1 for row in sample_rows if row.get("sample_type") != "sample"),
        "markers": len(marker_rows),
        "fastq_files": len(fastqs),
        "run_count": len(inventory_by_run),
        "require_fastqs": require_fastqs,
        "strict_metadata": strict_metadata,
        "strict_barcode_filenames": strict_barcode_filenames,
        "production_mode": production_mode,
        "barcode_kit": args.barcode_kit,
        "basecalling_model": args.basecalling_model,
        "barcode_filename_mismatches": len(barcode_filename_mismatches),
    }
    report["fastq_inventory"] = {
        "by_run": dict(sorted(inventory_by_run.items())),
        "by_run_barcode": dict(sorted(inventory_by_run_barcode.items())),
    }

    write_tsv(outdir / "validated_samplesheet.tsv", validated_rows, sample_fieldnames)
    write_tsv(outdir / "sample_fastqs.tsv", sample_fastq_rows, fastq_fieldnames)
    write_tsv(
        outdir / "marker_config.tsv",
        marker_rows,
        ["marker", "forward_primer", "reverse_primer", "expected_size", "min_len", "max_len"],
    )
    write_tsv(
        outdir / "control_manifest.tsv",
        control_rows,
        ["sample_uid", "run_id", "barcode_id", "sample_id", "sample_type"],
    )
    write_tsv(
        outdir / "barcode_filename_mismatches.tsv",
        barcode_filename_mismatches,
        ["run_id", "folder_barcode", "file_barcode", "fastq_path"],
    )
    with (outdir / "input_validation_report.json").open("w") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")

    if report["errors"]:
        for message in report["errors"]:
            print(f"ERROR: {message}", file=sys.stderr)
        for message in report["warnings"][:20]:
            print(f"WARNING: {message}", file=sys.stderr)
        sys.exit(1)

    for message in report["warnings"][:20]:
        print(f"WARNING: {message}", file=sys.stderr)


if __name__ == "__main__":
    main()
