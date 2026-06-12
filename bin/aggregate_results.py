#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def read_tsv(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def feature_id(marker, sequence):
    digest = hashlib.md5(sequence.encode()).hexdigest()[:12]
    return f"{marker}_{digest}"


def main():
    parser = argparse.ArgumentParser(description="Aggregate HÆMA blood-meal endpoint tables")
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--counts", nargs="*", default=[])
    parser.add_argument("--assignments", nargs="*", default=[])
    parser.add_argument("--preprocess-qc", nargs="*", default=[])
    parser.add_argument("--derep-summaries", nargs="*", default=[])
    parser.add_argument("--negative-control-multiplier", type=float, default=1.0)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    metadata_rows = read_tsv(args.samplesheet)
    metadata_by_uid = {row.get("sample_uid"): row for row in metadata_rows}

    count_rows = []
    for path in args.counts:
        count_rows.extend(read_tsv(path))

    assignment_rows = []
    for path in args.assignments:
        assignment_rows.extend(read_tsv(path))
    assignments_by_asv = {row.get("asv_id"): row for row in assignment_rows}

    master_rows = []
    negative_sequences = defaultdict(int)
    negative_hosts = defaultdict(int)
    for row in assignment_rows:
        meta = metadata_by_uid.get(row.get("sample_uid"), {})
        sample_type = meta.get("sample_type", row.get("sample_type", ""))
        try:
            count = int(float(row.get("count") or 0))
        except ValueError:
            count = 0
        if "negative" in sample_type.lower():
            negative_sequences[(row.get("marker"), row.get("sequence"))] = max(
                negative_sequences[(row.get("marker"), row.get("sequence"))],
                count,
            )
            host = row.get("host_assignment") or ""
            if host and host != "unassigned":
                negative_hosts[(row.get("marker"), host)] = max(negative_hosts[(row.get("marker"), host)], count)

    contamination_rows = []
    for row in assignment_rows:
        meta = metadata_by_uid.get(row.get("sample_uid"), {})
        combined = dict(meta)
        combined.update(row)
        combined["control_status"] = meta.get("sample_type", row.get("sample_type", "unknown"))
        try:
            count = int(float(row.get("count") or 0))
        except ValueError:
            count = 0
        marker = row.get("marker")
        seq_key = (marker, row.get("sequence"))
        host_key = (marker, row.get("host_assignment"))
        max_neg_seq = negative_sequences.get(seq_key, 0)
        max_neg_host = negative_hosts.get(host_key, 0)
        flag = False
        reasons = []
        if "negative" not in combined["control_status"].lower() and max_neg_seq:
            if count <= max_neg_seq * args.negative_control_multiplier:
                flag = True
                reasons.append(f"sequence_seen_in_negative_control_max_count={max_neg_seq}")
        if "negative" not in combined["control_status"].lower() and max_neg_host and row.get("host_assignment") != "unassigned":
            if count <= max_neg_host * args.negative_control_multiplier:
                flag = True
                reasons.append(f"host_seen_in_negative_control_max_count={max_neg_host}")
        combined["present_in_negative_control"] = str(bool(max_neg_seq or max_neg_host)).lower()
        combined["max_negative_sequence_count"] = max_neg_seq
        combined["max_negative_host_count"] = max_neg_host
        combined["contamination_flag"] = str(flag).lower()
        combined["contamination_reason"] = ";".join(reasons)
        master_rows.append(combined)
        if flag:
            contamination_rows.append(
                {
                    "sample_uid": row.get("sample_uid"),
                    "run_id": row.get("run_id"),
                    "barcode_id": row.get("barcode_id"),
                    "sample_id": row.get("sample_id"),
                    "marker": marker,
                    "cluster_id": row.get("cluster_id", ""),
                    "asv_id": row.get("asv_id"),
                    "host_assignment": row.get("host_assignment"),
                    "count": count,
                    "contamination_flag": "true",
                    "reason": ";".join(reasons),
                }
            )

    master_fieldnames = [
        "sample_uid",
        "run_id",
        "sample_id",
        "barcode_id",
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
        "bioclimatic_zone",
        "collection_region",
        "collection_cordinates",
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
        "marker",
        "cluster_id",
        "asv_id",
        "sequence",
        "count",
        "fraction",
        "host_assignment",
        "taxon_rank",
        "assignment_status",
        "confidence",
        "pident",
        "coverage",
        "evalue",
        "bitscore",
        "blast_source",
        "fallback_used",
        "primary_assignment_status",
        "assignment_method",
        "lca_taxid",
        "top_staxids",
        "control_status",
        "present_in_negative_control",
        "max_negative_sequence_count",
        "max_negative_host_count",
        "contamination_flag",
        "contamination_reason",
    ]
    write_tsv(outdir / "bloodmeal_master_endpoint.tsv", master_rows, master_fieldnames)
    write_tsv(outdir / "host_assignments.tsv", master_rows, master_fieldnames)

    features = {}
    sample_uids = [row.get("sample_uid") for row in metadata_rows if row.get("sample_uid")]
    for row in count_rows:
        if row.get("retained") != "true":
            continue
        fid = feature_id(row.get("marker", ""), row.get("sequence", ""))
        features.setdefault(
            fid,
            {
                "feature_id": fid,
                "marker": row.get("marker"),
                "sequence": row.get("sequence"),
                **{uid: 0 for uid in sample_uids},
            },
        )
        try:
            features[fid][row.get("sample_uid")] += int(float(row.get("count") or 0))
        except ValueError:
            pass
    write_tsv(outdir / "asv_count_table.tsv", list(features.values()), ["feature_id", "marker", "sequence"] + sample_uids)

    by_sample = defaultdict(lambda: Counter())
    by_marker = defaultdict(lambda: Counter())
    for row in master_rows:
        uid = row.get("sample_uid")
        marker = row.get("marker")
        try:
            count = int(float(row.get("count") or 0))
        except ValueError:
            count = 0
        by_sample[uid]["total_asv_reads"] += count
        by_sample[uid]["retained_asvs"] += 1
        if row.get("host_assignment") not in {"", "unassigned", "ambiguous"}:
            by_sample[uid]["assigned_hosts"] += 1
        if row.get("contamination_flag") == "true":
            by_sample[uid]["contamination_flags"] += 1
        key = (uid, marker)
        by_marker[key]["total_asv_reads"] += count
        by_marker[key]["retained_asvs"] += 1
        if row.get("host_assignment") not in {"", "unassigned", "ambiguous"}:
            by_marker[key]["assigned_hosts"] += 1

    sample_summary = []
    for uid, counts in sorted(by_sample.items()):
        meta = metadata_by_uid.get(uid, {})
        sample_summary.append(
            {
                "sample_uid": uid,
                "run_id": meta.get("run_id", ""),
                "sample_id": meta.get("sample_id", ""),
                "barcode_id": meta.get("barcode_id", ""),
                "sample_type": meta.get("sample_type", ""),
                "total_asv_reads": counts["total_asv_reads"],
                "retained_asvs": counts["retained_asvs"],
                "assigned_hosts": counts["assigned_hosts"],
                "contamination_flags": counts["contamination_flags"],
            }
        )
    write_tsv(
        outdir / "sample_level_summary.tsv",
        sample_summary,
        [
            "sample_uid",
            "run_id",
            "sample_id",
            "barcode_id",
            "sample_type",
            "total_asv_reads",
            "retained_asvs",
            "assigned_hosts",
            "contamination_flags",
        ],
    )

    derep_by_key = {}
    for path in args.derep_summaries:
        for row in read_tsv(path):
            derep_by_key[(row.get("sample_uid"), row.get("marker"))] = row
    marker_summary = []
    for (uid, marker), counts in sorted(by_marker.items()):
        derep = derep_by_key.get((uid, marker), {})
        marker_summary.append(
            {
                "sample_uid": uid,
                "marker": marker,
                "total_asv_reads": counts["total_asv_reads"],
                "retained_asvs": counts["retained_asvs"],
                "assigned_hosts": counts["assigned_hosts"],
                "mixed_template_warning": derep.get("mixed_template_warning", "false"),
            }
        )
    write_tsv(
        outdir / "marker_level_summary.tsv",
        marker_summary,
        ["sample_uid", "marker", "total_asv_reads", "retained_asvs", "assigned_hosts", "mixed_template_warning"],
    )

    qc_rows = []
    for path in args.preprocess_qc:
        for row in read_tsv(path):
            row = dict(row)
            row["source_file"] = Path(path).name
            qc_rows.append(row)
    qc_fields = sorted({key for row in qc_rows for key in row}) or ["source_file"]
    write_tsv(outdir / "qc_summary.tsv", qc_rows, qc_fields)

    write_tsv(
        outdir / "contamination_flags.tsv",
        contamination_rows,
        [
            "sample_uid",
            "run_id",
            "barcode_id",
            "sample_id",
            "marker",
            "cluster_id",
            "asv_id",
            "host_assignment",
            "count",
            "contamination_flag",
            "reason",
        ],
    )

    manifest = {
        "primary_endpoint": "bloodmeal_master_endpoint.tsv",
        "records": len(master_rows),
        "features": len(features),
        "samples": len(sample_uids),
        "contamination_flags": len(contamination_rows),
        "notes": [
            "Primary endpoint preserves controls and flags contamination; it does not silently remove rows.",
            "Mixed-host evidence is preserved at feature level; optional RAMBO-style host calls are reported separately.",
            "Optional Medaka polishing, phyloseq/decontam outputs, and taxid-backed LCA are recorded in the run manifest when enabled.",
        ],
    }
    with (outdir / "endpoint_manifest.json").open("w") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
