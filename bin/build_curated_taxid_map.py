#!/usr/bin/env python3
import argparse
import csv
import json
import sys
from pathlib import Path


def detect_dialect(path):
    text = Path(path).read_text().splitlines()
    sample = "\n".join(text[:5])
    if "\t" in sample:
        return "\t"
    return ","


def read_table(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return []
    delimiter = detect_dialect(path)
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def fasta_ids(path):
    ids = []
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def normalize_row(row):
    seqid = row.get("seqid") or row.get("sequence_id") or row.get("accession") or row.get("id") or ""
    taxid = row.get("taxid") or row.get("ncbi_taxid") or row.get("taxonomy_id") or ""
    return {
        "seqid": seqid.strip(),
        "taxid": str(taxid).strip(),
        "scientific_name": (row.get("scientific_name") or row.get("species") or row.get("name") or "").strip(),
        "rank": (row.get("rank") or "species").strip(),
        "common_name": (row.get("common_name") or "").strip(),
        "source_accession": (row.get("source_accession") or row.get("accession") or "").strip(),
        "provenance": (row.get("provenance") or row.get("source") or "").strip(),
    }


def write_tsv(path, rows, fieldnames):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description="Validate curated reference taxonomy metadata and build makeblastdb taxid map")
    parser.add_argument("--reference-fasta", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--require-taxids", default="false")
    parser.add_argument("--taxid-map", required=True)
    parser.add_argument("--taxonomy-table", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    require_taxids = str(args.require_taxids).strip().lower() in {"1", "true", "yes", "y", "on"}
    reference_ids = fasta_ids(args.reference_fasta)
    reference_id_set = set(reference_ids)
    rows = [normalize_row(row) for row in read_table(args.metadata)]

    errors = []
    warnings = []
    rows_by_id = {}
    for row in rows:
        if not row["seqid"]:
            warnings.append("Skipping taxonomy metadata row without seqid.")
            continue
        if row["seqid"] not in reference_id_set:
            warnings.append(f"Taxonomy metadata seqid not found in FASTA: {row['seqid']}")
            continue
        if not row["taxid"].isdigit() or row["taxid"] == "0":
            warnings.append(f"Taxonomy metadata row for {row['seqid']} lacks a valid positive integer taxid.")
            continue
        rows_by_id[row["seqid"]] = row

    missing = [seqid for seqid in reference_ids if seqid not in rows_by_id]
    if missing:
        message = f"{len(missing)} FASTA records lack usable curated taxid metadata."
        if require_taxids:
            errors.append(message + " Provide --curated_reference_metadata or disable --require_curated_taxids.")
        else:
            warnings.append(message)

    usable_rows = [rows_by_id[seqid] for seqid in reference_ids if seqid in rows_by_id]
    with Path(args.taxid_map).open("w") as handle:
        for row in usable_rows:
            handle.write(f"{row['seqid']}\t{row['taxid']}\n")

    write_tsv(
        args.taxonomy_table,
        usable_rows,
        ["seqid", "taxid", "scientific_name", "rank", "common_name", "source_accession", "provenance"],
    )

    report = {
        "reference_records": len(reference_ids),
        "metadata_records": len(rows),
        "usable_taxid_records": len(usable_rows),
        "missing_taxid_records": len(missing),
        "require_taxids": require_taxids,
        "errors": errors,
        "warnings": warnings,
    }
    with Path(args.report).open("w") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")

    for message in warnings[:20]:
        print(f"WARNING: {message}", file=sys.stderr)
    if errors:
        for message in errors:
            print(f"ERROR: {message}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
