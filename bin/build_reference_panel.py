#!/usr/bin/env python3
import argparse
import csv
import hashlib
from pathlib import Path


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def read_fasta(path):
    records = {}
    header = None
    parts = []
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(parts)
                header = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if header is not None:
            records[header] = "".join(parts)
    return records


def md5_sequence(seq):
    return hashlib.md5(seq.upper().encode()).hexdigest()


def accession_from_seqid(seqid):
    parts = seqid.split("|")
    return parts[1] if len(parts) > 1 else seqid


def write_fasta(records, path):
    with Path(path).open("w") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def write_tsv(path, rows, fieldnames):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description="Build a production-ready curated reference panel from accession/taxid/checksum-complete rows"
    )
    parser.add_argument("--source-fasta", required=True)
    parser.add_argument("--taxonomy-sidecar", required=True)
    parser.add_argument("--target-manifest", required=True)
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--output-taxonomy", required=True)
    parser.add_argument("--needs-accession-report", required=True)
    args = parser.parse_args()

    fasta = read_fasta(args.source_fasta)
    _tax_cols, taxonomy_rows = read_tsv(args.taxonomy_sidecar)
    _target_cols, targets = read_tsv(args.target_manifest)

    target_by_accession = {
        (row.get("accession") or "").strip(): row
        for row in targets
        if (row.get("accession") or "").strip()
    }
    taxonomy_by_accession = {
        (row.get("source_accession") or accession_from_seqid(row.get("seqid", ""))).strip(): row
        for row in taxonomy_rows
    }

    ready_records = []
    ready_taxonomy = []
    report_rows = []
    for target in targets:
        accession = (target.get("accession") or "").strip()
        taxid = (target.get("ncbi_taxid") or "").strip()
        status = (target.get("status") or "").strip().lower()
        base_report = {
            "priority": target.get("priority", ""),
            "group": target.get("group", ""),
            "scientific_name": target.get("scientific_name", ""),
            "common_name": target.get("common_name", ""),
            "ncbi_taxid": taxid,
            "accession": accession,
            "status": status or "missing_status",
            "included": "false",
            "reason": "",
        }
        if status != "ready":
            base_report["reason"] = status or "not_ready"
            report_rows.append(base_report)
            continue
        if not accession or not taxid:
            base_report["reason"] = "ready_row_missing_accession_or_taxid"
            report_rows.append(base_report)
            continue
        tax = taxonomy_by_accession.get(accession)
        if not tax:
            base_report["reason"] = "accession_missing_from_taxonomy_sidecar"
            report_rows.append(base_report)
            continue
        seqid = tax.get("seqid", "")
        seq = fasta.get(seqid)
        if not seq:
            base_report["reason"] = "seqid_missing_from_source_fasta"
            report_rows.append(base_report)
            continue
        expected_md5 = (tax.get("sequence_md5") or "").strip()
        observed_md5 = md5_sequence(seq)
        if not expected_md5 or expected_md5 != observed_md5:
            base_report["reason"] = "missing_or_mismatched_sequence_md5"
            base_report["observed_md5"] = observed_md5
            report_rows.append(base_report)
            continue
        ready_records.append((seqid, seq))
        ready_taxonomy.append(tax)
        base_report["included"] = "true"
        base_report["reason"] = "ready"
        base_report["observed_md5"] = observed_md5
        report_rows.append(base_report)

    write_fasta(ready_records, args.output_fasta)
    taxonomy_fields = list(taxonomy_rows[0].keys()) if taxonomy_rows else [
        "seqid", "taxid", "scientific_name", "rank", "common_name", "source_accession",
        "source_db", "retrieval_date", "provenance", "sequence_md5",
    ]
    write_tsv(args.output_taxonomy, ready_taxonomy, taxonomy_fields)
    report_fields = [
        "priority", "group", "scientific_name", "common_name", "ncbi_taxid", "accession",
        "status", "included", "reason", "observed_md5",
    ]
    write_tsv(args.needs_accession_report, report_rows, report_fields)


if __name__ == "__main__":
    main()
