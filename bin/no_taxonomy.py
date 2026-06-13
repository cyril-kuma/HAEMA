#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Emit explicit unassigned taxonomy rows when taxonomy is skipped")
    parser.add_argument("--counts", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with Path(args.counts).open(newline="") as inp, Path(args.output).open("w", newline="") as out:
        reader = csv.DictReader(inp, delimiter="\t")
        fieldnames = list(reader.fieldnames or []) + [
            "host_assignment",
            "taxon_rank",
            "assignment_status",
            "confidence",
            "n_top_hits",
            "top_sseqid",
            "top_stitle",
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
        ]
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            if row.get("retained") != "true":
                continue
            row.update(
                {
                    "host_assignment": "unassigned",
                    "taxon_rank": "none",
                    "assignment_status": "taxonomy_skipped",
                    "confidence": "not_evaluated",
                    "n_top_hits": 0,
                    "top_sseqid": "",
                    "top_stitle": "",
                    "pident": "",
                    "coverage": "",
                    "evalue": "",
                    "bitscore": "",
                    "blast_source": "",
                    "fallback_used": "false",
                    "primary_assignment_status": "",
                    "assignment_method": "taxonomy_disabled",
                    "lca_taxid": "",
                    "top_staxids": "",
                }
            )
            writer.writerow(row)


if __name__ == "__main__":
    main()
