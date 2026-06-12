#!/usr/bin/env python3
import argparse
import csv
from collections import Counter
from pathlib import Path

from fastq_utils import read_fastq


def main():
    parser = argparse.ArgumentParser(description="Dereplicate marker reads into first-pass ASV-like features")
    parser.add_argument("--sample-uid", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--barcode-id", required=True)
    parser.add_argument("--marker", required=True)
    parser.add_argument("--cluster-id", default="all")
    parser.add_argument("--input", required=True)
    parser.add_argument("--min-reads", type=int, default=3)
    parser.add_argument("--min-fraction", type=float, default=0.01)
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--output-counts", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    counts = Counter()
    total = 0
    for _header, seq, _plus, _qual in read_fastq(args.input):
        seq = seq.upper()
        if not seq:
            continue
        counts[seq] += 1
        total += 1

    retained = []
    for rank, (sequence, count) in enumerate(counts.most_common(), start=1):
        fraction = count / total if total else 0.0
        keep = count >= args.min_reads and fraction >= args.min_fraction
        asv_id = f"{args.sample_uid}|{args.marker}|{args.cluster_id}|ASV{rank:04d}"
        retained.append((asv_id, sequence, count, fraction, keep))

    with Path(args.output_fasta).open("w") as fasta:
        for asv_id, sequence, _count, _fraction, keep in retained:
            if keep:
                fasta.write(f">{asv_id}\n{sequence}\n")

    with Path(args.output_counts).open("w", newline="") as handle:
        fieldnames = [
            "sample_uid",
            "run_id",
            "sample_id",
            "barcode_id",
            "marker",
            "cluster_id",
            "asv_id",
            "sequence",
            "count",
            "fraction",
            "retained",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for asv_id, sequence, count, fraction, keep in retained:
            writer.writerow(
                {
                    "sample_uid": args.sample_uid,
                    "run_id": args.run_id,
                    "sample_id": args.sample_id,
                    "barcode_id": args.barcode_id,
                    "marker": args.marker,
                    "cluster_id": args.cluster_id,
                    "asv_id": asv_id,
                    "sequence": sequence,
                    "count": count,
                    "fraction": f"{fraction:.6f}",
                    "retained": str(keep).lower(),
                }
            )

    n_retained = sum(1 for item in retained if item[-1])
    with Path(args.summary).open("w", newline="") as handle:
        fieldnames = [
            "sample_uid",
            "run_id",
            "sample_id",
            "barcode_id",
            "marker",
            "cluster_id",
            "n_reads",
            "n_unique",
            "n_retained",
            "mixed_template_warning",
            "method",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "sample_uid": args.sample_uid,
                "run_id": args.run_id,
                "sample_id": args.sample_id,
                "barcode_id": args.barcode_id,
                "marker": args.marker,
                "cluster_id": args.cluster_id,
                "n_reads": total,
                "n_unique": len(counts),
                "n_retained": n_retained,
                "mixed_template_warning": str(n_retained > 1).lower(),
                "method": "dereplicate_no_consensus_collapse",
            }
        )


if __name__ == "__main__":
    main()
