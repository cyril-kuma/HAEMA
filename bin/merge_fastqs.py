#!/usr/bin/env python3
import argparse
import csv
import gzip
from pathlib import Path

from fastq_utils import mean_q, read_fastq, write_fastq_record


def main():
    parser = argparse.ArgumentParser(description="Merge FASTQ chunks and collect basic read stats")
    parser.add_argument("--sample-uid", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--barcode-id", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--stats", required=True)
    parser.add_argument("fastqs", nargs="+")
    args = parser.parse_args()

    reads = 0
    bases = 0
    q_sum = 0.0
    min_len = None
    max_len = 0
    with gzip.open(args.output, "wt") as out:
        for fastq in args.fastqs:
            for header, seq, _plus, qual in read_fastq(fastq):
                reads += 1
                length = len(seq)
                bases += length
                q_sum += mean_q(qual)
                min_len = length if min_len is None else min(min_len, length)
                max_len = max(max_len, length)
                write_fastq_record(out, header, seq, qual)

    mean_len = bases / reads if reads else 0.0
    mean_quality = q_sum / reads if reads else 0.0
    with Path(args.stats).open("w", newline="") as handle:
        fieldnames = [
            "sample_uid",
            "run_id",
            "sample_id",
            "barcode_id",
            "n_files",
            "reads",
            "bases",
            "min_len",
            "max_len",
            "mean_len",
            "mean_read_q",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "sample_uid": args.sample_uid,
                "run_id": args.run_id,
                "sample_id": args.sample_id,
                "barcode_id": args.barcode_id,
                "n_files": len(args.fastqs),
                "reads": reads,
                "bases": bases,
                "min_len": min_len or 0,
                "max_len": max_len,
                "mean_len": f"{mean_len:.2f}",
                # Unweighted mean of per-read mean Phred scores (not base-weighted); named explicitly.
                "mean_read_q": f"{mean_quality:.2f}",
            }
        )


if __name__ == "__main__":
    main()
