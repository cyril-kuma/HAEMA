#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path

from fastq_utils import read_fastq, write_fastq_record


BARCODE_PATTERNS = [
    re.compile(r"(?:barcode|barcode_id|bc)[=:_-]?(barcode\d{1,3})", re.IGNORECASE),
    re.compile(r"\b(barcode\d{1,3})\b", re.IGNORECASE),
]


def open_text(path, mode="wt"):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return path.open(mode)


def infer_barcode(header):
    for pattern in BARCODE_PATTERNS:
        match = pattern.search(header)
        if match:
            token = match.group(1).lower()
            number = int(re.search(r"\d+", token).group(0))
            return f"barcode{number:02d}"
    return "unclassified"


def read_fastq_list(path):
    with Path(path).open() as handle:
        return [line.strip() for line in handle if line.strip()]


def main():
    parser = argparse.ArgumentParser(
        description="Split pooled FASTQ by barcode tags in read headers into a MinKNOW-like fastq_pass layout"
    )
    parser.add_argument("--fastq-list", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--min-reads-per-barcode", type=int, default=1)
    args = parser.parse_args()

    output_root = Path(args.output_root)
    fastq_pass = output_root / args.run_id / "advanced_demux" / "fastq_pass"
    fastq_pass.mkdir(parents=True, exist_ok=True)

    writers = {}
    counts = Counter()
    source_counts = defaultdict(Counter)
    try:
        for fastq in read_fastq_list(args.fastq_list):
            source_name = Path(fastq).name
            for header, seq, plus, qual in read_fastq(fastq):
                barcode = infer_barcode(header)
                counts[barcode] += 1
                source_counts[source_name][barcode] += 1
                if barcode == "unclassified":
                    continue
                if barcode not in writers:
                    barcode_dir = fastq_pass / barcode
                    barcode_dir.mkdir(parents=True, exist_ok=True)
                    writers[barcode] = gzip.open(barcode_dir / f"{args.run_id}_{barcode}.advanced_demux.fastq.gz", "wt")
                write_fastq_record(writers[barcode], header, seq, qual)
    finally:
        for writer in writers.values():
            writer.close()

    # Remove low-support barcode files after writing. This avoids publishing
    # accidental one-read barcode bins unless the user explicitly allows them.
    retained = set()
    for barcode, count in counts.items():
        if barcode != "unclassified" and count >= args.min_reads_per_barcode:
            retained.add(barcode)
        elif barcode != "unclassified":
            target = fastq_pass / barcode
            for file_path in target.glob("*.fastq.gz"):
                file_path.unlink()
            try:
                target.rmdir()
            except OSError:
                pass

    with Path(args.summary).open("w", newline="") as handle:
        fieldnames = ["run_id", "source_file", "barcode_id", "reads", "retained", "method"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for source_file, per_barcode in sorted(source_counts.items()):
            for barcode, count in sorted(per_barcode.items()):
                writer.writerow(
                    {
                        "run_id": args.run_id,
                        "source_file": source_file,
                        "barcode_id": barcode,
                        "reads": count,
                        "retained": str(barcode in retained).lower(),
                        "method": "header_tag_split",
                    }
                )


if __name__ == "__main__":
    main()
