#!/usr/bin/env python3
import argparse
import csv
from collections import Counter
from pathlib import Path

from fastq_utils import read_fastq


def approx_identity(a, b):
    if not a or not b:
        return 0.0
    max_len = max(len(a), len(b))
    min_len = min(len(a), len(b))
    mismatches = abs(len(a) - len(b))
    mismatches += sum(1 for i in range(min_len) if a[i] != b[i])
    return 1.0 - (mismatches / max_len)


def consensus(sequences):
    if not sequences:
        return ""
    max_len = max(len(seq) for seq in sequences)
    bases = []
    for i in range(max_len):
        counts = Counter(seq[i] for seq in sequences if i < len(seq))
        if not counts:
            continue
        base, _count = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0]
        bases.append(base)
    return "".join(bases)


def main():
    parser = argparse.ArgumentParser(
        description="Greedy read clustering and consensus generation for ONT amplicon features"
    )
    parser.add_argument("--sample-uid", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--barcode-id", required=True)
    parser.add_argument("--marker", required=True)
    parser.add_argument("--cluster-id", default="all")
    parser.add_argument("--input", required=True)
    parser.add_argument("--min-identity", type=float, default=0.97)
    parser.add_argument("--min-reads", type=int, default=3)
    parser.add_argument("--min-fraction", type=float, default=0.01)
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--output-counts", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    raw_sequences = []
    for _header, seq, _plus, _qual in read_fastq(args.input):
        seq = seq.upper()
        if seq:
            raw_sequences.append(seq)

    clusters = []
    for seq in raw_sequences:
        best_idx = None
        best_identity = -1.0
        for idx, cluster in enumerate(clusters):
            identity = approx_identity(seq, cluster["representative"])
            if identity > best_identity:
                best_idx = idx
                best_identity = identity
        if best_idx is not None and best_identity >= args.min_identity:
            clusters[best_idx]["sequences"].append(seq)
            clusters[best_idx]["representative"] = consensus(clusters[best_idx]["sequences"])
        else:
            clusters.append({"representative": seq, "sequences": [seq]})

    clusters.sort(key=lambda cluster: len(cluster["sequences"]), reverse=True)
    total = len(raw_sequences)

    retained = []
    for rank, cluster in enumerate(clusters, start=1):
        count = len(cluster["sequences"])
        fraction = count / total if total else 0.0
        keep = count >= args.min_reads and fraction >= args.min_fraction
        feature_id = f"{args.sample_uid}|{args.marker}|{args.cluster_id}|CONS{rank:04d}"
        feature_sequence = consensus(cluster["sequences"])
        retained.append((feature_id, feature_sequence, count, fraction, keep))

    with Path(args.output_fasta).open("w") as fasta:
        for feature_id, sequence, _count, _fraction, keep in retained:
            if keep:
                fasta.write(f">{feature_id}\n{sequence}\n")

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
        for feature_id, sequence, count, fraction, keep in retained:
            writer.writerow(
                {
                    "sample_uid": args.sample_uid,
                    "run_id": args.run_id,
                    "sample_id": args.sample_id,
                    "barcode_id": args.barcode_id,
                    "marker": args.marker,
                    "cluster_id": args.cluster_id,
                    "asv_id": feature_id,
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
            "min_cluster_identity",
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
                "n_unique": len(set(raw_sequences)),
                "n_retained": n_retained,
                "mixed_template_warning": str(n_retained > 1).lower(),
                "method": "greedy_cluster_consensus_no_medaka",
                "min_cluster_identity": args.min_identity,
            }
        )


if __name__ == "__main__":
    main()
