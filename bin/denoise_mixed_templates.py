#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
from collections import Counter, defaultdict
from pathlib import Path

from fastq_utils import read_fastq, write_fastq_record


def read_records(path):
    records = []
    for header, seq, _plus, qual in read_fastq(path):
        if not seq:
            continue
        read_id = header.lstrip("@").split()[0]
        records.append({"read_id": read_id, "header": header, "sequence": seq.upper(), "quality": qual})
    return records


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
        if counts:
            bases.append(sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0])
    return "".join(bases)


def greedy_labels(records, min_identity):
    clusters = []
    labels = []
    for record in records:
        seq = record["sequence"]
        best_idx = None
        best_identity = -1.0
        for idx, cluster in enumerate(clusters):
            identity = approx_identity(seq, cluster["representative"])
            if identity > best_identity:
                best_idx = idx
                best_identity = identity
        if best_idx is not None and best_identity >= min_identity:
            clusters[best_idx]["sequences"].append(seq)
            clusters[best_idx]["representative"] = consensus(clusters[best_idx]["sequences"])
            labels.append(best_idx)
        else:
            clusters.append({"representative": seq, "sequences": [seq]})
            labels.append(len(clusters) - 1)
    return labels, "greedy_identity"


def kmer_vector(sequence, k, vocabulary):
    counts = Counter(sequence[i : i + k] for i in range(max(0, len(sequence) - k + 1)))
    total = sum(counts.values()) or 1
    return [counts.get(kmer, 0) / total for kmer in vocabulary]


def umap_hdbscan_labels(records, k, min_cluster_size):
    try:
        import hdbscan
        import umap
    except Exception as exc:
        raise RuntimeError(f"UMAP/HDBSCAN dependencies unavailable: {exc}") from exc

    vocabulary = sorted(
        {
            seq[i : i + k]
            for seq in (record["sequence"] for record in records)
            for i in range(max(0, len(seq) - k + 1))
            if set(seq[i : i + k]) <= set("ACGT")
        }
    )
    if not vocabulary:
        raise RuntimeError("No valid A/C/G/T k-mers were available for UMAP/HDBSCAN clustering")

    matrix = [kmer_vector(record["sequence"], k, vocabulary) for record in records]
    n_neighbors = min(15, max(2, len(records) - 1))
    embedding = umap.UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        min_dist=0.0,
        metric="cosine",
        random_state=42,
    ).fit_transform(matrix)
    clusterer = hdbscan.HDBSCAN(min_cluster_size=max(2, min_cluster_size), min_samples=1)
    labels = list(clusterer.fit_predict(embedding))
    if not any(label >= 0 for label in labels):
        raise RuntimeError("HDBSCAN assigned all reads to noise")
    return labels, "umap_hdbscan"


def retained_cluster_map(labels, total, min_cluster_size, min_cluster_fraction):
    counts = Counter(label for label in labels if label >= 0)
    retained = []
    for label, count in counts.most_common():
        fraction = count / total if total else 0.0
        if count >= min_cluster_size and fraction >= min_cluster_fraction:
            retained.append(label)
    if not retained and counts:
        retained = [counts.most_common(1)[0][0]]
    return {label: f"cluster{rank:03d}" for rank, label in enumerate(retained, start=1)}


def main():
    parser = argparse.ArgumentParser(
        description="Pre-consensus mixed-template denoising for ONT blood-meal marker reads"
    )
    parser.add_argument("--sample-uid", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--barcode-id", required=True)
    parser.add_argument("--marker", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--backend", default="umap_hdbscan")
    parser.add_argument("--kmer-size", type=int, default=5)
    parser.add_argument("--min-cluster-size", type=int, default=3)
    parser.add_argument("--min-cluster-fraction", type=float, default=0.01)
    parser.add_argument("--min-reads-for-umap", type=int, default=50)
    parser.add_argument("--greedy-min-identity", type=float, default=0.97)
    parser.add_argument("--allow-greedy-fallback", default="true")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--membership", required=True)
    args = parser.parse_args()

    allow_fallback = str(args.allow_greedy_fallback).lower() in {"1", "true", "yes", "y", "on"}
    records = read_records(args.input)
    total = len(records)
    backend_requested = args.backend
    fallback_reason = ""

    if total == 0:
        labels = []
        backend_used = "empty_input"
    elif backend_requested == "none":
        labels = [0 for _ in records]
        backend_used = "none_single_cluster"
    elif backend_requested == "greedy":
        labels, backend_used = greedy_labels(records, args.greedy_min_identity)
    elif backend_requested == "umap_hdbscan":
        if total < args.min_reads_for_umap:
            if not allow_fallback:
                raise SystemExit(
                    f"Only {total} reads are available, below --min-reads-for-umap {args.min_reads_for_umap}."
                )
            labels, backend_used = greedy_labels(records, args.greedy_min_identity)
            fallback_reason = "below_min_reads_for_umap"
        else:
            try:
                labels, backend_used = umap_hdbscan_labels(records, args.kmer_size, args.min_cluster_size)
            except Exception as exc:
                if not allow_fallback:
                    raise
                labels, backend_used = greedy_labels(records, args.greedy_min_identity)
                fallback_reason = str(exc)
    else:
        raise SystemExit(f"Unsupported --backend '{backend_requested}'")

    label_to_cluster = retained_cluster_map(labels, total, args.min_cluster_size, args.min_cluster_fraction)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cluster_handles = {}
    cluster_counts = defaultdict(int)
    membership_rows = []
    for record, label in zip(records, labels):
        cluster_id = label_to_cluster.get(label, "noise")
        retained = cluster_id != "noise"
        reason = "retained" if retained else "below_cluster_threshold_or_hdbscan_noise"
        if retained:
            if cluster_id not in cluster_handles:
                cluster_handles[cluster_id] = gzip.open(
                    output_dir / f"{args.sample_uid}.{args.marker}.{cluster_id}.fastq.gz", "wt"
                )
            write_fastq_record(cluster_handles[cluster_id], record["header"], record["sequence"], record["quality"])
            cluster_counts[cluster_id] += 1
        membership_rows.append(
            {
                "sample_uid": args.sample_uid,
                "run_id": args.run_id,
                "sample_id": args.sample_id,
                "barcode_id": args.barcode_id,
                "marker": args.marker,
                "read_id": record["read_id"],
                "raw_cluster_label": label,
                "cluster_id": cluster_id,
                "retained": str(retained).lower(),
                "reason": reason,
                "read_length": len(record["sequence"]),
                "backend_requested": backend_requested,
                "backend_used": backend_used,
            }
        )
    for handle in cluster_handles.values():
        handle.close()
    empty_placeholder_cluster = total == 0 or not cluster_counts
    if empty_placeholder_cluster:
        empty_path = output_dir / f"{args.sample_uid}.{args.marker}.cluster001.fastq.gz"
        with gzip.open(empty_path, "wt"):
            pass
        cluster_counts["cluster001"] = 0

    noise_reads = sum(1 for row in membership_rows if row["cluster_id"] == "noise")
    retained_reads = total - noise_reads
    n_clusters = sum(1 for value in cluster_counts.values() if value > 0)
    mixed_warning = n_clusters > 1

    with Path(args.membership).open("w", newline="") as handle:
        fieldnames = [
            "sample_uid",
            "run_id",
            "sample_id",
            "barcode_id",
            "marker",
            "read_id",
            "raw_cluster_label",
            "cluster_id",
            "retained",
            "reason",
            "read_length",
            "backend_requested",
            "backend_used",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(membership_rows)

    with Path(args.summary).open("w", newline="") as handle:
        fieldnames = [
            "sample_uid",
            "run_id",
            "sample_id",
            "barcode_id",
            "marker",
            "input_reads",
            "retained_reads",
            "noise_reads",
            "n_clusters",
            "cluster_read_counts",
            "backend_requested",
            "backend_used",
            "fallback_used",
            "fallback_reason",
            "mixed_template_warning",
            "min_cluster_size",
            "min_cluster_fraction",
            "kmer_size",
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
                "input_reads": total,
                "retained_reads": retained_reads,
                "noise_reads": noise_reads,
                "n_clusters": n_clusters,
                "cluster_read_counts": ";".join(f"{key}:{value}" for key, value in sorted(cluster_counts.items())),
                "backend_requested": backend_requested,
                "backend_used": backend_used,
                "fallback_used": str(bool(fallback_reason)).lower(),
                "fallback_reason": fallback_reason,
                "mixed_template_warning": str(mixed_warning).lower(),
                "min_cluster_size": args.min_cluster_size,
                "min_cluster_fraction": args.min_cluster_fraction,
                "kmer_size": args.kmer_size,
            }
        )


if __name__ == "__main__":
    main()
