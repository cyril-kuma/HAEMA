#!/usr/bin/env python3
import argparse
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path

from fastq_utils import find_terminal_pair, mean_q, read_fastq, revcomp, write_fastq_record


def str_to_bool(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def read_marker_config(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def length_marker(seq_len, markers):
    hits = []
    for marker in markers:
        min_len = int(marker["min_len"])
        max_len = int(marker["max_len"])
        if min_len <= seq_len <= max_len:
            midpoint = (min_len + max_len) / 2
            hits.append((abs(seq_len - midpoint), marker))
    if not hits:
        return None
    return sorted(hits, key=lambda item: item[0])[0][1]


def choose_primer_marker(seq, qual, markers, window, max_error_rate):
    candidates = []
    for marker in markers:
        pair = find_terminal_pair(
            seq,
            marker["forward_primer"],
            marker["reverse_primer"],
            window,
            max_error_rate,
        )
        if not pair:
            continue
        trimmed_seq = seq[pair["start"] : pair["end"]]
        trimmed_qual = qual[pair["start"] : pair["end"]]
        if pair["orientation"] == "reverse":
            trimmed_seq = revcomp(trimmed_seq)
            trimmed_qual = trimmed_qual[::-1]
        candidates.append((pair["mismatches"], -len(trimmed_seq), marker, pair, trimmed_seq, trimmed_qual))
    if not candidates:
        return None
    return sorted(candidates, key=lambda item: (item[0], item[1]))[0]


def main():
    parser = argparse.ArgumentParser(description="Trim primers, quality-filter, and split reads by marker")
    parser.add_argument("--sample-uid", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--marker-config", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--min-mean-q", type=float, required=True)
    parser.add_argument("--min-read-length", type=int, required=True)
    parser.add_argument("--primer-search-window", type=int, required=True)
    parser.add_argument("--primer-max-error-rate", type=float, required=True)
    parser.add_argument("--allow-length-fallback", default="true")
    args = parser.parse_args()

    allow_length_fallback = str_to_bool(args.allow_length_fallback)
    markers = read_marker_config(args.marker_config)
    marker_by_name = {m["marker"]: m for m in markers}
    writers = {}
    try:
        for marker in markers:
            writers[marker["marker"]] = gzip.open(f"{args.output_prefix}.{marker['marker']}.fastq.gz", "wt")

        summary = defaultdict(Counter)
        decisions = []
        raw_reads = 0

        for header, seq, _plus, qual in read_fastq(args.input):
            raw_reads += 1
            read_id = header[1:].split()[0] if header.startswith("@") else header.split()[0]
            q = mean_q(qual)
            raw_len = len(seq)
            for marker in markers:
                summary[marker["marker"]]["raw_reads"] += 1

            if q < args.min_mean_q:
                decisions.append(
                    [args.sample_uid, read_id, "filtered_low_quality", "", "", raw_len, 0, f"{q:.2f}"]
                )
                continue
            for marker in markers:
                summary[marker["marker"]]["pass_quality"] += 1

            if raw_len < args.min_read_length:
                decisions.append([args.sample_uid, read_id, "filtered_short", "", "", raw_len, 0, f"{q:.2f}"])
                continue

            primer_choice = choose_primer_marker(
                seq,
                qual,
                markers,
                args.primer_search_window,
                args.primer_max_error_rate,
            )
            assignment_method = ""
            chosen_marker = None
            trimmed_seq = seq
            trimmed_qual = qual
            if primer_choice:
                _mismatches, _neg_len, marker, _pair, trimmed_seq, trimmed_qual = primer_choice
                chosen_marker = marker
                assignment_method = "primer_pair"
                summary[marker["marker"]]["assigned_by_primer"] += 1
            elif allow_length_fallback:
                marker = length_marker(len(seq), markers)
                if marker:
                    chosen_marker = marker
                    assignment_method = "length_fallback"
                    summary[marker["marker"]]["assigned_by_length"] += 1

            if not chosen_marker:
                for marker in markers:
                    summary[marker["marker"]]["no_marker"] += 1
                decisions.append([args.sample_uid, read_id, "no_marker", "", "", raw_len, 0, f"{q:.2f}"])
                continue

            marker_name = chosen_marker["marker"]
            min_len = int(chosen_marker["min_len"])
            max_len = int(chosen_marker["max_len"])
            trimmed_len = len(trimmed_seq)
            if not (min_len <= trimmed_len <= max_len):
                summary[marker_name]["filtered_length"] += 1
                decisions.append(
                    [
                        args.sample_uid,
                        read_id,
                        "filtered_marker_length",
                        marker_name,
                        assignment_method,
                        raw_len,
                        trimmed_len,
                        f"{q:.2f}",
                    ]
                )
                continue

            out_header = f"{header} marker={marker_name};assignment={assignment_method}"
            write_fastq_record(writers[marker_name], out_header, trimmed_seq, trimmed_qual)
            summary[marker_name]["written_reads"] += 1
            decisions.append(
                [
                    args.sample_uid,
                    read_id,
                    "pass",
                    marker_name,
                    assignment_method,
                    raw_len,
                    trimmed_len,
                    f"{q:.2f}",
                ]
            )
    finally:
        for writer in writers.values():
            writer.close()

    with Path(f"{args.output_prefix}.trim_filter_split_summary.tsv").open("w", newline="") as handle:
        fieldnames = [
            "sample_uid",
            "marker",
            "raw_reads",
            "pass_quality",
            "assigned_by_primer",
            "assigned_by_length",
            "written_reads",
            "filtered_length",
            "no_marker",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for marker in marker_by_name:
            row = {"sample_uid": args.sample_uid, "marker": marker}
            for field in fieldnames:
                if field not in row:
                    row[field] = summary[marker].get(field, 0)
            writer.writerow(row)

    with Path(f"{args.output_prefix}.read_decisions.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample_uid",
                "read_id",
                "status",
                "marker",
                "assignment_method",
                "raw_length",
                "trimmed_length",
                "mean_q",
            ]
        )
        writer.writerows(decisions)


if __name__ == "__main__":
    main()
