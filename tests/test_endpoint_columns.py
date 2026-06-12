#!/usr/bin/env python3
"""Assert the core endpoint tables exist and carry their expected columns.

Usage: test_endpoint_columns.py <results_dir>
Exits non-zero with a clear message on any missing file or column.
"""
import csv
import sys
from pathlib import Path

# file (relative to results dir) -> columns that MUST be present
EXPECTED = {
    "05_endpoint_files/bloodmeal_master_endpoint.tsv": [
        "sample_uid", "marker", "host_assignment", "control_status", "contamination_flag",
    ],
    "05_endpoint_files/host_call_table.tsv": [
        "sample_uid", "marker", "host_assignment", "host_reads", "host_fraction",
        "mixed_status", "control_status",
    ],
    "05_endpoint_files/asv_count_table.tsv": ["feature_id", "marker", "sequence"],
    "05_endpoint_files/sample_level_summary.tsv": ["sample_uid", "sample_type", "total_asv_reads"],
    "05_endpoint_files/run_manifest.json": None,            # presence only
    "06_reports/rambo_model_summary.tsv": ["metric", "value"],
    "06_reports/positive_control_check.tsv": [
        "sample_uid", "control_kind", "expected_hosts", "observed_hosts",
        "n_expected", "n_recovered", "status",
    ],
}


def header(path):
    with open(path, newline="") as fh:
        return next(csv.reader(fh, delimiter="\t"), [])


def main():
    if len(sys.argv) != 2:
        print("usage: test_endpoint_columns.py <results_dir>", file=sys.stderr)
        sys.exit(2)
    root = Path(sys.argv[1])
    errors = []
    for rel, cols in EXPECTED.items():
        p = root / rel
        if not p.exists():
            errors.append(f"missing output: {rel}")
            continue
        if cols is None:
            continue
        have = set(header(p))
        for c in cols:
            if c not in have:
                errors.append(f"{rel}: missing column '{c}'")
    if errors:
        print("endpoint-column check FAILED:", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        sys.exit(1)
    print(f"endpoint-column check PASSED ({len(EXPECTED)} files)")


if __name__ == "__main__":
    main()
