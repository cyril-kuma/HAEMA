#!/usr/bin/env python3
"""Unit tests for input validation: malformed samplesheets/primers must fail clearly.

Run from the repo root:  python3 tests/test_validation.py
Uses only the standard library and the bundled test fixtures.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
VALIDATE = REPO / "bin" / "validate_inputs.py"
RAW = REPO / "assets" / "test_data" / "runs"
PRIMERS = REPO / "assets" / "test_data" / "primers.csv"

GOOD_HEADER = "run_id,barcode_id,sample_id,sample_type"
GOOD_ROW = "RUN_TEST,barcode01,TST001,sample"
FULL_HEADER = (
    "run_id,barcode_id,sample_id,specimen_id,sample_type,species,sibling_species,"
    "feeding_status,collection_date,collection_time,collection_location,bioclimatic_zone,"
    "collection_region,latitude,longitude,collection_context,collection_method,specimen_sex"
)
FULL_ROW = (
    "RUN_TEST,barcode01,TST001,TST001,sample,Anopheles_gambiae_s.l,Anopheles_coluzzii,"
    "Blood_fed,2026-01-01,00:00,Test_site,Forest,Test_region,0,0,Indoor,LTC,Female"
)


def run_validate(samplesheet, primers):
    with tempfile.TemporaryDirectory() as d:
        return subprocess.run(
            [sys.executable, str(VALIDATE), "--samplesheet", str(samplesheet),
             "--primers", str(primers), "--raw-data-dir", str(RAW), "--output-dir", d],
            capture_output=True, text=True,
        )


def write(tmp, name, text):
    p = Path(tmp) / name
    p.write_text(text)
    return p


def main():
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        # 1. Missing required column (barcode_id) -> must fail with a clear message.
        bad = write(tmp, "bad_missing_col.csv", "run_id,sample_id,sample_type\nRUN_TEST,TST001,sample\n")
        r = run_validate(bad, PRIMERS)
        if r.returncode == 0 or "missing required columns" not in r.stderr.lower():
            failures.append(f"missing-column samplesheet should fail clearly (rc={r.returncode})")

        # 2. Well-formed minimal samplesheet -> must pass (rc 0).
        good = write(tmp, "good.csv", f"{GOOD_HEADER}\n{GOOD_ROW}\n")
        r = run_validate(good, PRIMERS)
        if r.returncode != 0:
            failures.append(f"valid minimal samplesheet should pass (rc={r.returncode}; {r.stderr[-200:]})")

        # 3. Malformed primer file (missing Reverse_Primer column) -> must fail.
        bad_primers = write(tmp, "bad_primers.csv", "Gene,Forward_Primer,Size\ncyt_b,ACGT,450\n")
        r = run_validate(good, bad_primers)
        if r.returncode == 0 or "primer file missing required columns" not in r.stderr.lower():
            failures.append(f"malformed primer file should fail clearly (rc={r.returncode})")

        # 4. Current coordinate columns should not trigger the legacy typo-column warning.
        full = write(tmp, "full_metadata.csv", f"{FULL_HEADER}\n{FULL_ROW}\n")
        r = run_validate(full, PRIMERS)
        if r.returncode != 0:
            failures.append(f"full metadata samplesheet should pass (rc={r.returncode}; {r.stderr[-200:]})")
        if "collection_cordinates" in r.stderr:
            failures.append("latitude/longitude samplesheet should not warn about collection_cordinates")

    if failures:
        print("validation tests FAILED:", file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        sys.exit(1)
    print("validation tests PASSED (4 cases)")


if __name__ == "__main__":
    main()
