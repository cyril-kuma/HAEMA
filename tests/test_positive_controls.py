#!/usr/bin/env python3
"""Unit tests for positive-control checking logic in bin/rambo_mixed_model.py.

Run from the repo root:  python3 tests/test_positive_controls.py
"""
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rambo", REPO / "bin" / "rambo_mixed_model.py")
rambo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rambo)


def host_call(uid, host, control="positive_control", marker="cyt_b", reads=100):
    return {"sample_uid": uid, "control_status": control, "host_assignment": host,
            "host_reads": reads, "marker": marker}


def master_row(uid, expected, control="positive_control"):
    return {"sample_uid": uid, "control_status": control,
            "expected_host_scientific_name": expected, "sample_type": control}


def main():
    failures = []

    # taxon_match: species, genus, none
    cases = [
        ("Bos taurus", "Bos taurus", "species"),
        ("Bos taurus", "Bos indicus", "genus"),
        ("Bos taurus", "Homo sapiens", ""),
        ("bos_taurus", "Bos taurus", "species"),  # normalisation
    ]
    for exp, obs, want in cases:
        got = rambo.taxon_match(exp, obs)
        if got != want:
            failures.append(f"taxon_match({exp!r},{obs!r}) = {got!r}, want {want!r}")

    # split_expected_hosts: single + ';'/'+' separated
    if rambo.split_expected_hosts("Homo sapiens; Capra hircus") != ["Homo sapiens", "Capra hircus"]:
        failures.append("split_expected_hosts did not parse a ';'-separated mixture")

    # SINGLE-host controls: pass / genus / fail / indeterminate
    rows = [
        master_row("PC_PASS", "Bos taurus"),
        master_row("PC_GENUS", "Bos taurus"),
        master_row("PC_FAIL", "Mus musculus"),
        master_row("PC_INDET", ""),
        # MIXED-host lab-prepared controls (sample_type=sample, multi-host expected)
        master_row("MIX_ALL", "Homo sapiens; Capra hircus", control="sample"),
        master_row("MIX_PARTIAL", "Homo sapiens; Capra hircus", control="sample"),
        master_row("MIX_EXTRA", "Homo sapiens; Capra hircus", control="sample"),
        master_row("PLAIN", "", control="sample"),  # no expected -> ignored
    ]
    host_calls = [
        host_call("PC_PASS", "Bos taurus"),
        host_call("PC_GENUS", "Bos indicus"),
        host_call("PC_FAIL", "Homo sapiens"),
        host_call("MIX_ALL", "Homo sapiens", control="sample"),
        host_call("MIX_ALL", "Capra hircus", control="sample", marker="co1_short"),
        host_call("MIX_PARTIAL", "Homo sapiens", control="sample"),       # missing Capra
        host_call("MIX_EXTRA", "Homo sapiens", control="sample"),
        host_call("MIX_EXTRA", "Capra hircus", control="sample", marker="co1_short"),
        host_call("MIX_EXTRA", "Bos taurus", control="sample", marker="co1_long"),  # unexpected
        host_call("PLAIN", "Homo sapiens", control="sample"),
    ]
    check_rows, counts = rambo.check_positive_controls(rows, host_calls)
    by_uid = {r["sample_uid"]: r for r in check_rows}
    expected_status = {
        "PC_PASS": "pass",
        "PC_GENUS": "pass_genus",
        "PC_FAIL": "fail_unexpected_host",
        "PC_INDET": "indeterminate_no_expected_host",
        "MIX_ALL": "mixed_pass_all",
        "MIX_PARTIAL": "mixed_partial",
        "MIX_EXTRA": "mixed_pass_with_extra",
    }
    for uid, want in expected_status.items():
        got = by_uid.get(uid, {}).get("status")
        if got != want:
            failures.append(f"control {uid}: status {got!r}, want {want!r}")
    if "PLAIN" in by_uid:
        failures.append("sample without expected host should not be checked")
    # MIX_PARTIAL must record the missing host
    if "Capra hircus" not in (by_uid.get("MIX_PARTIAL", {}).get("missing_hosts") or ""):
        failures.append("MIX_PARTIAL should report Capra hircus as missing")
    # MIX_EXTRA must record the unexpected host
    if "Bos taurus" not in (by_uid.get("MIX_EXTRA", {}).get("unexpected_hosts") or ""):
        failures.append("MIX_EXTRA should report Bos taurus as unexpected")
    if by_uid.get("MIX_ALL", {}).get("control_kind") != "mixed_host_control":
        failures.append("MIX_ALL should be classified as mixed_host_control")

    if failures:
        print("positive-control tests FAILED:", file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        sys.exit(1)
    print(f"positive-control tests PASSED ({len(cases)} match cases + single + mixed control cases)")


if __name__ == "__main__":
    main()
