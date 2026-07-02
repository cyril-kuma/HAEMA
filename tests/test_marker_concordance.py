#!/usr/bin/env python3
"""Tests for bin/compute_marker_concordance.py against the real RAMBO host-call schema.

The RAMBO host-call table (rambo_host_call_table.tsv) columns are:
  sample_uid, run_id, sample_id, barcode_id, control_status, marker, best_cluster_id,
  host_assignment, host_rank, host_reads, host_fraction, n_supporting_features,
  best_feature_id, best_confidence, best_assignment_status, mixed_status,
  total_marker_reads, mixed_template_warning, model
"""
import csv
import json
import subprocess
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin" / "compute_marker_concordance.py"

HOST_CALL_FIELDS = [
    "sample_uid", "run_id", "sample_id", "barcode_id", "control_status", "marker",
    "best_cluster_id", "host_assignment", "host_rank", "host_reads", "host_fraction",
    "n_supporting_features", "best_feature_id", "best_confidence",
    "best_assignment_status", "mixed_status", "total_marker_reads",
    "mixed_template_warning", "model",
]


def _hc(sample_uid, marker, host, conf="high", control="sample", rank=1,
        status="assigned", mixed="single_host", sample_id="S1"):
    return {
        "sample_uid": sample_uid, "run_id": "R1", "sample_id": sample_id,
        "barcode_id": "bc01", "control_status": control, "marker": marker,
        "best_cluster_id": "", "host_assignment": host, "host_rank": rank,
        "host_reads": 50, "host_fraction": "0.9", "n_supporting_features": 1,
        "best_feature_id": "asv1", "best_confidence": conf,
        "best_assignment_status": status, "mixed_status": mixed,
        "total_marker_reads": 55, "mixed_template_warning": "false", "model": "rambo",
    }


def _write_tsv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def _run(tmp_path, host_calls):
    hc = tmp_path / "host_calls.tsv"
    ep = tmp_path / "endpoint.tsv"
    out = tmp_path / "concordance.tsv"
    summ = tmp_path / "concordance_summary.json"
    _write_tsv(hc, host_calls, HOST_CALL_FIELDS)
    _write_tsv(ep, [], ["sample_uid"])  # endpoint currently only needs to exist
    subprocess.run(
        [sys.executable, str(BIN), "--host-call-table", str(hc),
         "--endpoint-manifest", str(ep), "--output", str(out),
         "--summary-output", str(summ)],
        check=True, capture_output=True, text=True,
    )
    with open(out, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    summary = json.loads(summ.read_text())
    return {r["specimen_id"]: r for r in rows}, summary


def test_full_species_agreement(tmp_path):
    # cyt_b and co1_short both call Homo sapiens for the same specimen -> full agreement
    calls = [
        _hc("uidA", "cyt_b", "Homo sapiens"),
        _hc("uidA", "co1_short", "Homo sapiens"),
    ]
    by_spec, summary = _run(tmp_path, calls)
    assert by_spec["uidA"]["concordance_status"] == "full_species_agreement"
    assert by_spec["uidA"]["possible_numt_flag"] == "False"
    assert summary["full_species_agreement"] == 1


def test_discordant_flags_numt(tmp_path):
    # cyt_b says Bos taurus, co1_short says Homo sapiens -> discordant + NUMT caution
    calls = [
        _hc("uidB", "cyt_b", "Bos taurus"),
        _hc("uidB", "co1_short", "Homo sapiens"),
    ]
    by_spec, _ = _run(tmp_path, calls)
    assert by_spec["uidB"]["concordance_status"] == "discordant"
    assert by_spec["uidB"]["possible_numt_flag"] == "True"
    assert by_spec["uidB"]["possible_mixed_meal_flag"] == "True"


def test_genus_agreement(tmp_path):
    calls = [
        _hc("uidC", "co1_short", "Canis familiaris"),
        _hc("uidC", "co1_long", "Canis lupus"),
    ]
    by_spec, _ = _run(tmp_path, calls)
    assert by_spec["uidC"]["concordance_status"] == "genus_agreement"


def test_single_marker_only(tmp_path):
    calls = [_hc("uidD", "cyt_b", "Homo sapiens")]
    by_spec, _ = _run(tmp_path, calls)
    assert by_spec["uidD"]["concordance_status"] == "single_marker_only"


def test_controls_excluded(tmp_path):
    # A control specimen must never appear in concordance output.
    calls = [
        _hc("ctrlNeg", "cyt_b", "Homo sapiens", control="negative_control"),
        _hc("uidE", "cyt_b", "Homo sapiens"),
    ]
    by_spec, _ = _run(tmp_path, calls)
    assert "ctrlNeg" not in by_spec
    assert "uidE" in by_spec


def test_rank1_dominates_mixed(tmp_path):
    # Mixed-host marker: two rows same marker; rank-1 host is used for concordance.
    calls = [
        _hc("uidF", "cyt_b", "Homo sapiens", rank=1, mixed="mixed_host"),
        _hc("uidF", "cyt_b", "Bos taurus", rank=2, mixed="mixed_host"),
        _hc("uidF", "co1_short", "Homo sapiens", rank=1),
    ]
    by_spec, _ = _run(tmp_path, calls)
    assert by_spec["uidF"]["concordance_status"] == "full_species_agreement"


def test_empty_input_writes_header(tmp_path):
    by_spec, summary = _run(tmp_path, [])
    assert by_spec == {}
    assert summary["total_samples"] == 0


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-v"]))
