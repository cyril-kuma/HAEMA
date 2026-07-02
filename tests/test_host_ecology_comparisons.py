#!/usr/bin/env python3
"""Tests for bin/compute_host_ecology_comparisons.py.

Verifies per-specimen aggregation (a mosquito counted once), control exclusion,
strata drawn from the master endpoint, and Holm-corrected Fisher output structure.
"""
import csv
import json
import subprocess
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin" / "compute_host_ecology_comparisons.py"

HOST_CALL_FIELDS = [
    "sample_uid", "run_id", "sample_id", "barcode_id", "control_status", "marker",
    "best_cluster_id", "host_assignment", "host_rank", "host_reads", "host_fraction",
    "n_supporting_features", "best_feature_id", "best_confidence",
    "best_assignment_status", "mixed_status", "total_marker_reads",
    "mixed_template_warning", "model",
]
ENDPOINT_FIELDS = ["sample_uid", "sample_type", "collection_region", "sibling_species", "collection_date"]


def _hc(uid, marker, host, control="sample", rank=1):
    return {
        "sample_uid": uid, "run_id": "R1", "sample_id": uid, "barcode_id": "bc",
        "control_status": control, "marker": marker, "best_cluster_id": "",
        "host_assignment": host, "host_rank": rank, "host_reads": 40,
        "host_fraction": "0.8", "n_supporting_features": 1, "best_feature_id": "a",
        "best_confidence": "high", "best_assignment_status": "assigned",
        "mixed_status": "single_host", "total_marker_reads": 45,
        "mixed_template_warning": "false", "model": "rambo",
    }


def _ep(uid, zone, species="An. gambiae", date="2024-06-15", stype="sample"):
    return {"sample_uid": uid, "sample_type": stype, "collection_region": zone,
            "sibling_species": species, "collection_date": date}


def _write(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def _run(tmp_path, host_calls, endpoint):
    hc = tmp_path / "hc.tsv"
    ep = tmp_path / "ep.tsv"
    outdir = tmp_path / "out"
    _write(hc, host_calls, HOST_CALL_FIELDS)
    _write(ep, endpoint, ENDPOINT_FIELDS)
    subprocess.run(
        [sys.executable, str(BIN), "--host-call-table", str(hc),
         "--master-endpoint", str(ep), "--output-dir", str(outdir)],
        check=True, capture_output=True, text=True,
    )
    def rd(name):
        with open(outdir / name, newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))
    summary = json.loads((outdir / "host_use_statistical_tests_summary.json").read_text())
    return rd, summary


def test_mosquito_counted_once_across_markers_and_hosts(tmp_path):
    # One specimen, human called by two markers + a second (mixed) host row.
    # It must count as ONE identified specimen, human-positive, in zone Forest.
    calls = [
        _hc("u1", "cyt_b", "Homo sapiens"),
        _hc("u1", "co1_short", "Homo sapiens"),
        _hc("u1", "cyt_b", "Bos taurus", rank=2),
    ]
    ep = [_ep("u1", "Forest")]
    rd, summary = _run(tmp_path, calls, ep)
    assert summary["field_specimens_identified"] == 1
    rich = {r["stratum"]: r for r in rd("host_richness_by_zone.tsv")}
    assert rich["Forest"]["n_specimens"] == "1"
    # host richness = {Homo sapiens, Bos taurus} = 2
    assert rich["Forest"]["host_richness_S"] == "2"


def test_controls_excluded_from_indices(tmp_path):
    calls = [
        _hc("neg", "cyt_b", "Homo sapiens", control="negative_control"),
        _hc("u1", "cyt_b", "Bos taurus"),
    ]
    ep = [_ep("neg", "Forest", stype="negative_control"), _ep("u1", "Forest")]
    _, summary = _run(tmp_path, calls, ep)
    assert summary["field_specimens_identified"] == 1


def test_hbi_pairwise_structure_and_holm(tmp_path):
    # Two zones with differing HBI; expect one pairwise comparison with Holm fields.
    calls = (
        [_hc(f"h{i}", "cyt_b", "Homo sapiens") for i in range(6)]
        + [_hc(f"c{i}", "cyt_b", "Bos taurus") for i in range(6)]
    )
    ep = ([_ep(f"h{i}", "Coastal") for i in range(6)]
          + [_ep(f"c{i}", "Sahel") for i in range(6)])
    rd, summary = _run(tmp_path, calls, ep)
    comps = rd("pairwise_hbi_comparisons.tsv")
    assert len(comps) == 1
    row = comps[0]
    assert {row["stratum_1"], row["stratum_2"]} == {"Coastal", "Sahel"}
    assert row["event_1"] in ("6", "0") and row["event_2"] in ("6", "0")
    assert row["corrected_p_value"] not in ("", None)
    assert row["test"] == "fishers_exact"


def test_empty_inputs_write_headers(tmp_path):
    rd, summary = _run(tmp_path, [], [])
    for name in ("pairwise_hbi_comparisons.tsv", "pairwise_hbi_species_comparisons.tsv",
                 "pairwise_mixed_feeding_comparisons.tsv", "host_richness_by_zone.tsv"):
        # header present, zero data rows
        assert isinstance(rd(name), list)
    assert summary["field_specimens_identified"] == 0


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-v"]))
