#!/usr/bin/env python3
"""Regression tests for the in-silico mixed-host DETECTION logic (bin/classify_reads.py).

Protects the classify-then-count mixed-host decision (ONT-adapted Logue 2016) against future
changes: a single-host input must not be called mixed, a two-host mixture above threshold must be
called mixed, a minor host below the calibrated fraction floor must be dropped, and the vector must
never be counted as a host. Fast and deterministic — no minimap2, no containers, no Nextflow.

These are COMPUTATIONAL detection tests only. They do NOT validate blood-volume quantification or
any biological read-fraction<->blood-proportion relationship; `host_fractions_benchmarked` stays
false (asserted below via the fixture metadata).
"""
import importlib.util
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("cr", REPO / "bin" / "classify_reads.py")
cr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cr)

DEFAULT_MIN_READS = 3
CAL_FRACTION = 0.02  # calibrated default (docs/denoising_calibration)


def _call(counts, min_fraction=CAL_FRACTION, global_min=1):
    total = sum(counts.values())
    return cr.denoise_and_call(counts, total, DEFAULT_MIN_READS, min_fraction, global_min,
                               cr.NON_HOST_DEFAULT)


def test_single_host_not_called_mixed():
    kept, call = _call({"Homo sapiens": 300})
    assert call == "single"
    assert set(kept) == {"Homo sapiens"}


def test_two_host_mixture_called_mixed():
    # 90:10 — minor host (10%) is well above the 2% floor.
    kept, call = _call({"Homo sapiens": 270, "Bos taurus": 30})
    assert call == "mixed"
    assert set(kept) == {"Homo sapiens", "Bos taurus"}


def test_minor_host_below_threshold_dropped():
    # 1% minor sits below the calibrated 2% floor -> not called (single).
    kept, call = _call({"Homo sapiens": 297, "Bos taurus": 3})
    assert call == "single"
    assert "Bos taurus" not in kept


def test_minor_host_detected_at_calibrated_floor():
    # ~2% minor with enough reads is retained.
    kept, call = _call({"Homo sapiens": 490, "Capra hircus": 10})
    assert call == "mixed"
    assert "Capra hircus" in kept


def test_vector_self_hit_excluded():
    # Anopheles (vector) mtDNA must never be counted as a host.
    kept, call = _call({"Homo sapiens": 280, "Anopheles gambiae": 20})
    assert call == "single"
    assert not any(t.startswith("Anopheles") for t in kept)


def test_global_min_count_filter():
    # A taxon seen only once is dropped regardless of fraction (Logue global filter).
    kept, _ = _call({"Homo sapiens": 50, "Ovis aries": 1}, global_min=1)
    assert "Ovis aries" not in kept


def test_paf_identity_coverage_filtering():
    # PAF: qname qlen qstart qend strand tname tlen tstart tend nmatch alnlen mapq
    paf = [
        "r1\t100\t0\t100\t+\tHomo_sapiens|NC_012920.1\t16000\t0\t100\t98\t100\t60",   # 98% id, 100% cov -> keep
        "r2\t100\t0\t100\t+\tBos_taurus|NC_006853.1\t16000\t0\t100\t50\t100\t60",     # 50% id -> drop (<80)
        "r3\t100\t0\t40\t+\tHomo_sapiens|NC_012920.1\t16000\t0\t40\t40\t40\t60",       # 40% cov -> drop (<50)
    ]
    hits = cr.best_hits_from_paf(paf, min_identity=80.0, min_coverage=50.0)
    assert hits == {"r1": "Homo sapiens"}


def test_fixture_metadata_asserts_not_benchmarked():
    meta_path = REPO / "tests" / "data" / "insilico_mixtures" / "run_metadata.json"
    assert meta_path.exists(), "in-silico fixture metadata missing"
    meta = json.loads(meta_path.read_text())
    assert meta["host_fractions_benchmarked"] is False
    assert meta["wetlab_known_ratio_validation"] is False
    assert meta["in_silico_mixture_validation"] is True
