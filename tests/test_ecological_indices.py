#!/usr/bin/env python3
"""Unit tests for vector-host ecological-index logic in bin/compute_ecological_indices.py.

Run from the repo root:  python3 tests/test_ecological_indices.py
"""
import importlib.util
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("eco", REPO / "bin" / "compute_ecological_indices.py")
eco = importlib.util.module_from_spec(spec)
spec.loader.exec_module(eco)


def row(uid, host, control="sample", marker="cyt_b"):
    return {"sample_uid": uid, "control_status": control, "host_assignment": host, "marker": marker}


def approx(a, b, tol=1e-3):
    return abs(a - b) <= tol


def get(rows, metric, col="value"):
    for r in rows:
        if r["metric"] == metric:
            return r[col]
    return None


def test_marker_union_and_control_exclusion():
    rows = [
        row("m1", "Homo sapiens", marker="cyt_b"),
        row("m1", "Homo sapiens", marker="co1_short"),       # same host, 2 markers -> one set
        row("m2", "Homo sapiens", marker="cyt_b"),
        row("m2", "Capra hircus", marker="co1_short"),       # mixed across markers
        row("m3", "Ovis aries", marker="co1_short"),
        row("m4", "unassigned", marker="cyt_b"),              # tested, unidentified
        row("c1", "Homo sapiens", control="negative_control"),  # control -> excluded
    ]
    hs = eco.per_mosquito_hosts(rows)
    assert set(hs.keys()) == {"m1", "m2", "m3", "m4"}, "controls must be excluded"
    assert hs["m1"] == {"Homo sapiens"}, "same host on 2 markers collapses to one"
    assert hs["m2"] == {"Homo sapiens", "Capra hircus"}, "hosts union across markers"
    assert hs["m4"] == set(), "unassigned -> empty (tested but unidentified)"
    print("marker-union + control-exclusion test PASSED")


def test_hbi_mixed_counting_and_partition():
    rows = [
        row("m1", "Homo sapiens"),                      # human-only
        row("m2", "Homo sapiens"), row("m2", "Capra hircus", marker="co1_short"),  # mixed
        row("m3", "Ovis aries"),                        # animal-only
        row("m4", "unassigned"),                        # unidentified
    ]
    hs = eco.per_mosquito_hosts(rows)
    out = eco.index_block(list(hs.keys()), hs, "overall", "all")

    assert get(out, "n_tested") == 4
    assert get(out, "n_identified") == 3
    # HBI counts the mixed meal as human-positive: m1 + m2 = 2/3
    assert approx(get(out, "human_blood_index"), 2 / 3), "mixed human+animal must count toward HBI"
    # zoophily: m2 + m3 = 2/3 ; HBI + zoophily > 1 because the mixed meal is in both
    assert approx(get(out, "animal_blood_index_zoophily"), 2 / 3)
    # feeding-type partition sums to exactly 1
    parts = (get(out, "human_only_fraction") + get(out, "mixed_human_animal_fraction")
             + get(out, "animal_only_fraction"))
    assert approx(parts, 1.0), "human-only/mixed/animal-only must partition to 1"
    assert approx(get(out, "mixed_feeding_rate"), 1 / 3)
    # host-specific indices do NOT sum to 1 (mixed meal counted under each host)
    hsum = (get(out, "host_blood_index::Homo sapiens") + get(out, "host_blood_index::Capra hircus")
            + get(out, "host_blood_index::Ovis aries"))
    assert hsum > 1.0, "host-specific blood indices should exceed 1 when mixed meals exist"
    print("HBI mixed-counting + partition test PASSED")


def test_wilson_ci_bounds():
    for x, n in [(0, 5), (2, 3), (5, 5), (1, 14)]:
        p, lo, hi = eco.wilson_ci(x, n)
        assert 0.0 <= lo <= p <= hi <= 1.0, f"Wilson CI out of bounds for {x}/{n}"
    # zero events -> lower bound 0, upper bound > 0 (rule-of-three-ish), never negative
    _, lo, hi = eco.wilson_ci(0, 10)
    assert lo == 0.0 and hi > 0.0
    print("Wilson CI bounds test PASSED")


def test_diversity_single_vs_even():
    single = eco.diversity({"Homo sapiens": 10})
    assert single["host_richness"] == 1 and approx(single["shannon_h"], 0.0) \
        and approx(single["gini_simpson"], 0.0), "single host -> zero diversity (no -0.0)"
    even = eco.diversity({"a": 5, "b": 5})
    assert even["host_richness"] == 2 and approx(even["pielou_evenness"], 1.0), \
        "two equal hosts -> max evenness"
    print("diversity single/even test PASSED")


if __name__ == "__main__":
    test_marker_union_and_control_exclusion()
    test_hbi_mixed_counting_and_partition()
    test_wilson_ci_bounds()
    test_diversity_single_vs_even()
    print("ecological-index tests PASSED (union, HBI mixed-counting, partition, Wilson CI, diversity)")
