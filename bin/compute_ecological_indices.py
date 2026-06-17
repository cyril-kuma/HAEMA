#!/usr/bin/env python3
"""Compute vector-host ecological indices from HÆMA blood-meal host calls.

Standard malaria-vector blood-meal-ecology indices, computed strictly from molecular host
detections (no host-availability census is assumed — so availability-dependent indices such as the
forage ratio and Kay feeding index are deliberately *not* computed; see the module docs).

Indices (per mosquito = sample, hosts unioned across the three markers, controls excluded):
  - Human Blood Index (HBI): proportion of host-identified mosquitoes whose blood meal contains
    human blood. A mixed human+animal meal counts as human-positive — it represents a human bite
    (Garrett-Jones 1964, Bull. WHO 30:241-261; vectorial capacity scales as HBI^2).
  - Animal blood index / zoophily: proportion containing >=1 non-human host (overlaps HBI for
    mixed meals, so HBI + zoophily can exceed 1).
  - Feeding-type partition: human-only / mixed (human+animal) / animal-only — these DO sum to 1.
  - Mixed-feeding (multiple-blood-meal) rate: proportion feeding on >=2 distinct host taxa.
  - Host-specific blood indices: proportion of meals *containing* each host taxon (do not sum to 1).
  - Host-community diversity by mosquito incidence: richness S, Shannon H', Gini-Simpson (1-D),
    Pielou evenness J'.

Proportions carry Wilson score 95% confidence intervals (valid at small n, unlike Wald). Estimates
are reported overall and stratified by ecological zone and vector sibling species; strata are
descriptive (small n), not formal tests.

References: Garrett-Jones (1964) Bull WHO; Orsborne et al. (2018) Malar J 17:479 (HBI meta-
regression: HBI tracks location > sibling species); Hess et al. (1968) forage ratio (excluded:
needs host census); Levins (1968) niche breadth; Pianka (1973) niche overlap; Shannon (1948);
Simpson (1949); Pielou (1966).
"""
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

HUMAN = "Homo sapiens"
NON_HOST = {"", "unassigned"}
Z = 1.959963984540054  # 95% normal quantile


def read_tsv(path):
    p = Path(path) if path else None
    if not p or not p.exists() or p.stat().st_size == 0:
        return []
    with p.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def wilson_ci(x, n, z=Z):
    """Wilson score interval for a binomial proportion (robust at small n)."""
    if n == 0:
        return (None, None, None)
    p = x / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (p, max(0.0, centre - half), min(1.0, centre + half))


def diversity(counts):
    """Shannon H', Gini-Simpson (1-D), Pielou J', richness S from incidence counts (a dict)."""
    vals = [c for c in counts.values() if c > 0]
    total = sum(vals)
    s = len(vals)
    if total == 0 or s == 0:
        return {"host_richness": 0, "shannon_h": 0.0, "gini_simpson": 0.0, "pielou_evenness": 0.0}
    ps = [c / total for c in vals]
    h = max(0.0, -sum(p * math.log(p) for p in ps))  # clamp -0.0 / fp noise for single-host groups
    simpson = max(0.0, 1.0 - sum(p * p for p in ps))
    pielou = (h / math.log(s)) if s > 1 else 0.0
    return {"host_richness": s, "shannon_h": h, "gini_simpson": simpson, "pielou_evenness": pielou}


def per_mosquito_hosts(host_rows):
    """Collapse host_call_table (sample x marker x host) to one host *set* per field mosquito.

    Returns {sample_uid: set(host_taxa)} over control_status == 'sample' only; a sample present
    with no real host call maps to an empty set (tested-but-unidentified)."""
    hosts = defaultdict(set)
    for r in host_rows:
        if r.get("control_status") != "sample":
            continue
        uid = r.get("sample_uid", "")
        hosts.setdefault(uid, set())
        h = (r.get("host_assignment") or "").strip()
        if h not in NON_HOST:
            hosts[uid].add(h)
    return hosts


def index_block(uids, host_sets, stratum_type, stratum):
    """Compute the full index set for a group of mosquito sample_uids. Yields tidy rows."""
    n_tested = len(uids)
    identified = [u for u in uids if host_sets.get(u)]
    n_id = len(identified)
    rows = []

    def emit(metric, value, ci=(None, None), n=n_id, detail=""):
        lo, hi = ci
        rows.append({
            "stratum_type": stratum_type, "stratum": stratum, "metric": metric,
            "value": "" if value is None else round(value, 4),
            "ci_low": "" if lo is None else round(lo, 4),
            "ci_high": "" if hi is None else round(hi, 4),
            "n": n, "detail": detail,
        })

    small = "small_n (<5); interpret with caution" if 0 < n_id < 5 else ""
    emit("n_tested", n_tested, n=n_tested)
    emit("n_identified", n_id)
    if n_tested:
        p, lo, hi = wilson_ci(n_id, n_tested)
        emit("identification_rate", p, (lo, hi), n=n_tested,
             detail="host-identified / tested; differential ID bias possible")
    if n_id == 0:
        return rows

    contains_human = sum(1 for u in identified if HUMAN in host_sets[u])
    contains_animal = sum(1 for u in identified if any(h != HUMAN for h in host_sets[u]))
    human_only = sum(1 for u in identified if host_sets[u] == {HUMAN})
    animal_only = sum(1 for u in identified if HUMAN not in host_sets[u])
    mixed_ha = sum(1 for u in identified if HUMAN in host_sets[u] and any(h != HUMAN for h in host_sets[u]))
    multi = sum(1 for u in identified if len(host_sets[u]) >= 2)

    p, lo, hi = wilson_ci(contains_human, n_id); emit("human_blood_index", p, (lo, hi), detail="any human-containing meal; " + small if small else "any human-containing meal")
    p, lo, hi = wilson_ci(contains_animal, n_id); emit("animal_blood_index_zoophily", p, (lo, hi), detail="any non-human host")
    p, lo, hi = wilson_ci(human_only, n_id); emit("human_only_fraction", p, (lo, hi))
    p, lo, hi = wilson_ci(mixed_ha, n_id); emit("mixed_human_animal_fraction", p, (lo, hi))
    p, lo, hi = wilson_ci(animal_only, n_id); emit("animal_only_fraction", p, (lo, hi))
    p, lo, hi = wilson_ci(multi, n_id); emit("mixed_feeding_rate", p, (lo, hi), detail=">=2 distinct host taxa")

    # host-specific blood indices (proportion of identified meals containing host X; do not sum to 1)
    host_incidence = defaultdict(int)
    for u in identified:
        for h in host_sets[u]:
            host_incidence[h] += 1
    for h in sorted(host_incidence, key=lambda k: -host_incidence[k]):
        p, lo, hi = wilson_ci(host_incidence[h], n_id)
        emit(f"host_blood_index::{h}", p, (lo, hi))

    # host-community diversity by mosquito incidence (not read fractions)
    div = diversity(host_incidence)
    for k, v in div.items():
        emit(k, v, detail="by mosquito incidence")
    return rows


def main():
    ap = argparse.ArgumentParser(description="Compute vector-host ecological indices from host calls.")
    ap.add_argument("--host-calls", required=True, help="host_call_table.tsv")
    ap.add_argument("--master-endpoint", default="", help="bloodmeal_master_endpoint.tsv (for strata metadata)")
    ap.add_argument("--zone-column", default="collection_region", help="ecological-zone metadata column")
    ap.add_argument("--species-column", default="sibling_species", help="vector sibling-species column")
    ap.add_argument("--output-tsv", required=True)
    ap.add_argument("--output-json", default="")
    args = ap.parse_args()

    host_rows = read_tsv(args.host_calls)
    host_sets = per_mosquito_hosts(host_rows)
    all_uids = list(host_sets.keys())

    # map sample_uid -> strata metadata from the master endpoint (one value per uid)
    zone_of, species_of = {}, {}
    for r in read_tsv(args.master_endpoint):
        if r.get("sample_type") != "sample":
            continue
        uid = r.get("sample_uid", "")
        if uid:
            zone_of.setdefault(uid, (r.get(args.zone_column) or "").strip())
            species_of.setdefault(uid, (r.get(args.species_column) or "").strip())

    rows = []
    rows += index_block(all_uids, host_sets, "overall", "all_field_samples")
    for label, mapping in (("ecological_zone", zone_of), ("sibling_species", species_of)):
        levels = sorted({v for v in mapping.values() if v})
        for lv in levels:
            uids = [u for u in all_uids if mapping.get(u) == lv]
            if uids:
                rows += index_block(uids, host_sets, label, lv)

    fields = ["stratum_type", "stratum", "metric", "value", "ci_low", "ci_high", "n", "detail"]
    with open(args.output_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    if args.output_json:
        nested = defaultdict(lambda: defaultdict(dict))
        for r in rows:
            nested[f"{r['stratum_type']}:{r['stratum']}"][r["metric"]] = {
                "value": r["value"], "ci_low": r["ci_low"], "ci_high": r["ci_high"], "n": r["n"]}
        with open(args.output_json, "w") as fh:
            json.dump({"indices": nested,
                       "excluded": ["forage_ratio", "kay_feeding_index", "manly_selection",
                                    "vectorial_capacity"],
                       "excluded_reason": "require host-availability census or survival/biting-rate "
                                          "data not collected by this pipeline"}, fh, indent=2)

    overall = [r for r in rows if r["stratum_type"] == "overall"]
    hbi = next((r for r in overall if r["metric"] == "human_blood_index"), None)
    print(f"Ecological indices written to {args.output_tsv} ({len(rows)} rows). "
          + (f"Overall HBI = {hbi['value']} (95% CI {hbi['ci_low']}-{hbi['ci_high']}, n={hbi['n']})."
             if hbi else "No host-identified field samples."))


if __name__ == "__main__":
    main()
