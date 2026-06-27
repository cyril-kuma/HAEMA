#!/usr/bin/env python3
"""
HÆMA figure data-prep: turn pipeline output tables into small, plot-ready tables so the
plotting scripts only pivot/sort/label (no filtering or derivation at plot time).

Primary input is the rich per-ASV table `host_assignments.tsv` (host_assignment, marker,
pident, coverage, assignment_status, control_status, count + full sample metadata) plus the
per-read `*.read_decisions.tsv` QC tables and the (validated) samplesheet for vector
composition over ALL collected mosquitoes.

Emits to --outdir:
  marker_recovery.tsv, assignment_confidence.tsv, read_length_distribution.tsv,
  sibling_species_by_site.tsv, mixed_combinations.tsv, feeding_multiplicity.tsv,
  rarefaction.tsv, host_overall.tsv
"""
import argparse, csv, glob, math, os, random
from collections import defaultdict, Counter

# A feature is an accepted host call iff host_assignment is a real taxon (not unassigned).
# (assignment_status is assigned_taxid_exact_lca vs no_confident_blast_hit -> unassigned.)
NONHOST = {"", "unassigned", "no_host", "none"}
MARKERS = ["cyt_b", "co1_short", "co1_long"]
AMPLICON_SIZE = {"cyt_b": 350, "co1_short": 220, "co1_long": 650}  # nominal target sizes (bp)
COMMON = {"Homo sapiens": "Human", "Bos taurus": "Cattle", "Ovis aries": "Sheep",
          "Capra hircus": "Goat", "Canis lupus familiaris": "Dog", "Sus scrofa": "Pig",
          "Equus asinus": "Donkey", "Gallus gallus": "Chicken", "Felis catus": "Cat"}


def w(path, header, rows):
    with open(path, "w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t")
        wr.writerow(header)
        wr.writerows(rows)
    print(f"  wrote {os.path.basename(path)} ({len(rows)} rows)")


def load_assignments(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(x, d=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def per_sample_hosts(rows, count_field):
    """field sample_uid -> {host: summed reads}. Use the rambo host_call_table (count_field
    'host_reads') for ecological aggregations so figures match the published indices."""
    agg = defaultdict(lambda: defaultdict(float))
    for r in rows:
        if r.get("control_status") != "sample":
            continue
        h = (r.get("host_assignment") or "").strip()
        if h.lower() in NONHOST:
            continue
        agg[r["sample_uid"]][h] += fnum(r.get(count_field))
    return agg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host-assignments", required=True, help="host_assignments.tsv (per-ASV, for Fig 2)")
    ap.add_argument("--host-calls", required=True, help="host_call_table.tsv (rambo-filtered, ecology)")
    ap.add_argument("--read-decisions-glob", default="", help="glob for *.read_decisions.tsv")
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    out = lambda n: os.path.join(a.outdir, n)
    rows = load_assignments(a.host_assignments)
    hc_rows = load_assignments(a.host_calls)
    field = [r for r in rows if r.get("control_status") == "sample"]
    accepted = [r for r in field if (r.get("host_assignment") or "").lower() not in NONHOST]

    # 1. marker recovery: host x marker -> n distinct samples
    rec = defaultdict(set)
    for r in accepted:
        rec[(r["host_assignment"], r.get("marker", ""))].add(r["sample_uid"])
    hosts = sorted({h for h, _ in rec}, key=lambda h: -len({s for (hh, m), ss in rec.items()
                   if hh == h for s in ss}))
    w(out("marker_recovery.tsv"), ["host", "common", "marker", "n_samples"],
      [[h, COMMON.get(h, h), m, len(rec[(h, m)])] for h in hosts for m in MARKERS if (h, m) in rec])

    # 2. assignment confidence: per accepted feature projection
    w(out("assignment_confidence.tsv"), ["host", "common", "marker", "pident", "coverage"],
      [[r["host_assignment"], COMMON.get(r["host_assignment"], r["host_assignment"]),
        r.get("marker", ""), fnum(r.get("pident")), fnum(r.get("coverage"))]
       for r in accepted if r.get("pident") not in (None, "")])

    # 3. read-length distribution by marker x assigned-status (binned)
    if a.read_decisions_glob:
        binw = 25
        dist = Counter()  # (marker_or_unassigned, bin_mid) -> count
        for fp in glob.glob(a.read_decisions_glob):
            with open(fp) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    L = fnum(r.get("trimmed_length") or r.get("raw_length"))
                    if L <= 0:
                        continue
                    mk = (r.get("marker") or "").strip()
                    key = mk if mk in MARKERS else "unassigned"
                    dist[(key, int(L // binw) * binw + binw // 2)] += 1
        w(out("read_length_distribution.tsv"), ["series", "length_bp", "count"],
          [[k, b, c] for (k, b), c in sorted(dist.items())])
        w(out("amplicon_sizes.tsv"), ["marker", "expected_bp"],
          [[m, AMPLICON_SIZE[m]] for m in MARKERS])

    # 4. sibling-species composition + site coordinates (ALL collected mosquitoes, samplesheet)
    site_sp = Counter()
    site_n = Counter()
    site_ll = {}
    with open(a.samplesheet) as fh:
        for r in csv.DictReader(fh):
            site = (r.get("collection_location") or "").strip()
            zone = (r.get("bioclimatic_zone") or "").strip()
            sp = (r.get("sibling_species") or "").strip()
            if not site or zone in ("", "#N/A"):
                continue
            site_n[(zone, site)] += 1
            if sp:
                site_sp[(zone, site, sp)] += 1
            lat, lon = fnum(r.get("latitude"), None), fnum(r.get("longitude"), None)
            if lat and lon and site not in site_ll:
                site_ll[site] = (zone, lat, lon)
    w(out("sibling_species_by_site.tsv"), ["zone", "site", "sibling_species", "n"],
      [[z, s, sp, n] for (z, s, sp), n in sorted(site_sp.items())])
    w(out("sites.tsv"), ["site", "zone", "latitude", "longitude", "n_mosquitoes"],
      [[s, z, lat, lon, site_n[(z, s)]] for s, (z, lat, lon) in sorted(site_ll.items())])

    # 5/6. mixed-host combinations + feeding multiplicity (rambo host calls -> per-sample host sets)
    agg = per_sample_hosts(hc_rows, "host_reads")
    mult = Counter(len(v) for v in agg.values())
    w(out("feeding_multiplicity.tsv"), ["n_hosts", "n_samples"],
      [[k, mult[k]] for k in sorted(mult)])
    combos = Counter(tuple(sorted(v)) for v in agg.values() if len(v) >= 2)
    crow = []
    for combo, c in combos.most_common():
        label = " + ".join(COMMON.get(x, x) for x in combo)
        typ = "animal_animal" if "Homo sapiens" not in combo else "human_animal"
        crow.append([label, len(combo), typ, c])
    w(out("mixed_combinations.tsv"), ["combination", "n_hosts", "type", "count"], crow)

    # 7. rarefaction / host accumulation (mean over permutations)
    sample_sets = [set(v) for v in agg.values() if v]
    R, n = 100, len(sample_sets)
    means, sds = [], []
    rng = random.Random(42)
    for k in range(1, n + 1):
        vals = []
        for _ in range(R):
            idx = rng.sample(range(n), k)
            seen = set().union(*(sample_sets[i] for i in idx))
            vals.append(len(seen))
        mu = sum(vals) / R
        means.append(mu)
        sds.append((sum((v - mu) ** 2 for v in vals) / R) ** 0.5)
    w(out("rarefaction.tsv"), ["n_samples", "mean_taxa", "sd_taxa"],
      [[k + 1, round(means[k], 3), round(sds[k], 3)] for k in range(n)])

    # 8. overall host detection frequency (field, accepted)
    n_id = len([u for u, v in agg.items() if v])
    host_n = Counter()
    for v in agg.values():
        for h in v:
            host_n[h] += 1
    w(out("host_overall.tsv"), ["host", "common", "n_detected", "proportion"],
      [[h, COMMON.get(h, h), host_n[h], round(host_n[h] / max(n_id, 1), 4)]
       for h in sorted(host_n, key=lambda x: -host_n[x])])
    print(f"figure_data_prep: {len(field)} field features, {n_id} host-identified samples.")


if __name__ == "__main__":
    main()
