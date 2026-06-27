#!/usr/bin/env python3
"""Build host-use ecology tables for the curated HÆMA publication figures.

The main ecological-index step already reports standard vector indices in a tidy table. This
helper produces the additional plot-ready matrices needed by the Objective 1 publication figure
suite: vector-host incidence, Levins niche breadth, Pianka overlap, Bray-Curtis turnover, and a
simple network summary. The unit is a mosquito sample; a mixed meal contributes one incidence to
each detected host.
"""
import argparse
import csv
import math
import random
import re
from collections import Counter, defaultdict
from itertools import combinations_with_replacement
from pathlib import Path

HUMAN = "Homo sapiens"
NON_HOST = {"", "unassigned", "ambiguous", "none", "no_host"}
HOST_ORDER = [
    "Homo sapiens", "Bos taurus", "Ovis aries", "Capra hircus",
    "Canis lupus familiaris", "Sus scrofa", "Equus asinus", "Gallus gallus",
    "Felis catus",
]
SPECIES_ORDER = ["Anopheles_coluzzii", "Anopheles_gambiae_s.s", "Anopheles_arabiensis"]
ZONE_ORDER = ["Coastal_Savannah", "Forest", "Northern_Savannah"]
Z = 1.959963984540054


def sanitize(value):
    value = str(value or "").strip()
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("_") or "missing"


def read_delimited(path, delimiter="\t"):
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return []
    with p.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_tsv(path, rows, fieldnames):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fnum(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def canonical_species(value):
    value = (value or "").strip()
    aliases = {
        "An_coluzzii": "Anopheles_coluzzii",
        "An. coluzzii": "Anopheles_coluzzii",
        "An_gambiae_s.s": "Anopheles_gambiae_s.s",
        "An. gambiae s.s.": "Anopheles_gambiae_s.s",
        "An_arabiensis": "Anopheles_arabiensis",
        "An. arabiensis": "Anopheles_arabiensis",
    }
    return aliases.get(value, value or "Unknown")


def ordered(values, preferred):
    present = set(values)
    return [v for v in preferred if v in present] + sorted(present - set(preferred))


def wilson_ci(x, n, z=Z):
    if n == 0:
        return (None, None, None)
    p = x / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (p, max(0.0, centre - half), min(1.0, centre + half))


def sample_metadata(samplesheet):
    rows = read_delimited(samplesheet, delimiter=",")
    meta = {}
    for row in rows:
        uid = sanitize(f"{row.get('run_id', '')}__{row.get('barcode_id', '')}__{row.get('sample_id', '')}")
        if (row.get("sample_type") or "").strip() != "sample":
            continue
        meta[uid] = {
            "species": canonical_species(row.get("sibling_species", "")),
            "zone": (row.get("bioclimatic_zone") or row.get("collection_region") or "Unknown").strip() or "Unknown",
        }
    return meta


def per_sample_hosts(host_calls):
    hosts = defaultdict(set)
    for row in read_delimited(host_calls):
        if row.get("control_status") != "sample":
            continue
        uid = row.get("sample_uid", "")
        if not uid:
            continue
        hosts.setdefault(uid, set())
        host = (row.get("host_assignment") or "").strip()
        if host.lower() not in NON_HOST:
            hosts[uid].add(host)
    return hosts


def host_universe(host_sets):
    hosts = {h for hs in host_sets.values() for h in hs}
    return ordered(hosts, HOST_ORDER)


def incidence_counts(uids, host_sets, hosts):
    counts = Counter()
    for uid in uids:
        for host in host_sets.get(uid, set()):
            counts[host] += 1
    return [counts[h] for h in hosts]


def levins_ba(counts):
    total = sum(counts)
    n_hosts = len(counts)
    if total <= 0 or n_hosts <= 1:
        return 0.0
    ps = [c / total for c in counts]
    denom = sum(p * p for p in ps)
    if denom <= 0:
        return 0.0
    breadth = 1.0 / denom
    return max(0.0, min(1.0, (breadth - 1.0) / (n_hosts - 1.0)))


def pianka(counts_a, counts_b):
    sa = sum(counts_a)
    sb = sum(counts_b)
    if sa <= 0 or sb <= 0:
        return 0.0
    pa = [c / sa for c in counts_a]
    pb = [c / sb for c in counts_b]
    denom = math.sqrt(sum(p * p for p in pa) * sum(p * p for p in pb))
    return 0.0 if denom == 0 else sum(a * b for a, b in zip(pa, pb)) / denom


def bray_curtis(counts_a, counts_b):
    denom = sum(counts_a) + sum(counts_b)
    if denom <= 0:
        return 0.0
    return sum(abs(a - b) for a, b in zip(counts_a, counts_b)) / denom


def bootstrap(values, stat_fn, n_iter=400, seed=42):
    if len(values) < 2:
        v = stat_fn(values)
        return (v, v, v)
    rng = random.Random(seed)
    obs = stat_fn(values)
    stats = []
    for _ in range(n_iter):
        sample = [values[rng.randrange(len(values))] for _ in values]
        stats.append(stat_fn(sample))
    stats.sort()
    lo = stats[int(0.025 * (len(stats) - 1))]
    hi = stats[int(0.975 * (len(stats) - 1))]
    return (obs, lo, hi)


def network_h2prime_approx(matrix):
    """Entropy-scaled network specialisation proxy on [0, 1].

    True H2' is a null-model-corrected bipartite specialisation statistic. For a dependency-free
    pipeline helper we report a transparent entropy-scaled proxy: 0 means edge weights are evenly
    spread over all possible vector-host cells; 1 means all incidence is concentrated in one cell.
    """
    weights = [v for row in matrix for v in row if v > 0]
    n_cells = sum(1 for row in matrix for _ in row)
    total = sum(weights)
    if total <= 0 or n_cells <= 1:
        return 0.0
    entropy = -sum((w / total) * math.log(w / total) for w in weights)
    return max(0.0, min(1.0, 1.0 - entropy / math.log(n_cells)))


def main():
    parser = argparse.ArgumentParser(description="Compute publication host-use ecology tables.")
    parser.add_argument("--host-calls", required=True)
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    meta = sample_metadata(args.samplesheet)
    host_sets = per_sample_hosts(args.host_calls)
    identified = sorted(uid for uid, hosts in host_sets.items() if hosts and uid in meta)
    hosts = host_universe({uid: host_sets[uid] for uid in identified})

    species_levels = ordered((meta[uid]["species"] for uid in identified), SPECIES_ORDER)
    zone_levels = ordered((meta[uid]["zone"] for uid in identified), ZONE_ORDER)

    # Species x host incidence matrix.
    matrix_rows = []
    species_counts = {}
    for species in species_levels:
        uids = [uid for uid in identified if meta[uid]["species"] == species]
        counts = incidence_counts(uids, host_sets, hosts)
        species_counts[species] = counts
        matrix_rows.append({"species": species, **{host: counts[i] for i, host in enumerate(hosts)}})
    write_tsv(outdir / "vector_host_matrix.tsv", matrix_rows, ["species"] + hosts)

    # Overall + zone zooprophylaxis, species niche breadth, and network summary.
    index_rows = []

    def add_index(index, stratum, value, ci_low="", ci_high="", n="", detail=""):
        index_rows.append({
            "index": index,
            "stratum": stratum,
            "value": "" if value is None else round(value, 4),
            "ci_low": "" if ci_low in (None, "") else round(ci_low, 4),
            "ci_high": "" if ci_high in (None, "") else round(ci_high, 4),
            "n": n,
            "detail": detail,
        })

    for stratum, uids in [("Overall", identified)] + [
        (zone, [uid for uid in identified if meta[uid]["zone"] == zone]) for zone in zone_levels
    ]:
        n = len(uids)
        animal_positive = sum(1 for uid in uids if any(host != HUMAN for host in host_sets[uid]))
        value, lo, hi = wilson_ci(animal_positive, n)
        add_index("zooprophylaxis_index", stratum, value, lo, hi, n, "any non-human host among host-identified meals")

    for species in species_levels:
        sample_host_sets = [host_sets[uid] for uid in identified if meta[uid]["species"] == species]

        def stat(sample_sets):
            counts = Counter()
            for hs in sample_sets:
                for host in hs:
                    counts[host] += 1
            return levins_ba([counts[h] for h in hosts])

        value, lo, hi = bootstrap(sample_host_sets, stat)
        add_index("niche_breadth_BA", species, value, lo, hi, len(sample_host_sets), "Levins standardized niche breadth")

    matrix = [species_counts[s] for s in species_levels]
    possible = len(species_levels) * len(hosts)
    realised = sum(1 for row in matrix for value in row if value > 0)
    add_index("network_connectance", "Overall", realised / possible if possible else 0.0, n=possible, detail="nonzero vector-host cells / possible cells")
    add_index("network_H2prime_proxy", "Overall", network_h2prime_approx(matrix), detail="entropy-scaled H2prime proxy")
    write_tsv(outdir / "host_ecology_indices.tsv", index_rows, ["index", "stratum", "value", "ci_low", "ci_high", "n", "detail"])

    # Pianka niche overlap among sibling species.
    overlap_rows = []
    for species in species_levels:
        row = {"species": species}
        for other in species_levels:
            row[other] = round(pianka(species_counts.get(species, []), species_counts.get(other, [])), 4)
        overlap_rows.append(row)
    write_tsv(outdir / "niche_overlap_matrix.tsv", overlap_rows, ["species"] + species_levels)

    # Bray-Curtis host-use turnover among bioclimatic zones, including diagonal rows so a single
    # zone still yields a valid square matrix for the figure script.
    zone_counts = {}
    for zone in zone_levels:
        uids = [uid for uid in identified if meta[uid]["zone"] == zone]
        zone_counts[zone] = incidence_counts(uids, host_sets, hosts)
    beta_rows = []
    for zone_a, zone_b in combinations_with_replacement(zone_levels, 2):
        beta_rows.append({
            "zone_a": zone_a,
            "zone_b": zone_b,
            "bray_curtis": round(bray_curtis(zone_counts[zone_a], zone_counts[zone_b]), 4),
        })
    write_tsv(outdir / "beta_diversity_matrix.tsv", beta_rows, ["zone_a", "zone_b", "bray_curtis"])

    print(
        f"Host-use ecology tables written to {outdir} "
        f"({len(identified)} host-identified samples, {len(species_levels)} species, {len(hosts)} hosts)."
    )


if __name__ == "__main__":
    main()
