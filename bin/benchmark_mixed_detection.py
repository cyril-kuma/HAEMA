#!/usr/bin/env python3
"""Aggregate an in-silico known-ratio mixture experiment into detection/threshold curves.

Given a combined minimap2 PAF over all synthetic mixtures (read ids prefixed `<mixture_id>|...`) and
the mixture truth manifest, computes, per true-minor-fraction rung and per candidate minor-host
read-fraction threshold: minor-host detection rate (sensitivity), false-positive rate (an
unexpected third taxon kept), and observed-vs-true minor fraction. This calibrates
`min_host_fraction` for the read-classification denoising mode WITHOUT wet-lab controls.

See docs/denoising_redesign_plan.md §6. Scope: bioinformatic detection floor only — not
amplification/mtDNA/blood-volume bias.
"""
import argparse
from collections import defaultdict

from classify_reads import best_hits_from_paf, taxon_from_target  # noqa: F401


def load_manifest(path):
    rows = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            rows[r["mixture_id"]] = r
    return rows


def norm(t):
    return t.replace("_", " ").strip()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--paf", required=True, help="combined PAF, read ids prefixed '<mixture_id>|...'")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--min-identity", type=float, default=80.0)
    ap.add_argument("--min-coverage", type=float, default=50.0)
    ap.add_argument("--min-reads-per-host", type=int, default=3)
    ap.add_argument("--thresholds", default="0.005,0.01,0.02,0.05,0.10")
    ap.add_argument("--output-prefix", required=True)
    args = ap.parse_args()

    manifest = load_manifest(args.manifest)
    thresholds = [float(t) for t in args.thresholds.split(",")]

    # Split PAF lines by mixture prefix.
    by_mix = defaultdict(list)
    with open(args.paf) as fh:
        for line in fh:
            if not line.strip():
                continue
            mix = line.split("\t", 1)[0].split("|", 1)[0]
            by_mix[mix].append(line)

    # Per mixture: taxon read counts (best hit per read).
    per_mix_counts = {}
    for mix, lines in by_mix.items():
        read_taxon = best_hits_from_paf(lines, args.min_identity, args.min_coverage)
        counts = defaultdict(int)
        for t in read_taxon.values():
            counts[norm(t)] += 1
        per_mix_counts[mix] = counts

    # Detailed per-mixture observed-vs-true (at the finest threshold for reporting).
    detail_path = f"{args.output_prefix}_per_mixture.tsv"
    with open(detail_path, "w") as fh:
        fh.write("mixture_id\tmajor_host\tminor_host\ttrue_minor_fraction\tclassified_reads\t"
                 "observed_minor_fraction\tminor_detected_at_0.01\tn_taxa_kept_at_0.01\n")
        for mix, m in sorted(manifest.items()):
            counts = per_mix_counts.get(mix, {})
            total = sum(counts.values())
            minor = norm(m["minor_host"])
            obs_minor = counts.get(minor, 0) / total if total else 0.0
            kept = [t for t, c in counts.items()
                    if c >= args.min_reads_per_host and (c / total if total else 0) >= 0.01]
            fh.write(f"{mix}\t{m['major_host']}\t{m['minor_host']}\t{m['true_minor_fraction']}\t"
                     f"{total}\t{obs_minor:.4f}\t{minor in kept}\t{len(kept)}\n")

    # Sweep: detection + false-positive rate by (true rung, threshold).
    rung_thr = defaultdict(lambda: {"n": 0, "det": defaultdict(int), "fp": defaultdict(int)})
    for mix, m in manifest.items():
        counts = per_mix_counts.get(mix, {})
        total = sum(counts.values())
        if not total:
            continue
        rung = m["true_minor_fraction"]
        major, minor = norm(m["major_host"]), norm(m["minor_host"])
        cell = rung_thr[rung]
        cell["n"] += 1
        for thr in thresholds:
            kept = {t for t, c in counts.items()
                    if c >= args.min_reads_per_host and c / total >= thr}
            if minor in kept:
                cell["det"][thr] += 1
            # false positive: any kept taxon that is neither the intended major nor minor
            if kept - {major, minor}:
                cell["fp"][thr] += 1

    sweep_path = f"{args.output_prefix}_sweep.tsv"
    with open(sweep_path, "w") as fh:
        fh.write("true_minor_fraction\tn_mixtures\tmin_host_fraction_threshold\t"
                 "minor_detection_rate\tfalse_positive_rate\n")
        for rung in sorted(rung_thr, key=float):
            cell = rung_thr[rung]
            n = cell["n"]
            for thr in thresholds:
                det = cell["det"][thr] / n if n else 0.0
                fp = cell["fp"][thr] / n if n else 0.0
                fh.write(f"{rung}\t{n}\t{thr}\t{det:.3f}\t{fp:.3f}\n")

    print(f"[benchmark] {len(manifest)} mixtures; per-mixture -> {detail_path}; sweep -> {sweep_path}")
    # Console summary at the default 0.01 threshold.
    print("[benchmark] minor-host detection rate @ threshold=0.01:")
    for rung in sorted(rung_thr, key=float):
        cell = rung_thr[rung]
        n = cell["n"]
        det = cell["det"].get(0.01, 0) / n if n else 0
        fp = cell["fp"].get(0.01, 0) / n if n else 0
        print(f"    true_minor={float(rung)*100:4.0f}%  n={n:3d}  detected={det:5.1%}  false_pos={fp:5.1%}")


if __name__ == "__main__":
    main()
