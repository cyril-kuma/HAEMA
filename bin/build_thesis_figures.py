#!/usr/bin/env python3
"""
HÆMA publication figures (thesis / manuscript / conference).

Rebuilds the host-preference figures from the run endpoint tables using a single,
fixed, colour-blind-safe host colour key and aggregated, statistically-framed views
(Wilson 95% CIs) instead of illegible per-sample bars.

Run inside haema-figures:0.3.0. Two style modes: manuscript (compact, 300/600 dpi
PDF+PNG) and presentation (large fonts for slides).
"""
import argparse, csv, math, os
from collections import defaultdict, Counter

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ------------------------------------------------------------------ colours ---
# Fixed host colour key (Okabe-Ito + Paul Tol), used identically in EVERY figure.
HOST_COLORS = {
    "Homo sapiens":            "#0072B2",  # blue  (focal host)
    "Bos taurus":              "#D55E00",  # vermillion
    "Ovis aries":              "#009E73",  # bluish green
    "Capra hircus":            "#E69F00",  # orange
    "Canis lupus familiaris":  "#CC79A7",  # reddish purple
    "Sus scrofa":              "#56B4E9",  # sky blue
    "Equus asinus":            "#F0E442",  # yellow
    "Gallus gallus":           "#882255",  # wine
    "Felis catus":             "#999933",  # olive
    "unassigned":              "#BBBBBB",
    "unresolved":              "#BBBBBB",
}
COMMON = {
    "Homo sapiens": "Human", "Bos taurus": "Cattle", "Ovis aries": "Sheep",
    "Capra hircus": "Goat", "Canis lupus familiaris": "Dog", "Sus scrofa": "Pig",
    "Equus asinus": "Donkey", "Gallus gallus": "Chicken", "Felis catus": "Cat",
}
# Feeding-type / categorical accents (colour-blind safe)
C_HUMAN_ONLY = "#0072B2"
C_MIXED      = "#9467bd"   # purple — mixed human+animal
C_ANIMAL_ONLY= "#D55E00"
C_UNID       = "#BBBBBB"
C_SINGLE     = "#7f9bb3"
C_HA         = "#5e3c99"   # human-animal mix (transmission-relevant)
C_AA         = "#1b7837"   # animal-animal mix


def common(sci):
    return COMMON.get(sci, sci)


def two_line(sci):
    """Common name + italic scientific (mathtext) for axis ticks."""
    c = COMMON.get(sci)
    if not c:
        return sci
    safe = sci.replace(" ", r"\ ")
    return f"{c}\n" + r"$\it{" + safe + r"}$"


def host_color(sci):
    return HOST_COLORS.get(sci, "#444444")


# ------------------------------------------------------------------- style ---
def set_style(mode):
    base = {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.8, "axes.titleweight": "bold",
        "xtick.direction": "out", "ytick.direction": "out",
        "legend.frameon": False, "figure.dpi": 120,
        "savefig.bbox": "tight", "savefig.facecolor": "white",
        "axes.grid": False,
    }
    if mode == "presentation":
        base.update({
            "font.size": 15, "axes.titlesize": 18, "axes.labelsize": 16,
            "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 13,
            "axes.linewidth": 1.2, "lines.linewidth": 2.2,
        })
    else:
        base.update({
            "font.size": 8, "axes.titlesize": 9.5, "axes.labelsize": 8.5,
            "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
        })
    matplotlib.rcParams.update(base)


def panel_label(ax, s, dx=-0.16, dy=1.04):
    ax.text(dx, dy, s, transform=ax.transAxes, fontweight="bold",
            va="top", ha="left",
            fontsize=matplotlib.rcParams["axes.titlesize"] + 1.5)


def save(fig, outdir, name, dpi):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), dpi=dpi)
    plt.close(fig)
    print(f"  wrote {name}.pdf / .png")


# -------------------------------------------------------------------- data ---
def load_eco(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    for c in ("value", "ci_low", "ci_high", "n"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def eco_get(df, stype, metric_prefix=None, metric=None):
    sub = df[df["stratum_type"] == stype]
    if metric is not None:
        sub = sub[sub["metric"] == metric]
    if metric_prefix is not None:
        sub = sub[sub["metric"].str.startswith(metric_prefix)]
    return sub


BAD_ZONES = {"", "#N/A", "NA", "N/A", "Unknown", "nan", "None"}


def clean_zone(v):
    v = (v or "").strip()
    return "Unknown" if v in BAD_ZONES else v


def sample_meta(samplesheet, zone_column="collection_region"):
    """sample_uid -> {zone, species} from the samplesheet (zone from zone_column)."""
    zone, species = {}, {}
    with open(samplesheet) as fh:
        r = csv.DictReader(fh)
        for row in r:
            uid = f"{row['run_id']}__{row['barcode_id']}__{row['sample_id']}"
            zone[uid] = clean_zone(row.get(zone_column, ""))
            species[uid] = row.get("sibling_species", "") or "Unknown"
    return zone, species


def per_sample_hosts(host_calls):
    """field sample_uid -> {host: summed reads} (assigned hosts only)."""
    agg = defaultdict(lambda: defaultdict(float))
    with open(host_calls) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("control_status") != "sample":
                continue
            h = row.get("host_assignment", "")
            if not h or h.lower() == "unassigned":
                continue
            try:
                reads = float(row.get("host_reads") or 0)
            except ValueError:
                reads = 0.0
            agg[row["sample_uid"]][h] += reads
    return agg


def shannon(counts):
    tot = sum(counts)
    if tot <= 0:
        return 0.0
    ps = [c / tot for c in counts if c > 0]
    return -sum(p * math.log(p) for p in ps)


# ----------------------------------------------------------------- figures ---
def fig05_host_composition(eco, outdir, dpi, figsize):
    hbi = eco_get(eco, "overall", metric_prefix="host_blood_index::").copy()
    hbi["host"] = hbi["metric"].str.split("::").str[1]
    hbi = hbi.sort_values("value", ascending=True)
    n_overall = int(eco_get(eco, "overall", metric="n_identified")["value"].iloc[0])

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.45)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    # Panel A: host blood indices with Wilson CI
    y = np.arange(len(hbi))
    vals = hbi["value"].values
    lo = vals - hbi["ci_low"].values
    hi = hbi["ci_high"].values - vals
    cols = [host_color(h) for h in hbi["host"]]
    axA.barh(y, vals, color=cols, edgecolor="black", linewidth=0.4, height=0.72)
    axA.errorbar(vals, y, xerr=[lo, hi], fmt="none", ecolor="black",
                 elinewidth=0.8, capsize=2)
    axA.set_yticks(y)
    axA.set_yticklabels([two_line(h) for h in hbi["host"]])
    axA.set_xlabel("Host blood index\n(proportion of blood meals containing host)")
    axA.set_xlim(0, 1.0)
    for yi, v in zip(y, vals):
        txt = "<0.01" if v < 0.005 else f"{v:.2f}"
        axA.text(min(v + max(hi) * 0.4 + 0.02, 0.97), yi, txt,
                 va="center", ha="left", fontsize=matplotlib.rcParams["ytick.labelsize"])
    axA.set_title("Host blood indices")
    panel_label(axA, "A", dx=-0.42)

    # Panel B: feeding-type composition (stacked, among host-identified)
    g = lambda m: float(eco_get(eco, "overall", metric=m)["value"].iloc[0])
    parts = [("Human only", g("human_only_fraction"), C_HUMAN_ONLY),
             ("Mixed human + animal", g("mixed_human_animal_fraction"), C_MIXED),
             ("Animal only", g("animal_only_fraction"), C_ANIMAL_ONLY)]
    left = 0.0
    for label, frac, c in parts:
        axB.barh(0, frac, left=left, color=c, edgecolor="white", linewidth=0.8,
                 height=0.62, label=f"{label} ({frac*100:.0f}%)")
        if frac > 0.04:
            axB.text(left + frac / 2, 0, f"{frac*100:.0f}%", va="center",
                     ha="center", color="white", fontweight="bold",
                     fontsize=matplotlib.rcParams["ytick.labelsize"])
        left += frac
    axB.set_xlim(0, 1)
    axB.set_ylim(-0.5, 1.05)
    axB.set_yticks([])
    axB.set_xlabel("Proportion of host-identified meals")
    axB.legend(loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=1,
               handlelength=1.1, borderaxespad=0)
    axB.spines["left"].set_visible(False)
    axB.set_title("Feeding type")
    panel_label(axB, "B", dx=-0.08)

    fig.suptitle(f"Blood-meal host composition  (n = {n_overall} host-identified field samples)",
                 fontweight="bold", y=1.02,
                 fontsize=matplotlib.rcParams["axes.titlesize"] + 1)
    save(fig, outdir, "figure_05_host_composition", dpi)


def fig06_mixed_feeding(eco, host_calls, samplesheet, outdir, dpi, figsize):
    region, _ = sample_meta(samplesheet)
    agg = per_sample_hosts(host_calls)
    nhosts = Counter(len(v) for v in agg.values())
    combos = Counter(tuple(sorted(v)) for v in agg.values() if len(v) >= 2)

    fig = plt.figure(figsize=figsize)
    # Wide inter-panel gap so Panel B's long combination y-labels never reach into Panel A.
    gs = fig.add_gridspec(1, 2, width_ratios=[0.7, 1.5], wspace=1.0)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    # Panel A: number of distinct hosts per meal
    ks = sorted(nhosts)
    counts = [nhosts[k] for k in ks]
    total = sum(counts)
    barcols = [C_SINGLE if k == 1 else C_HA for k in ks]
    bars = axA.bar([str(k) for k in ks], counts, color=barcols,
                   edgecolor="black", linewidth=0.5, width=0.7)
    for b, c in zip(bars, counts):
        axA.text(b.get_x() + b.get_width() / 2, c, f"{c}\n({c/total*100:.0f}%)",
                 ha="center", va="bottom",
                 fontsize=matplotlib.rcParams["ytick.labelsize"])
    axA.set_xlabel("Distinct host taxa per blood meal")
    axA.set_ylabel("Field samples")
    axA.set_ylim(0, max(counts) * 1.18)
    mixed_pct = (total - nhosts[1]) / total * 100
    axA.set_title(f"Feeding multiplicity\n({mixed_pct:.0f}% mixed)")
    panel_label(axA, "A", dx=-0.30)

    # Panel B: top mixed-host combinations
    top = combos.most_common(10)[::-1]
    labels = [" + ".join(common(x) for x in combo) for combo, _ in top]
    vals = [c for _, c in top]
    def is_aa(combo):  # animal-animal (no human)
        return "Homo sapiens" not in combo
    cols = [C_AA if is_aa(combo) else C_HA for combo, _ in top]
    yy = np.arange(len(top))
    axB.barh(yy, vals, color=cols, edgecolor="black", linewidth=0.4, height=0.74)
    for yi, v in zip(yy, vals):
        axB.text(v + 0.3, yi, str(v), va="center", ha="left",
                 fontsize=matplotlib.rcParams["ytick.labelsize"])
    axB.set_yticks(yy)
    axB.set_yticklabels(labels)
    axB.set_xlabel("Number of mixed blood meals")
    axB.set_xlim(0, max(vals) * 1.15)
    axB.legend(handles=[Patch(facecolor=C_HA, label="Human + animal"),
                        Patch(facecolor=C_AA, label="Animal + animal")],
               loc="lower right")
    axB.set_title("Mixed-host blood-meal combinations")
    panel_label(axB, "B", dx=-0.42)
    save(fig, outdir, "figure_06_mixed_feeding", dpi)


def fig11_by_zone(eco, outdir, dpi, figsize, zone_label="ecological zone", fname="figure_11_feeding_by_zone"):
    zt = eco_get(eco, "ecological_zone")
    hbi = zt[zt["metric"] == "human_blood_index"].set_index("stratum")
    mix = zt[zt["metric"] == "mixed_feeding_rate"].set_index("stratum")
    # n_identified is the denominator of the indices/CIs plotted here (matches figure_09 convention)
    ntab = zt[zt["metric"] == "n_identified"].set_index("stratum")["value"]
    order = [z for z in hbi["value"].sort_values(ascending=False).index.tolist()
             if z not in BAD_ZONES]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=figsize)

    def zbar(ax, tab, color, xlabel, title, pl):
        x = np.arange(len(order))
        vals = tab.loc[order, "value"].values
        lo = vals - tab.loc[order, "ci_low"].values
        hi = tab.loc[order, "ci_high"].values - vals
        ax.bar(x, vals, color=color, edgecolor="black", linewidth=0.5, width=0.72)
        ax.errorbar(x, vals, yerr=[lo, hi], fmt="none", ecolor="black",
                    elinewidth=0.8, capsize=2.5)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{z} (n={int(ntab[z])})" for z in order],
                           rotation=40, ha="right")
        ax.set_ylabel(xlabel)
        ax.set_ylim(0, 1.05)
        ax.set_title(title)
        panel_label(ax, pl, dx=-0.20)

    zbar(axA, hbi, C_HUMAN_ONLY, "Human blood index", "Human feeding by zone", "A")
    zbar(axB, mix, C_MIXED, "Mixed-feeding rate", "Mixed feeding by zone", "B")
    fig.suptitle(f"Spatial variation in host feeding across {zone_label}s",
                 fontweight="bold", y=1.03,
                 fontsize=matplotlib.rcParams["axes.titlesize"] + 1)
    fig.subplots_adjust(wspace=0.32)
    save(fig, outdir, fname, dpi)


def fig12_alpha_diversity(host_calls, samplesheet, outdir, dpi, figsize,
                          zone_label="ecological zone", zone_column="collection_region",
                          fname="figure_12_alpha_diversity"):
    zone, _ = sample_meta(samplesheet, zone_column)
    agg = per_sample_hosts(host_calls)
    rows = []
    for uid, hosts in agg.items():
        rows.append({"zone": zone.get(uid, "Unknown"),
                     "richness": len(hosts),
                     "shannon": shannon(list(hosts.values()))})
    d = pd.DataFrame(rows)
    d = d[~d["zone"].isin(BAD_ZONES)]
    order = d["zone"].value_counts().index.tolist()  # by sample size
    ncol = d["zone"].value_counts()

    fig, axes = plt.subplots(1, 2, figsize=figsize, sharex=True)
    import seaborn as sns
    for ax, metric, title, ylab in (
        (axes[0], "richness", "Host richness", "Distinct host taxa per meal"),
        (axes[1], "shannon", "Shannon diversity", "Shannon index (H')")):
        sns.boxplot(data=d, x="zone", y=metric, order=order, ax=ax,
                    color="#cfd8e2", width=0.6, fliersize=0,
                    linewidth=0.8, showcaps=True)
        sns.stripplot(data=d, x="zone", y=metric, order=order, ax=ax,
                      color="#26456e", alpha=0.45, size=2.4, jitter=0.22)
        ax.set_title(title)
        ax.set_ylabel(ylab)
        ax.set_xlabel("")
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([f"{z} (n={int(ncol[z])})" for z in order],
                           rotation=40, ha="right")
    fig.suptitle(f"Blood-meal host diversity by {zone_label} (host-identified field samples)",
                 fontweight="bold", y=1.02,
                 fontsize=matplotlib.rcParams["axes.titlesize"] + 1)
    fig.subplots_adjust(wspace=0.28)
    save(fig, outdir, fname, dpi)


def fig13_by_species(eco, outdir, dpi, figsize):
    st = eco_get(eco, "sibling_species")
    sp_order = ["An_coluzzii", "An_gambiae_s.s", "An_arabiensis"]
    sp_order = [s for s in sp_order if s in set(st["stratum"])]
    sp_disp = {"An_coluzzii": "An. coluzzii", "An_gambiae_s.s": "An. gambiae s.s.",
               "An_arabiensis": "An. arabiensis"}
    metrics = [("human_blood_index", "Human blood index", C_HUMAN_ONLY),
               ("animal_blood_index_zoophily", "Animal blood index\n(zoophily)", C_ANIMAL_ONLY),
               ("mixed_feeding_rate", "Mixed-feeding rate", C_MIXED)]
    ntab = st[st["metric"] == "n_identified"].set_index("stratum")["value"]

    fig, ax = plt.subplots(figsize=figsize)
    nsp = len(sp_order)
    group_w = 0.8
    bw = group_w / nsp
    x = np.arange(len(metrics))
    for i, sp in enumerate(sp_order):
        vals, los, his = [], [], []
        for m, _, _ in metrics:
            r = st[(st["stratum"] == sp) & (st["metric"] == m)]
            v = r["value"].iloc[0]
            vals.append(v)
            los.append(v - r["ci_low"].iloc[0])
            his.append(r["ci_high"].iloc[0] - v)
        off = (i - (nsp - 1) / 2) * bw
        grey = ["#4477AA", "#EE6677", "#228833"][i % 3]
        disp = sp_disp.get(sp, sp)
        ital = r"$\it{" + disp.replace(" ", r"\ ") + r"}$"
        ax.bar(x + off, vals, width=bw * 0.92, color=grey, edgecolor="black",
               linewidth=0.5, label=f"{ital} (n={int(ntab[sp])})")
        ax.errorbar(x + off, vals, yerr=[los, his], fmt="none", ecolor="black",
                    elinewidth=0.7, capsize=2)
    ax.set_xticks(x)
    ax.set_xticklabels([lbl for _, lbl, _ in metrics])
    ax.set_ylabel("Index value (proportion of meals)")
    ax.set_ylim(0, 1.0)
    ax.legend(title="Vector sibling species", loc="upper right")
    ax.set_title("Host-feeding behaviour by vector sibling species")
    save(fig, outdir, "figure_13_feeding_by_species", dpi)


def fig14_host_zone_heatmap(eco, outdir, dpi, figsize, zone_label="ecological zone",
                            fname="figure_14_host_zone_heatmap"):
    zt = eco_get(eco, "ecological_zone")
    hb = zt[zt["metric"].str.startswith("host_blood_index::")].copy()
    hb["host"] = hb["metric"].str.split("::").str[1]
    mat = hb.pivot_table(index="host", columns="stratum", values="value", aggfunc="first")
    # order rows by overall blood index, cols by n
    overall = eco_get(eco, "overall", metric_prefix="host_blood_index::").copy()
    overall["host"] = overall["metric"].str.split("::").str[1]
    row_order = overall.sort_values("value", ascending=False)["host"].tolist()
    row_order = [h for h in row_order if h in mat.index]
    # n_identified is the denominator of the indices/CIs plotted here (matches figure_09 convention)
    ntab = zt[zt["metric"] == "n_identified"].set_index("stratum")["value"]
    col_order = ntab.sort_values(ascending=False).index.tolist()
    col_order = [c for c in col_order if c in mat.columns and c not in BAD_ZONES]
    mat = mat.reindex(index=row_order, columns=col_order).fillna(0.0)

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(mat.values, cmap="viridis", aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(col_order)))
    ax.set_xticklabels([f"{z}\n(n={int(ntab[z])})" for z in col_order],
                       rotation=35, ha="right")
    ax.set_yticks(np.arange(len(row_order)))
    ax.set_yticklabels([two_line(h) for h in row_order])
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            if v > 0.005:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        color="white" if v < 0.6 else "black",
                        fontsize=matplotlib.rcParams["ytick.labelsize"] - 0.5)
    cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("Host blood index")
    ax.set_title(f"Host blood index by {zone_label}")
    save(fig, outdir, fname, dpi)


def fig15_decontam(count_table, decontam, samplesheet, outdir, dpi, figsize):
    # decontam_results.tsv prevalence columns are empty in this run; reconstruct
    # field/neg-control prevalence from the ASV count table + control labels.
    neg = set()
    with open(samplesheet) as fh:
        for row in csv.DictReader(fh):
            if row.get("sample_type") == "negative_control":
                neg.add(f"{row['run_id']}__{row['barcode_id']}__{row['sample_id']}")
    flag = {}
    with open(decontam) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            flag[row["feature_id"]] = str(row.get("contaminant", "")).lower() in ("true", "1", "yes")

    df = pd.read_csv(count_table, sep="\t")
    meta = ["feature_id", "marker", "sequence"]
    samples = [c for c in df.columns if c not in meta]
    neg_cols = [c for c in samples if c in neg]
    field_cols = [c for c in samples if c not in neg]
    F = df[field_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    N = df[neg_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    d = pd.DataFrame({
        "feature_id": df["feature_id"],
        "pos_prevalence": (F > 0).sum(axis=1) / max(len(field_cols), 1),
        "neg_prevalence": (N > 0).sum(axis=1) / max(len(neg_cols), 1),
        "reads": F.sum(axis=1) + N.sum(axis=1),
    })
    d["contaminant"] = d["feature_id"].map(lambda x: flag.get(x, False))
    flagm = d["contaminant"].values
    size = 10 + 26 * np.log10(d["reads"].clip(lower=1).values) / \
        max(np.log10(max(d["reads"].max(), 10)), 1)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(d.loc[~flagm, "pos_prevalence"], d.loc[~flagm, "neg_prevalence"],
               s=size[~flagm], color="#4575b4", alpha=0.5, edgecolor="none",
               label=f"Retained (n={int((~flagm).sum())})")
    ax.scatter(d.loc[flagm, "pos_prevalence"], d.loc[flagm, "neg_prevalence"],
               s=np.maximum(size[flagm], 60), color="#d73027", edgecolor="black",
               linewidth=0.6, label=f"Flagged contaminant (n={int(flagm.sum())})", zorder=5)
    xmax = float(max(d["pos_prevalence"].max(), 0.05)) * 1.15
    ymax = float(max(d["neg_prevalence"].max(), 0.05)) * 1.12
    # annotate flagged features to the LEFT so labels stay clear of the legend
    for _, r in d[flagm].iterrows():
        ax.annotate(r["feature_id"], (r["pos_prevalence"], r["neg_prevalence"]),
                    xytext=(-8, 0), textcoords="offset points", ha="right", va="center",
                    fontsize=matplotlib.rcParams["ytick.labelsize"] - 0.5, color="#333333")
    ax.set_xlabel("Prevalence in field samples")
    ax.set_ylabel("Prevalence in negative controls")
    ax.set_xlim(-0.01 * xmax, xmax)
    ax.set_ylim(-0.02 * ymax, ymax)
    ax.legend(loc="center right", title="point size ∝ total reads")
    ax.set_title("Contaminant identification (decontam prevalence)")
    save(fig, outdir, "figure_15_decontam", dpi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--endpoint-dir", required=True)
    ap.add_argument("--reports-dir", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mode", choices=["manuscript", "presentation"], default="manuscript")
    ap.add_argument("--eco-indices", default="",
                    help="override ecological_indices.tsv (e.g. a bioclimatic-zone-stratified table)")
    ap.add_argument("--zone-column", default="collection_region",
                    help="samplesheet column for per-sample zone grouping (fig 12)")
    ap.add_argument("--zone-label", default="ecological zone",
                    help="display label for the zone stratum (titles/axes)")
    ap.add_argument("--only", default="",
                    help="comma-separated figure ids to build, e.g. 11,12,14 (default: all)")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    set_style(a.mode)
    dpi = 600 if a.mode == "manuscript" else 200
    only = {s.strip() for s in a.only.split(",") if s.strip()} or None
    want = lambda fid: only is None or fid in only

    eco = load_eco(a.eco_indices or os.path.join(a.endpoint_dir, "ecological_indices.tsv"))
    hc = os.path.join(a.endpoint_dir, "host_call_table.tsv")
    cnt = os.path.join(a.endpoint_dir, "asv_count_table.tsv")
    dec = os.path.join(a.reports_dir, "decontam_results.tsv")

    # figure sizes (inches) per mode
    if a.mode == "manuscript":
        sz = dict(f5=(7.2, 3.2), f6=(7.8, 3.4), f11=(7.0, 3.4), f12=(7.0, 3.4),
                  f13=(5.2, 3.6), f14=(6.6, 4.0), f15=(5.4, 4.0))
    else:
        sz = dict(f5=(12, 5.6), f6=(13.5, 6), f11=(12, 6), f12=(12, 6),
                  f13=(9.5, 6.2), f14=(11, 6.6), f15=(9.5, 6.4))

    zl = a.zone_label
    print(f"[{a.mode}] writing to {a.outdir}  (zone='{zl}')")
    if want("5"):  fig05_host_composition(eco, a.outdir, dpi, sz["f5"])
    if want("6"):  fig06_mixed_feeding(eco, hc, a.samplesheet, a.outdir, dpi, sz["f6"])
    if want("11"): fig11_by_zone(eco, a.outdir, dpi, sz["f11"], zone_label=zl)
    if want("12"): fig12_alpha_diversity(hc, a.samplesheet, a.outdir, dpi, sz["f12"],
                                         zone_label=zl, zone_column=a.zone_column)
    if want("13"): fig13_by_species(eco, a.outdir, dpi, sz["f13"])
    if want("14"): fig14_host_zone_heatmap(eco, a.outdir, dpi, sz["f14"], zone_label=zl)
    if want("15"): fig15_decontam(cnt, dec, a.samplesheet, a.outdir, dpi, sz["f15"])
    print("done.")


if __name__ == "__main__":
    main()
