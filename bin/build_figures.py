#!/usr/bin/env python3
"""Generate publication-quality manuscript figures from HÆMA endpoint tables.

Every figure is built strictly from real pipeline outputs (the TSVs under
``05_endpoint_files`` and ``06_reports`` plus ``run_manifest.json``); nothing is
simulated. Each figure is written as vector PDF + SVG and a 300 dpi PNG, and the
script emits a ``figure_captions.md`` (draft captions naming the input files) and a
``figure_manifest.tsv`` recording what was produced.

Run standalone:
    python build_figures.py --endpoint-dir RESULTS/05_endpoint_files \
        --reports-dir RESULTS/06_reports --manifest RESULTS/05_endpoint_files/run_manifest.json \
        --outdir RESULTS/07_figures

It is also driven by the BUILD_FIGURES Nextflow process.
"""
import argparse
import json
import sys
import traceback
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless / container-safe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.lines import Line2D

# --------------------------------------------------------------------------------------
# Publication style
# --------------------------------------------------------------------------------------
plt.rcParams.update({
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "pdf.fonttype": 42,   # embed editable TrueType text in PDF
    "ps.fonttype": 42,
    "svg.fonttype": "none",  # keep text as text in SVG (editable in Illustrator/Inkscape)
})

# Colour-blind-safe Okabe–Ito palette, with a stable taxon ordering.
HOST_COLORS = {
    "Homo sapiens": "#0072B2",
    "Ovis aries": "#E69F00",
    "Bos taurus": "#009E73",
    "Capra hircus": "#CC79A7",
    "Canis lupus familiaris": "#D55E00",
    "unassigned": "#BBBBBB",
}
OTHER_HOST_COLOR = "#56B4E9"
HOST_ORDER = list(HOST_COLORS.keys())

SAMPLE_TYPE_COLORS = {
    "sample": "#4C72B0",
    "positive_control": "#55A868",
    "negative_control": "#C44E52",
}
SAMPLE_TYPE_LABELS = {
    "sample": "Field sample",
    "positive_control": "Positive control",
    "negative_control": "Negative control",
}
MARKER_ORDER = ["cyt_b", "co1_short", "co1_long"]
MARKER_LABELS = {"cyt_b": "cyt $b$", "co1_short": "short COI", "co1_long": "long COI"}
MARKER_COLORS = {"cyt_b": "#8172B3", "co1_short": "#CCB974", "co1_long": "#64B5CD"}

ITALIC_TAXA = set(HOST_COLORS) - {"unassigned"}


def host_color(name):
    return HOST_COLORS.get(name, OTHER_HOST_COLOR)


def italicise_taxon(name):
    """Wrap a binomial in mathtext italics, leaving 'unassigned'/multiword labels sane."""
    if name in ITALIC_TAXA:
        return "$\\it{" + name.replace(" ", "\\ ") + "}$"
    return name


# --------------------------------------------------------------------------------------
# IO helpers
# --------------------------------------------------------------------------------------
def read_tsv(path):
    """Read a TSV into a DataFrame; return an empty frame if missing/empty."""
    p = Path(path) if path else None
    if not p or not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    except Exception as exc:  # pragma: no cover - defensive
        print(f"  [warn] could not read {p}: {exc}", file=sys.stderr)
        return pd.DataFrame()


def num(series):
    return pd.to_numeric(series, errors="coerce")


def safe_float(value):
    """Parse a scalar to float, tolerating ''/'NA'/None (returns None)."""
    try:
        if value in (None, "", "NA", "nan"):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def panel_label(ax, letter, dx=-0.07, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=12, fontweight="bold",
            va="top", ha="right")


def save_figure(fig, outdir, name, formats, manifest, caption, inputs):
    outdir = Path(outdir)
    written = []
    for fmt in formats:
        path = outdir / f"{name}.{fmt}"
        fig.savefig(path, format=fmt)
        written.append(path.name)
    plt.close(fig)
    manifest.append({
        "figure": name,
        "files": ";".join(written),
        "inputs": ";".join(inputs),
        "caption": caption.replace("\n", " ").strip(),
    })
    print(f"  [ok] {name}: {', '.join(written)}")
    return written


# --------------------------------------------------------------------------------------
# Figure 1 — pipeline workflow schematic (architecture, not data)
# --------------------------------------------------------------------------------------
def fig_workflow(outdir, formats, manifest):
    stages = [
        ("01  Input & metadata validation", "samplesheet + MIEM/MIMARKS", "#264653"),
        ("02  Merge · primer-trim · Q/length filter", "per-sample FASTQ QC", "#287271"),
        ("03  Marker split: cyt $b$ / short COI / long COI", "amplicon length + primer", "#2A9D8F"),
        ("04  Mixed-template denoising", "UMAP + HDBSCAN (greedy fallback)", "#8AB17D"),
        ("05  Cluster consensus + Medaka polish", "per-cluster representative", "#E9C46A"),
        ("06  Taxonomy: curated BLASTn + LCA", "vertebrate mitogenome panel", "#F4A261"),
        ("07  Contamination assessment", "negative-control background · decontam", "#EE8959"),
        ("08  Mixed-host evidence & host calls", "RAMBO-style abundance model", "#E76F51"),
        ("09  Endpoints, figures & reports", "master table · phyloseq · MultiQC", "#BC4749"),
    ]
    fig, ax = plt.subplots(figsize=(7.2, 9.0))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, len(stages) * 1.16 + 0.4)
    ax.axis("off")
    box_w, box_h = 7.4, 0.78
    x0 = 1.3
    centres = []
    for i, (title, sub, color) in enumerate(reversed(stages)):
        y = i * 1.16 + 0.5
        cy = y + box_h / 2
        centres.append(cy)
        box = FancyBboxPatch((x0, y), box_w, box_h, boxstyle="round,pad=0.02,rounding_size=0.12",
                             linewidth=0, facecolor=color, alpha=0.92)
        ax.add_patch(box)
        ax.text(x0 + 0.25, cy + 0.13, title, color="white", fontsize=9.5, fontweight="bold",
                va="center", ha="left")
        ax.text(x0 + 0.25, cy - 0.20, sub, color="white", fontsize=7.8, va="center", ha="left",
                alpha=0.95)
    # downward arrows between consecutive boxes
    for upper, lower in zip(centres[1:], centres[:-1]):
        ax.add_patch(FancyArrowPatch((x0 + box_w / 2, upper - box_h / 2),
                                     (x0 + box_w / 2, lower + box_h / 2),
                                     arrowstyle="-|>", mutation_scale=12, linewidth=1.1,
                                     color="#555555"))
    ax.text(x0 + box_w / 2, len(stages) * 1.16 + 0.3,
            "HÆMA — ONT blood-meal metabarcoding workflow",
            ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.text(x0 + box_w / 2, 0.12,
            "Oxford Nanopore R10.4.1 · tri-marker (cyt $b$, short COI, long COI) · containerised Nextflow DSL2",
            ha="center", va="bottom", fontsize=7.6, color="#444444", style="italic")
    caption = ("Figure 1. HÆMA pipeline architecture. Already-basecalled MinKNOW FASTQ are validated, "
               "primer-trimmed and quality/length filtered, split into three mitochondrial markers, "
               "denoised to resolve mixed templates (UMAP/HDBSCAN with a greedy fallback), reduced to "
               "per-cluster consensus sequences (optionally Medaka-polished), assigned taxonomy by "
               "curated BLASTn with a conservative LCA, screened against negative-control background, "
               "and summarised into host calls and publication endpoints.")
    return save_figure(fig, outdir, "figure_01_workflow", formats, manifest, caption,
                        ["pipeline architecture (modules/local/*)"])


# --------------------------------------------------------------------------------------
# Figure 2 — sequencing & QC
# --------------------------------------------------------------------------------------
def _qc_per_sample(qc):
    """Collapse qc_summary (7 rows/sample at mixed granularity) to tidy per-sample/per-marker views."""
    src = qc["source_file"] if "source_file" in qc else pd.Series([""] * len(qc))
    merge = qc[src.str.contains("merge_stats", na=False)].copy()
    trim = qc[src.str.contains("trim_filter_split", na=False)].copy()
    denoise = qc[src.str.contains("mixed_denoise", na=False)].copy()
    return merge, trim, denoise


def fig_qc(qc, master, outdir, formats, manifest):
    merge, trim, denoise = _qc_per_sample(qc)
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.4))
    axA, axB, axC, axD = axes.ravel()

    # ---- A: read-processing funnel (totals across all samples) ----
    raw = num(merge["reads"]).sum() if "reads" in merge else np.nan
    # pass_quality is a per-sample value repeated across the marker rows -> take one per sample
    quality = num(trim.groupby("sample_uid")["pass_quality"].first()).sum() if "pass_quality" in trim else np.nan
    assigned = num(trim["written_reads"]).sum() if "written_reads" in trim else np.nan
    retained = num(denoise["retained_reads"]).sum() if "retained_reads" in denoise else np.nan
    stages = ["Raw\nmerged", "Quality\npassed", "Marker\nassigned", "Retained\n(denoised)"]
    vals = [raw, quality, assigned, retained]
    colors = ["#264653", "#2A9D8F", "#E9C46A", "#E76F51"]
    ypos = np.arange(len(stages))[::-1]
    axA.barh(ypos, vals, color=colors, height=0.62)
    for y, v in zip(ypos, vals):
        pct = f"  {v/raw*100:.1f}%" if raw and not np.isnan(v) else ""
        axA.text(v, y, f"  {int(v):,}{pct}", va="center", ha="left", fontsize=7.6)
    axA.set_yticks(ypos)
    axA.set_yticklabels(stages, fontsize=7.8)
    axA.set_xlabel("Reads (summed over all samples)")
    axA.set_xlim(0, raw * 1.18 if raw else None)
    axA.set_title("Read-processing funnel")
    panel_label(axA, "A")

    # ---- B: per-sample mean read quality ----
    if "mean_read_q" in merge:
        q = num(merge["mean_read_q"]).dropna()
        axB.hist(q, bins=np.arange(np.floor(q.min()), np.ceil(q.max()) + 1, 1.0),
                 color="#4C72B0", edgecolor="white", linewidth=0.6)
        axB.axvline(20, color="#C44E52", linestyle="--", linewidth=1.2)
        axB.text(20.2, axB.get_ylim()[1] * 0.92, "Q20 filter", color="#C44E52", fontsize=7.5)
        axB.set_xlabel("Mean read quality (Phred)")
        axB.set_ylabel(f"Samples (n = {len(q)})")
    axB.set_title("Per-sample read quality")
    panel_label(axB, "B")

    # ---- C: consensus length by marker (real per-feature lengths) ----
    if not master.empty and "sequence" in master and "marker" in master:
        m = master[(master["sequence"] != "")].copy()
        m["seqlen"] = m["sequence"].str.len()
        data, labels, cols = [], [], []
        for mk in MARKER_ORDER:
            vals_m = m.loc[m["marker"] == mk, "seqlen"].astype(int)
            if len(vals_m):
                data.append(vals_m.values)
                labels.append(f"{MARKER_LABELS[mk]}\n(n={len(vals_m)})")
                cols.append(MARKER_COLORS[mk])
        if data:
            parts = axC.boxplot(data, vert=True, patch_artist=True, widths=0.55,
                                showfliers=False, medianprops=dict(color="black", linewidth=1.2))
            for patch, c in zip(parts["boxes"], cols):
                patch.set_facecolor(c)
                patch.set_alpha(0.85)
            for i, (vals_m, c) in enumerate(zip(data, cols), start=1):
                jitter = np.random.default_rng(i).normal(0, 0.05, size=len(vals_m))
                axC.scatter(np.full(len(vals_m), i) + jitter, vals_m, s=6, color="black",
                            alpha=0.35, zorder=3, linewidths=0)
            axC.set_xticks(range(1, len(labels) + 1))
            axC.set_xticklabels(labels, fontsize=7.6)
            axC.set_ylabel("Consensus length (bp)")
    axC.set_title("Amplicon size by marker")
    panel_label(axC, "C")

    # ---- D: reads assigned per marker, by assignment evidence ----
    if not trim.empty and "written_reads" in trim:
        by_primer, by_length = [], []
        for mk in MARKER_ORDER:
            sub = trim[trim["marker"] == mk]
            by_primer.append(num(sub["assigned_by_primer"]).sum() if "assigned_by_primer" in sub else 0)
            by_length.append(num(sub["assigned_by_length"]).sum() if "assigned_by_length" in sub else 0)
        x = np.arange(len(MARKER_ORDER))
        axD.bar(x, by_primer, color="#2A9D8F", label="by primer")
        axD.bar(x, by_length, bottom=by_primer, color="#E9C46A", label="by length")
        axD.set_xticks(x)
        axD.set_xticklabels([MARKER_LABELS[m] for m in MARKER_ORDER])
        axD.set_ylabel("Reads assigned")
        axD.legend(loc="upper right")
    axD.set_title("Marker assignment yield")
    panel_label(axD, "D")

    fig.tight_layout(w_pad=2.5, h_pad=3.0)
    caption = ("Figure 2. Sequencing and quality control. (A) Read-processing funnel summed over all "
               "libraries: raw merged reads, reads passing the Q20 mean-quality filter, reads assigned "
               "to a marker, and reads retained after mixed-template denoising. (B) Distribution of "
               "per-sample mean read quality; dashed line is the Q20 threshold. (C) Consensus-sequence "
               "length per marker, recovering the three expected amplicon size classes (points are "
               "individual consensus features). (D) Reads assigned to each marker, partitioned by "
               "primer-based versus length-based assignment. Inputs: qc_summary.tsv, "
               "bloodmeal_master_endpoint.tsv.")
    return save_figure(fig, outdir, "figure_02_sequencing_qc", formats, manifest, caption,
                        ["qc_summary.tsv", "bloodmeal_master_endpoint.tsv"])


# --------------------------------------------------------------------------------------
# Figure 3 — sample depth & denoising
# --------------------------------------------------------------------------------------
def fig_depth_denoise(sample_sum, qc, outdir, formats, manifest):
    fig = plt.figure(figsize=(9.4, 7.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], hspace=0.42, wspace=0.30)
    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1, 0])
    axC = fig.add_subplot(gs[1, 1])

    # ---- A: per-sample retained read depth, coloured by sample type ----
    if not sample_sum.empty:
        d = sample_sum.copy()
        d["depth"] = num(d["total_asv_reads"])
        d = d.sort_values("depth", ascending=False).reset_index(drop=True)
        colors = [SAMPLE_TYPE_COLORS.get(t, "#999999") for t in d["sample_type"]]
        axA.bar(np.arange(len(d)), d["depth"].clip(lower=0.6), color=colors, width=0.82)
        axA.set_yscale("log")
        axA.set_ylabel("Retained reads (log)")
        axA.set_xlim(-0.8, len(d) - 0.2)
        axA.set_xticks(np.arange(len(d)))
        axA.set_xticklabels(d["sample_id"], rotation=90, fontsize=5.6)
        axA.set_xlabel("Sample (sorted by depth)")
        handles = [plt.Rectangle((0, 0), 1, 1, color=SAMPLE_TYPE_COLORS[t])
                   for t in SAMPLE_TYPE_COLORS]
        axA.legend(handles, [SAMPLE_TYPE_LABELS[t] for t in SAMPLE_TYPE_COLORS],
                   loc="upper right", ncol=3)
        axA.set_title("Per-sample retained sequencing depth")
    panel_label(axA, "A", dx=-0.04, dy=1.10)

    # ---- B: retained ASVs vs assigned hosts ----
    if not sample_sum.empty:
        x = num(sample_sum["retained_asvs"])
        y = num(sample_sum["assigned_hosts"])
        for t in SAMPLE_TYPE_COLORS:
            m = sample_sum["sample_type"] == t
            axB.scatter(x[m], y[m], s=34, color=SAMPLE_TYPE_COLORS[t], alpha=0.85,
                        edgecolor="white", linewidth=0.5, label=SAMPLE_TYPE_LABELS[t])
        lim = max(x.max(), y.max()) * 1.08 if len(x) else 1
        axB.plot([0, lim], [0, lim], color="#999999", linestyle=":", linewidth=0.9)
        axB.set_xlabel("Retained ASVs / clusters")
        axB.set_ylabel("Assigned host features")
        axB.set_title("Clustering vs taxonomic yield")
        axB.legend(loc="lower right", fontsize=6.8)
    panel_label(axB, "B")

    # ---- C: denoising backend usage per marker ----
    _, _, denoise = _qc_per_sample(qc)
    if not denoise.empty and "backend_used" in denoise:
        backends = ["umap_hdbscan", "greedy_identity"]
        bcolors = {"umap_hdbscan": "#2A9D8F", "greedy_identity": "#E9C46A"}
        counts = {b: [] for b in backends}
        for mk in MARKER_ORDER:
            sub = denoise[denoise["marker"] == mk]
            for b in backends:
                counts[b].append(int((sub["backend_used"] == b).sum()))
        x = np.arange(len(MARKER_ORDER))
        bottom = np.zeros(len(MARKER_ORDER))
        for b in backends:
            axC.bar(x, counts[b], bottom=bottom, color=bcolors[b],
                    label=b.replace("_", " "))
            bottom += np.array(counts[b])
        axC.set_xticks(x)
        axC.set_xticklabels([MARKER_LABELS[m] for m in MARKER_ORDER])
        axC.set_ylabel("Sample × marker groups")
        axC.set_title("Denoising backend used")
        axC.legend(loc="upper right", fontsize=6.8)
    panel_label(axC, "C")

    caption = ("Figure 3. Library performance and denoising. (A) Retained read depth per sample on a "
               "log scale, coloured by control class; negative controls sit far below field samples. "
               "(B) Retained ASV/cluster count versus number of assigned host features per sample "
               "(dotted line = parity). (C) Number of sample × marker groups resolved by the UMAP/"
               "HDBSCAN backend versus the deterministic greedy-identity fallback (used for low-read "
               "groups). Inputs: sample_level_summary.tsv, qc_summary.tsv.")
    return save_figure(fig, outdir, "figure_03_depth_denoising", formats, manifest, caption,
                        ["sample_level_summary.tsv", "qc_summary.tsv"])


# --------------------------------------------------------------------------------------
# Figure 4 — host assignment overview
# --------------------------------------------------------------------------------------
def fig_host_overview(host_calls, master, outdir, formats, manifest):
    fig = plt.figure(figsize=(9.6, 7.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], hspace=0.45, wspace=0.42)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, :])

    field = host_calls[host_calls["control_status"] == "sample"] if not host_calls.empty else host_calls
    real = field[~field["host_assignment"].isin(["", "unassigned"])] if not field.empty else field

    # ---- A: host taxon frequency (number of field samples) ----
    if not real.empty:
        freq = real.groupby("host_assignment")["sample_uid"].nunique().sort_values()
        ypos = np.arange(len(freq))
        axA.barh(ypos, freq.values, color=[host_color(h) for h in freq.index])
        axA.set_yticks(ypos)
        axA.set_yticklabels([italicise_taxon(h) for h in freq.index])
        for y, v in zip(ypos, freq.values):
            axA.text(v + 0.1, y, str(int(v)), va="center", fontsize=7.5)
        axA.set_xlabel("Field samples with host detected")
        axA.set_title("Detected blood-meal hosts")
        axA.set_xlim(0, freq.max() * 1.15)
    panel_label(axA, "A")

    # ---- B: host × marker detection heatmap ----
    if not real.empty:
        hosts = [h for h in HOST_ORDER if h in set(real["host_assignment"])]
        hosts += sorted(set(real["host_assignment"]) - set(hosts) - {"unassigned"})
        mat = np.zeros((len(hosts), len(MARKER_ORDER)))
        for i, h in enumerate(hosts):
            for j, mk in enumerate(MARKER_ORDER):
                mat[i, j] = real[(real["host_assignment"] == h) & (real["marker"] == mk)]["sample_uid"].nunique()
        im = axB.imshow(mat, cmap="YlGnBu", aspect="auto")
        axB.set_xticks(range(len(MARKER_ORDER)))
        axB.set_xticklabels([MARKER_LABELS[m] for m in MARKER_ORDER])
        axB.set_yticks(range(len(hosts)))
        axB.set_yticklabels([italicise_taxon(h) for h in hosts])
        for i in range(len(hosts)):
            for j in range(len(MARKER_ORDER)):
                if mat[i, j] > 0:
                    axB.text(j, i, int(mat[i, j]), ha="center", va="center",
                             color="white" if mat[i, j] > mat.max() * 0.6 else "black", fontsize=7.5)
        cbar = fig.colorbar(im, ax=axB, fraction=0.046, pad=0.04)
        cbar.set_label("samples", fontsize=7)
        axB.set_title("Host recovery by marker")
    panel_label(axB, "B")

    # ---- C: BLAST identity vs coverage of accepted assignments ----
    if not master.empty and "pident" in master and "coverage" in master:
        m = master.copy()
        m["pid"] = num(m["pident"])
        m["cov"] = num(m["coverage"])
        m = m[(~m["pid"].isna()) & (~m["cov"].isna()) &
              (~m["host_assignment"].isin(["", "unassigned"]))]
        for h in HOST_ORDER + sorted(set(m["host_assignment"]) - set(HOST_ORDER)):
            sub = m[m["host_assignment"] == h]
            if sub.empty:
                continue
            axC.scatter(sub["pid"], sub["cov"], s=26, color=host_color(h), alpha=0.8,
                        edgecolor="white", linewidth=0.4, label=italicise_taxon(h))
        axC.set_xlim(96.5, 100.4)
        axC.set_ylim(78, 101.8)
        axC.axvline(97, color="#C44E52", linestyle="--", linewidth=1.0)
        axC.axhline(80, color="#C44E52", linestyle="--", linewidth=1.0)
        axC.text(96.93, 101.6, "min identity 97% ", color="#C44E52", fontsize=7,
                 rotation=90, va="top", ha="center")
        axC.text(100.35, 80.3, "min coverage 80%", color="#C44E52", fontsize=7,
                 va="bottom", ha="right")
        axC.set_xlabel("BLAST identity (%)")
        axC.set_ylabel("Query coverage (%)")
        axC.set_title("Assignment confidence (accepted features)")
        axC.legend(loc="center", ncol=3, fontsize=6.6, columnspacing=0.9, handletextpad=0.3)
    panel_label(axC, "C", dx=-0.05, dy=1.08)

    caption = ("Figure 4. Host assignment overview. (A) Number of field samples in which each "
               "vertebrate host was detected. (B) Host recovery by marker (cell = number of samples "
               "with that host called on that marker), showing the complementarity of the three loci. "
               "(C) BLAST percent identity versus query coverage for every accepted consensus feature, "
               "coloured by assigned host; dashed lines mark the 97% identity / 80% coverage acceptance "
               "thresholds. Inputs: host_call_table.tsv, bloodmeal_master_endpoint.tsv.")
    return save_figure(fig, outdir, "figure_04_host_assignment", formats, manifest, caption,
                        ["host_call_table.tsv", "bloodmeal_master_endpoint.tsv"])


# --------------------------------------------------------------------------------------
# Figure 5 — blood-meal host composition per sample
# --------------------------------------------------------------------------------------
def fig_host_composition(host_calls, outdir, formats, manifest):
    field = host_calls[(host_calls["control_status"] == "sample")].copy() if not host_calls.empty else host_calls
    field = field[~field["host_assignment"].isin(["", "unassigned"])] if not field.empty else field
    if field.empty:
        print("  [skip] figure 05: no field-sample host calls")
        return []
    field["frac"] = num(field["host_fraction"])
    present_markers = [m for m in MARKER_ORDER if m in set(field["marker"])]
    hosts_present = [h for h in HOST_ORDER if h in set(field["host_assignment"])]
    hosts_present += sorted(set(field["host_assignment"]) - set(hosts_present))

    fig, axes = plt.subplots(len(present_markers), 1,
                             figsize=(9.6, 2.5 * len(present_markers) + 0.6), squeeze=False)
    for ax, mk in zip(axes.ravel(), present_markers):
        sub = field[field["marker"] == mk]
        samples = sorted(sub["sample_id"].unique())
        x = np.arange(len(samples))
        bottoms = np.zeros(len(samples))
        for h in hosts_present:
            vals = []
            for s in samples:
                r = sub[(sub["sample_id"] == s) & (sub["host_assignment"] == h)]
                vals.append(num(r["host_fraction"]).sum() if not r.empty else 0.0)
            vals = np.array(vals)
            if vals.sum() > 0:
                ax.bar(x, vals, bottom=bottoms, color=host_color(h), width=0.82,
                       edgecolor="white", linewidth=0.4)
                bottoms += vals
        # fill each bar to 1.0 with an "unresolved" (noise/unassigned) remainder so bars are
        # true compositions and low-retention single-host bars do not appear to float
        remainder = np.clip(1.0 - bottoms, 0.0, None)
        if remainder.sum() > 1e-9:
            ax.bar(x, remainder, bottom=bottoms, color="#E5E5E5", width=0.82,
                   edgecolor="white", linewidth=0.4)
        # mark mixed-host samples with a caret
        for xi, s in enumerate(samples):
            n_hosts = sub[sub["sample_id"] == s]["host_assignment"].nunique()
            if n_hosts > 1:
                ax.text(xi, 1.02, "▲", ha="center", va="bottom", color="#222222", fontsize=7)
        ax.set_xticks(x)
        ax.set_xticklabels(samples, rotation=90, fontsize=6)
        ax.set_ylim(0, 1.12)
        ax.set_ylabel("Host read fraction")
        ax.set_title(f"{MARKER_LABELS[mk]}  (n = {len(samples)} samples)", loc="left")
        ax.margins(x=0.01)
    # shared legend
    handles = [plt.Rectangle((0, 0), 1, 1, color=host_color(h)) for h in hosts_present]
    labels = [italicise_taxon(h) for h in hosts_present]
    handles.append(plt.Rectangle((0, 0), 1, 1, color="#E5E5E5"))
    labels.append("unresolved")
    handles.append(Line2D([0], [0], marker="^", color="none", markerfacecolor="#222",
                          markeredgecolor="#222", markersize=6))
    labels.append("mixed-host feed")
    axes.ravel()[0].legend(handles, labels, loc="upper left", bbox_to_anchor=(1.005, 1.0),
                           fontsize=7, title="Host")
    fig.suptitle("Blood-meal host composition per field sample", fontsize=11, fontweight="bold",
                 x=0.07, ha="left", y=0.995)
    fig.tight_layout(rect=[0, 0, 0.86, 0.98])
    caption = ("Figure 5. Blood-meal host composition. Stacked bars give the read fraction of each "
               "vertebrate host within every field sample, shown separately for each marker. Bars that "
               "stack more than one colour (▲) are mixed-host blood meals. Host fractions are abundance "
               "evidence summaries, not validated quantitative estimates (see Limitations). Input: "
               "host_call_table.tsv.")
    return save_figure(fig, outdir, "figure_05_host_composition", formats, manifest, caption,
                        ["host_call_table.tsv"])


# --------------------------------------------------------------------------------------
# Figure 6 — mixed-host feeding evidence
# --------------------------------------------------------------------------------------
def fig_mixed_host(host_calls, rambo, outdir, formats, manifest):
    fig = plt.figure(figsize=(9.4, 4.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.9, 1.4], wspace=0.32)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    field = host_calls[host_calls["control_status"] == "sample"] if not host_calls.empty else host_calls

    # ---- A: feeding-type composition (sample × marker groups) ----
    if not field.empty and "mixed_status" in field:
        order = ["single_host", "mixed_host", "no_host_signal"]
        labels = {"single_host": "single-host", "mixed_host": "mixed-host", "no_host_signal": "no signal"}
        cols = {"single_host": "#4C72B0", "mixed_host": "#E76F51", "no_host_signal": "#BBBBBB"}
        groups = field.drop_duplicates(["sample_uid", "marker"])
        counts = [int((groups["mixed_status"] == s).sum()) for s in order]
        wedges, _, autotexts = axA.pie(counts, colors=[cols[s] for s in order],
                                       autopct=lambda p: f"{p*sum(counts)/100:.0f}", startangle=90,
                                       wedgeprops=dict(width=0.42, edgecolor="white"))
        for t in autotexts:
            t.set_fontsize(8)
        axA.legend(wedges, [labels[s] for s in order], loc="center", fontsize=7.5,
                   frameon=False, bbox_to_anchor=(0.5, -0.08), ncol=1)
        axA.set_title("Feeding type\n(sample × marker groups)")
    panel_label(axA, "A", dx=0.02, dy=1.10)

    # ---- B: detected mixed feeds, host fractions ----
    if not field.empty:
        mixed = field[field["mixed_status"] == "mixed_host"].copy()
        if not mixed.empty:
            mixed["frac"] = num(mixed["host_fraction"])
            mixed["key"] = mixed["sample_id"] + "  ·  " + mixed["marker"].map(
                lambda m: MARKER_LABELS.get(m, m))
            keys = list(dict.fromkeys(mixed["key"]))
            ypos = np.arange(len(keys))[::-1]
            for k, y in zip(keys, ypos):
                sub = mixed[mixed["key"] == k].sort_values("frac", ascending=False)
                left = 0.0
                for _, r in sub.iterrows():
                    axB.barh(y, r["frac"], left=left, color=host_color(r["host_assignment"]),
                             edgecolor="white", linewidth=0.6)
                    if r["frac"] > 0.06:
                        axB.text(left + r["frac"] / 2, y,
                                 f"{italicise_taxon(r['host_assignment'])}\n{r['frac']*100:.0f}%",
                                 ha="center", va="center", fontsize=6.6,
                                 color="white" if r["host_assignment"] != "unassigned" else "black")
                    left += r["frac"]
                if left < 0.999:  # unresolved (noise/unassigned) remainder to 1.0
                    axB.barh(y, 1.0 - left, left=left, color="#E5E5E5",
                             edgecolor="white", linewidth=0.6)
            axB.set_yticks(ypos)
            axB.set_yticklabels(keys, fontsize=7.5)
            axB.set_xlim(0, 1.0)
            axB.set_xlabel("Host read fraction within sample × marker")
            axB.set_title("Resolved mixed-host blood meals")
            axB.legend([plt.Rectangle((0, 0), 1, 1, color="#E5E5E5")],
                       ["unresolved (noise/unassigned)"], loc="lower right", fontsize=6.5)
        else:
            axB.text(0.5, 0.5, "No mixed-host field samples detected", ha="center", va="center")
            axB.axis("off")
    panel_label(axB, "B", dx=-0.02)

    # annotate mixed-host control recovery rate from the RAMBO rollup
    note = ""
    if not rambo.empty:
        r = dict(zip(rambo["metric"], rambo["value"]))
        rate = safe_float(r.get("mixed_host_recovery_rate"))
        exp = r.get("mixed_expected_hosts")
        rec = r.get("mixed_expected_hosts_recovered")
        if rate is not None:
            note = (f"Mixed-host control sensitivity: {rec}/{exp} declared hosts recovered "
                    f"({rate*100:.0f}%). Thresholds not yet benchmarked.")
    if note:
        fig.text(0.5, -0.02, note, ha="center", fontsize=7.4, style="italic", color="#444")

    caption = ("Figure 6. Mixed-host feeding evidence. (A) Composition of sample × marker groups by "
               "feeding type (single-host, mixed-host, or no host signal). (B) The mixed-host blood "
               "meals detected in field samples, each shown as the within-group read fractions of its "
               "co-occurring hosts. Mixed-host detection thresholds are not yet benchmarked and host "
               "fractions are evidence summaries (see Limitations). Inputs: host_call_table.tsv, "
               "rambo_model_summary.tsv.")
    return save_figure(fig, outdir, "figure_06_mixed_host", formats, manifest, caption,
                        ["host_call_table.tsv", "rambo_model_summary.tsv"])


# --------------------------------------------------------------------------------------
# Figure 7 — controls & contamination
# --------------------------------------------------------------------------------------
def fig_controls_contam(control_check, contam_summary, qc_bg, outdir, formats, manifest):
    fig = plt.figure(figsize=(9.6, 4.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 1.0], wspace=0.34)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    # ---- A: positive-control recovery ----
    if not control_check.empty:
        cc = control_check.copy()
        cc["n_exp"] = num(cc["n_expected"]).fillna(0)
        cc["n_rec"] = num(cc["n_recovered"]).fillna(0)
        # barcode is the middle field of the compound sample_uid (run__barcode__sample)
        cc["barcode"] = cc["sample_uid"].map(
            lambda u: u.split("__")[1] if isinstance(u, str) and "__" in u else u)
        cc = cc.sort_values(["control_kind", "barcode"])
        labels = [f"{r.barcode} · {r.control_kind.replace('_control','').replace('_',' ')}"
                  for r in cc.itertuples()]
        ypos = np.arange(len(cc))[::-1]
        axA.barh(ypos, cc["n_exp"], color="#DDDDDD", height=0.62, label="expected hosts")
        axA.barh(ypos, cc["n_rec"], color="#55A868", height=0.40, label="recovered")
        status_color = {"pass": "#2A9D8F", "mixed_pass_all": "#2A9D8F", "pass_genus": "#2A9D8F",
                        "mixed_partial": "#E9C46A", "fail_no_host_signal": "#C44E52",
                        "fail_unexpected_host": "#C44E52", "mixed_fail": "#C44E52"}
        for y, r in zip(ypos, cc.itertuples()):
            c = status_color.get(r.status, "#888888")
            axA.text(max(r.n_exp, r.n_rec) + 0.07, y, r.status.replace("_", " "),
                     va="center", fontsize=6.8, color=c, fontweight="bold")
        axA.set_yticks(ypos)
        axA.set_yticklabels(labels, fontsize=7.2)
        axA.set_xlabel("Number of host taxa")
        axA.set_xlim(0, cc["n_exp"].max() + 1.4 if len(cc) else 1)
        axA.set_title("Positive-control host recovery")
        axA.legend(loc="lower right", fontsize=6.8)
    panel_label(axA, "A", dx=-0.03)

    # ---- B: contamination / decontam summary ----
    axB.axis("off")
    lines = ["Contamination screening", ""]
    if not contam_summary.empty:
        cm = dict(zip(contam_summary["metric"], contam_summary["value"]))
        lines += [
            f"phyloseq:  {cm.get('phyloseq_status','?')}",
            f"decontam:  {cm.get('decontam_status','?')}",
            f"decontam threshold:  {cm.get('decontam_threshold','?')}",
            f"features tested:  {cm.get('n_features','?')}",
            f"flagged contaminants:  {cm.get('n_contaminants','?')}",
        ]
    if not qc_bg.empty:
        lines += ["", "Negative-control background:"]
        for r in qc_bg.itertuples():
            host = getattr(r, "host_assignment", "?")
            mk = getattr(r, "marker", "?")
            cnt = getattr(r, "max_negative_control_count", "?")
            lines.append(f"  • {italic_plain(host)} ({mk}): max {cnt} reads")
    box_txt = "\n".join(lines)
    axB.text(0.0, 1.0, box_txt, va="top", ha="left", fontsize=8.2, family="DejaVu Sans",
             linespacing=1.6,
             bbox=dict(boxstyle="round,pad=0.6", facecolor="#F4F6F8", edgecolor="#CCCCCC"))
    axB.set_title("Laboratory background", loc="left")
    panel_label(axB, "B", dx=0.02)

    caption = ("Figure 7. Controls and contamination. (A) For each positive control, the number of "
               "declared (expected) host taxa versus the number recovered, annotated with the recovery "
               "status; mixed-host controls are lab-prepared mixtures. (B) Laboratory-background "
               "screening: phyloseq/decontam status, the prevalence threshold, the number of features "
               "tested and flagged, and the per-marker negative-control read background. Inputs: "
               "positive_control_check.tsv, contamination_model_summary.tsv, qc_background_thresholds.tsv.")
    return save_figure(fig, outdir, "figure_07_controls_contamination", formats, manifest, caption,
                        ["positive_control_check.tsv", "contamination_model_summary.tsv",
                         "qc_background_thresholds.tsv"])


def italic_plain(name):
    """Plain-text label (used inside a text box where mathtext is undesirable)."""
    return name


# --------------------------------------------------------------------------------------
# Figure 8 — ecological stratification
# --------------------------------------------------------------------------------------
def fig_ecology(master, outdir, formats, manifest):
    if master.empty:
        print("  [skip] figure 08: empty master endpoint")
        return []
    field = master[(master["sample_type"] == "sample") &
                   (~master["host_assignment"].isin(["", "unassigned"]))].copy()
    if field.empty:
        print("  [skip] figure 08: no field host assignments")
        return []

    def host_by(col, ax, title):
        levels = [l for l in field[col].dropna().unique() if l != ""]
        # order levels by sample count
        levels = sorted(levels, key=lambda l: -field[field[col] == l]["sample_uid"].nunique())
        hosts = [h for h in HOST_ORDER if h in set(field["host_assignment"])]
        hosts += sorted(set(field["host_assignment"]) - set(hosts))
        x = np.arange(len(levels))
        bottoms = np.zeros(len(levels))
        for h in hosts:
            vals = []
            for l in levels:
                vals.append(field[(field[col] == l) & (field["host_assignment"] == h)]["sample_uid"].nunique())
            vals = np.array(vals)
            if vals.sum() > 0:
                ax.bar(x, vals, bottom=bottoms, color=host_color(h), width=0.7,
                       edgecolor="white", linewidth=0.5)
                bottoms += vals
        n_by_level = [field[field[col] == l]["sample_uid"].nunique() for l in levels]
        ax.set_xticks(x)
        ax.set_xticklabels([f"{l.replace('_',' ')}\n(n={n})" for l, n in zip(levels, n_by_level)],
                           fontsize=7.4)
        ax.set_ylabel("Host detections (samples)")
        ax.set_title(title)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.6, 4.6))
    host_by("collection_region", axA, "Hosts by ecological zone")
    panel_label(axA, "A")
    host_by("sibling_species", axB, "Hosts by vector sibling species")
    # italicise vector species tick labels
    for lbl in axB.get_xticklabels():
        lbl.set_fontstyle("italic")
    panel_label(axB, "B")

    hosts = [h for h in HOST_ORDER if h in set(field["host_assignment"])]
    hosts += sorted(set(field["host_assignment"]) - set(hosts))
    handles = [plt.Rectangle((0, 0), 1, 1, color=host_color(h)) for h in hosts]
    fig.legend(handles, [italicise_taxon(h) for h in hosts], loc="lower center",
               ncol=min(len(hosts), 6), fontsize=7.5, bbox_to_anchor=(0.5, -0.04), title="Host")
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    caption = ("Figure 8. Host detections across ecologies and vectors. Stacked bars give the number of "
               "field samples in which each host was detected, broken down by (A) ecological zone and "
               "(B) Anopheles gambiae s.l. sibling species (per-group sample sizes shown). These "
               "describe the demonstration cohort and are not formal association tests. Input: "
               "bloodmeal_master_endpoint.tsv.")
    return save_figure(fig, outdir, "figure_08_ecology", formats, manifest, caption,
                        ["bloodmeal_master_endpoint.tsv"])


# --------------------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Build HÆMA publication figures from endpoint tables.")
    ap.add_argument("--endpoint-dir", required=True, help="05_endpoint_files directory")
    ap.add_argument("--reports-dir", default="", help="06_reports directory")
    ap.add_argument("--manifest", default="", help="run_manifest.json (for caption provenance)")
    ap.add_argument("--outdir", required=True, help="output directory for figures")
    ap.add_argument("--formats", default="pdf,svg,png", help="comma-separated: pdf,svg,png")
    args = ap.parse_args()

    ep = Path(args.endpoint_dir)
    rp = Path(args.reports_dir) if args.reports_dir else ep.parent / "06_reports"
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    formats = [f.strip() for f in args.formats.split(",") if f.strip()]

    # Load all available tables (missing ones degrade gracefully).
    master = read_tsv(ep / "bloodmeal_master_endpoint.tsv")
    host_calls = read_tsv(ep / "host_call_table.tsv")
    sample_sum = read_tsv(ep / "sample_level_summary.tsv")
    qc = read_tsv(ep / "qc_summary.tsv")
    control_check = read_tsv(rp / "positive_control_check.tsv")
    rambo = read_tsv(rp / "rambo_model_summary.tsv")
    contam_summary = read_tsv(rp / "contamination_model_summary.tsv")
    qc_bg = read_tsv(rp / "qc_background_thresholds.tsv")

    print(f"Loaded: master={len(master)} host_calls={len(host_calls)} samples={len(sample_sum)} "
          f"qc={len(qc)} controls={len(control_check)}")

    manifest = []
    figures = [
        ("figure_01_workflow", lambda: fig_workflow(outdir, formats, manifest)),
        ("figure_02_sequencing_qc", lambda: fig_qc(qc, master, outdir, formats, manifest)),
        ("figure_03_depth_denoising", lambda: fig_depth_denoise(sample_sum, qc, outdir, formats, manifest)),
        ("figure_04_host_assignment", lambda: fig_host_overview(host_calls, master, outdir, formats, manifest)),
        ("figure_05_host_composition", lambda: fig_host_composition(host_calls, outdir, formats, manifest)),
        ("figure_06_mixed_host", lambda: fig_mixed_host(host_calls, rambo, outdir, formats, manifest)),
        ("figure_07_controls_contamination",
         lambda: fig_controls_contam(control_check, contam_summary, qc_bg, outdir, formats, manifest)),
        ("figure_08_ecology", lambda: fig_ecology(master, outdir, formats, manifest)),
    ]
    failures = 0
    for name, fn in figures:
        try:
            fn()
        except Exception:
            failures += 1
            print(f"  [FAIL] {name}:\n{traceback.format_exc()}", file=sys.stderr)

    # captions + manifest
    if manifest:
        man_df = pd.DataFrame(manifest)
        man_df.to_csv(outdir / "figure_manifest.tsv", sep="\t", index=False)
        with (outdir / "figure_captions.md").open("w") as fh:
            fh.write("# HÆMA figure captions (draft)\n\n")
            run = ""
            if args.manifest and Path(args.manifest).exists():
                try:
                    run = json.load(open(args.manifest)).get("pipeline_version", "")
                except Exception:
                    run = ""
            if run:
                fh.write(f"Generated from HÆMA endpoint outputs (pipeline {run}).\n\n")
            for row in manifest:
                fh.write(f"**{row['figure']}** — files: `{row['files']}`\n\n")
                fh.write(f"{row['caption']}\n\n")
                fh.write(f"_Input files:_ {row['inputs']}\n\n---\n\n")
        print(f"\nWrote {len(manifest)} figures + figure_captions.md + figure_manifest.tsv to {outdir}")
    if failures:
        # A figure-rendering step should not sink the whole pipeline: warn loudly but only fail
        # hard when *no* figure could be produced (a genuinely broken step / missing inputs).
        print(f"{failures} figure(s) failed to render (see traceback above)", file=sys.stderr)
    if not manifest:
        print("No figures were produced — check that endpoint tables exist.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
