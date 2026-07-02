#!/usr/bin/env python3
"""Render the 3 candidate replacements for Figure 3 Panel B side-by-side for a decision."""
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import figure_style as S


def zget(zt, z, m, col="value"):
    r = zt[(zt.stratum == z) & (zt.metric == m)]
    return float(r[col].iloc[0]) if len(r) else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eco-bioclim", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    S.apply_house_style()
    eco = pd.read_csv(a.eco_bioclim, sep="\t")
    zt = eco[eco.stratum_type == "ecological_zone"].copy()
    zones = [z for z in S.ZONE_ORDER if z in set(zt.stratum)]
    zlab = [S.zone_label(z) for z in zones]

    fig = plt.figure(figsize=(8.4, 3.0))
    gs = fig.add_gridspec(1, 3, wspace=0.5)
    axA, axB, axC = (fig.add_subplot(gs[0, i]) for i in range(3))

    # ---- Option A: grouped bars + Wilson CIs (HBI / BBI / zoophily) ----
    metrics = [("human_blood_index", "HBI (human)", "#0072B2"),
               ("host_blood_index::Bos taurus", "BBI (cattle)", "#D55E00"),
               ("animal_blood_index_zoophily", "Zoophily", "#117733")]
    x = np.arange(len(zones)); bw = 0.26
    for i, (m, lbl, c) in enumerate(metrics):
        vals = [zget(zt, z, m) for z in zones]
        lo = [vals[j] - zget(zt, z, m, "ci_low") for j, z in enumerate(zones)]
        hi = [zget(zt, z, m, "ci_high") - vals[j] for j, z in enumerate(zones)]
        axA.bar(x + (i - 1) * bw, vals, width=bw, color=c, label=lbl)
        axA.errorbar(x + (i - 1) * bw, vals, yerr=[lo, hi], fmt="none",
                     ecolor="black", elinewidth=0.5, capsize=1.5)
    axA.set_xticks(x); axA.set_xticklabels(zlab, fontsize=5.5, rotation=45, ha="right")
    axA.set_ylim(0, 1.0); axA.set_ylabel("Index value")
    axA.legend(fontsize=5, loc="upper center", ncol=1)
    axA.set_title("Option A — grouped bars + CI", fontsize=7)

    # ---- Option B: 100% stacked feeding-type composition ----
    parts = [("human_only_fraction", "human only", "#0072B2"),
             ("mixed_human_animal_fraction", "mixed human+animal", "#9467bd"),
             ("animal_only_fraction", "animal only", "#D55E00")]
    for j, z in enumerate(zones):
        bottom = 0.0
        for m, lbl, c in parts:
            v = zget(zt, z, m)
            axB.bar(j, v, bottom=bottom, color=c, width=0.7,
                    edgecolor="white", linewidth=0.5)
            if v > 0.08:
                axB.text(j, bottom + v / 2, f"{v*100:.0f}", ha="center", va="center",
                         color="white", fontsize=5)
            bottom += v
    axB.set_xticks(range(len(zones))); axB.set_xticklabels(zlab, fontsize=5.5, rotation=45, ha="right")
    axB.set_ylim(0, 1.0); axB.set_ylabel("Proportion of meals")
    axB.legend(handles=[Patch(facecolor=c, label=l) for _, l, c in parts],
               fontsize=5, loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=1)
    axB.set_title("Option B — stacked feeding-type", fontsize=7)

    # ---- Option C: dot / forest plot (HBI + zoophily, CIs), zones on y ----
    yy = np.arange(len(zones))[::-1]
    for j, z in zip(yy, zones):
        for m, c, mk, off in (("human_blood_index", "#0072B2", "o", 0.12),
                              ("animal_blood_index_zoophily", "#117733", "o", -0.12)):
            v = zget(zt, z, m); lo = zget(zt, z, m, "ci_low"); hi = zget(zt, z, m, "ci_high")
            face = c if m == "human_blood_index" else "white"
            axC.plot([lo, hi], [j + off, j + off], color=c, lw=1.0, zorder=1)
            axC.plot(v, j + off, mk, ms=5, color=c, mfc=face, mec=c, zorder=2)
    axC.set_yticks(yy[::-1]); axC.set_yticklabels(zlab[::-1], fontsize=5.5)
    axC.set_xlim(0, 1.0); axC.set_xlabel("Index value (95% CI)")
    axC.legend(handles=[Line2D([0], [0], marker="o", color="#0072B2", lw=0, label="anthropophily (HBI)"),
                        Line2D([0], [0], marker="o", color="#117733", mfc="white", lw=0, label="zoophily")],
               fontsize=5, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    axC.set_title("Option C — dot/forest plot", fontsize=7)

    fig.suptitle("Figure 3 Panel B — candidate replacements (same data, +95% Wilson CI)",
                 fontweight="bold", fontsize=8, y=1.04)
    for ext in ("png",):
        import os
        os.makedirs(a.outdir, exist_ok=True)
        fig.savefig(os.path.join(a.outdir, f"panelB_options.{ext}"), dpi=170, bbox_inches="tight")
    print("wrote panelB_options.png")


if __name__ == "__main__":
    main()
