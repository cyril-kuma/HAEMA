#!/usr/bin/env python3
"""
HÆMA supplementary figures (Objective 1), house style. Renders plot-ready tables only.
S1 host accumulation (rarefaction); S2 host blood index x bioclimatic zone matrix;
B1 (Appendix B) Pianka niche overlap + Bray-Curtis zonal host-use turnover (moved out of Fig 5).
(Methods/QC supplementary figures — workflow, sequencing QC, controls — are the pipeline's
auto-generated figures, kept as-is.)
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import figure_style as S


def s1_rarefaction(fd, outdir):
    d = pd.read_csv(os.path.join(fd, "rarefaction.tsv"), sep="\t")
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    ax.fill_between(d.n_samples, d.mean_taxa - d.sd_taxa, d.mean_taxa + d.sd_taxa,
                    color="#4477AA", alpha=0.25)
    ax.plot(d.n_samples, d.mean_taxa, color="#4477AA", lw=1.2)
    ax.set_xlabel("Host-identified samples")
    ax.set_ylabel("Vertebrate host taxa")
    ax.set_title("Host accumulation", fontsize=7)
    ax.set_xlim(0, d.n_samples.max())
    ax.set_ylim(0, d.mean_taxa.max() * 1.08)
    S.panel_label(ax, "")
    S.save(fig, outdir, "figure_S1_rarefaction")


def s2_host_zone_matrix(eco_bioclim, outdir):
    zt = eco_bioclim[eco_bioclim.stratum_type == "ecological_zone"].copy()
    hb = zt[zt.metric.str.startswith("host_blood_index::")].copy()
    hb["host"] = hb.metric.str.split("::").str[1]
    zones = [z for z in S.ZONE_ORDER if z in set(hb.stratum)]
    hosts = [h for h in S.HOST_ORDER if h in set(hb.host)]
    mat = (hb.pivot_table(index="host", columns="stratum", values="value", aggfunc="first")
           .reindex(index=hosts, columns=zones).fillna(0.0))
    fig, ax = plt.subplots(figsize=(3.6, 3.2))
    im = ax.imshow(mat.values, cmap=S.SEQ_CMAP, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(zones)))
    ax.set_xticklabels([S.zone_label(z) for z in zones], fontsize=5.5, rotation=45, ha="right")
    ax.set_yticks(range(len(hosts)))
    ax.set_yticklabels([S.two_line(h) for h in hosts], fontsize=5)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            if v > 0.005:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5,
                        color="white" if v > 0.6 else "black")
    cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("Host blood index", fontsize=6)
    ax.set_title("Host blood index by bioclimatic zone", fontsize=7)
    S.save(fig, outdir, "figure_S2_host_zone_matrix")


def _heat(ax, mat, labels_x, labels_y, cmap, vmin, vmax, title, fmt="{:.2f}"):
    im = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")
    ax.set_xticks(range(len(labels_x)))
    ax.set_xticklabels(labels_x, fontsize=5.5, rotation=45, ha="right")
    ax.set_yticks(range(len(labels_y)))
    ax.set_yticklabels(labels_y, fontsize=5.5)
    span = vmax - vmin
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=5,
                    color="white" if (v - vmin) / span > 0.55 else "black")
    ax.set_title(title, fontsize=7)
    return im


def figureB1_overlap_turnover(fd, outdir):
    """Appendix B: (A) Pianka niche-overlap between sibling species on the FULL 0-1 scale
    (values cluster near 1 -> near-complete dietary overlap), (B) Bray-Curtis host-use
    turnover between bioclimatic zones. Single sequential colourmap; moved out of main Fig 5."""
    ov = pd.read_csv(os.path.join(fd, "niche_overlap_matrix.tsv"), sep="\t").set_index("species")
    sp = [s for s in S.SPECIES_HEX if s in ov.index]
    ovm = ov.reindex(index=sp, columns=sp).values.astype(float)

    bd = pd.read_csv(os.path.join(fd, "beta_diversity_matrix.tsv"), sep="\t")
    zones = [z for z in S.ZONE_ORDER if z in set(bd.zone_a) | set(bd.zone_b)]
    bc = pd.DataFrame(0.0, index=zones, columns=zones)
    for _, r in bd.iterrows():
        if r.zone_a in zones and r.zone_b in zones:
            bc.loc[r.zone_a, r.zone_b] = r.bray_curtis
            bc.loc[r.zone_b, r.zone_a] = r.bray_curtis

    fig = plt.figure(figsize=(7.0, 3.1))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.55)
    axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])

    splbl = [S.italic(S.SPECIES_LABEL[s]) for s in sp]
    imA = _heat(axA, ovm, splbl, splbl, S.SEQ_CMAP, 0.0, 1.0,
                "Pianka niche overlap (0–1)")
    cbA = fig.colorbar(imA, ax=axA, shrink=0.78, pad=0.03)
    cbA.set_label("Pianka overlap", fontsize=6)
    axA.text(0.5, -0.40, "full 0–1 scale: all pairs ≈1, i.e. near-complete host-use overlap",
             transform=axA.transAxes, ha="center", va="top", fontsize=4.8, color="0.4",
             style="italic")
    S.panel_label(axA, "A", dx=-0.16)

    zlbl = [S.zone_label(z) for z in zones]
    imB = _heat(axB, bc.values, zlbl, zlbl, S.SEQ_CMAP, 0.0, 1.0,
                "Bray–Curtis host-use turnover")
    cbB = fig.colorbar(imB, ax=axB, shrink=0.78, pad=0.03)
    cbB.set_label("Bray–Curtis dissimilarity", fontsize=6)
    S.panel_label(axB, "B", dx=-0.18)

    fig.suptitle("Appendix B — sibling-species niche overlap and zonal host-use turnover",
                 fontweight="bold", fontsize=8, y=1.04)
    S.save(fig, outdir, "figure_B1_niche_overlap_turnover")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figure-data", required=True)
    ap.add_argument("--eco-bioclim", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    S.apply_house_style()
    print(f"supplementary figures -> {a.outdir}")
    s1_rarefaction(a.figure_data, a.outdir)
    s2_host_zone_matrix(pd.read_csv(a.eco_bioclim, sep="\t"), a.outdir)
    figureB1_overlap_turnover(a.figure_data, a.outdir)
    print("done.")


if __name__ == "__main__":
    main()
