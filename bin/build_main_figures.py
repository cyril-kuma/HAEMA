#!/usr/bin/env python3
"""
HÆMA main publication figures (Objective 1) — 5 self-contained figures for thesis/manuscript/
conference. Renders ONLY the plot-ready tables from figure_data_prep.py + the ecological indices
(overall + bioclimatic) + GADM Ghana GeoJSON. House style from figure_style.py.

Fig 1  Study system: (A) Ghana bioclimatic-zone map + sites, (B) sibling-species composition by site.
Fig 2  Metabarcoding efficacy: (A) read length by marker, (B) marker recovery, (C) assignment confidence.
Fig 3  Host use across zones: (A) overall host blood indices, (B) HBI/BBI/zoophily slopegraph by zone.
Fig 4  Mixed feeding: (A) UpSet-style host combinations, (B) mixed-feeding rate by zone.
Fig 5  Host-use ecology: (A) zooprophylaxis by zone, (B) Levins' niche breadth by sibling species,
       (C) vector-host bipartite feeding network. (Pianka overlap + Bray-Curtis turnover -> Appendix B.)
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import figure_style as S


def load(fd, name):
    return pd.read_csv(os.path.join(fd, name), sep="\t")


def eco_rows(df, stype):
    return df[df["stratum_type"] == stype].copy()


# ============================================================ Figure 1 =========
def figure1(fd, geo_dir, outdir):
    import geopandas as gpd
    from shapely.geometry import Point
    sites = load(fd, "sites.tsv")
    sp = load(fd, "sibling_species_by_site.tsv")

    fig = plt.figure(figsize=(7.0, 3.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.05)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    # --- A: map ---
    g1 = gpd.read_file(os.path.join(geo_dir, "gadm41_GHA_1.json"))
    g1["zone"] = g1["NAME_1"].map(S.REGION_ZONE)
    zones = g1.dissolve(by="zone").reset_index().to_crs(32630)
    country = gpd.read_file(os.path.join(geo_dir, "gadm41_GHA_0.json")).to_crs(32630)
    for _, row in zones.iterrows():
        gpd.GeoSeries([row.geometry]).plot(ax=axA, color=S.zone_color(row["zone"]),
                                           alpha=0.55, edgecolor="white", linewidth=0.5)
    country.boundary.plot(ax=axA, color="0.3", linewidth=0.7)
    pts = gpd.GeoDataFrame(sites, geometry=[Point(xy) for xy in zip(sites.longitude, sites.latitude)],
                           crs=4326).to_crs(32630)
    for _, r in pts.iterrows():
        axA.plot(r.geometry.x, r.geometry.y, "o", ms=4, mec="black", mew=0.5,
                 color=S.zone_color(r["zone"]), zorder=5)
        axA.annotate(r["site"], (r.geometry.x, r.geometry.y), xytext=(3, 2),
                     textcoords="offset points", fontsize=5, zorder=6,
                     rotation=45, rotation_mode="anchor", ha="left", va="bottom")
    # scale bar (100 km) + north arrow (projected CRS -> metres)
    x0, y0 = axA.get_xlim()[0], axA.get_ylim()[0]
    xr = axA.get_xlim()[1] - x0
    bx = x0 + 0.06 * xr
    by = y0 + 0.05 * (axA.get_ylim()[1] - y0)
    axA.plot([bx, bx + 100000], [by, by], color="black", lw=1.5)
    axA.text(bx + 50000, by + 8000, "100 km", ha="center", fontsize=5)
    axA.annotate("N", xy=(x0 + 0.06 * xr, axA.get_ylim()[1] - 0.10 * (axA.get_ylim()[1] - y0)),
                 xytext=(x0 + 0.06 * xr, axA.get_ylim()[1] - 0.22 * (axA.get_ylim()[1] - y0)),
                 ha="center", fontsize=7, fontweight="bold",
                 arrowprops=dict(arrowstyle="-|>", color="black", lw=1.2))
    axA.legend(handles=[Patch(facecolor=S.zone_color(z), alpha=0.55, edgecolor="white",
               label=S.zone_label(z)) for z in S.ZONE_ORDER], loc="lower right",
               fontsize=5, title="Bioclimatic zone", title_fontsize=5.5)
    axA.set_axis_off()
    axA.set_aspect("equal")
    S.panel_label(axA, "A", dx=0.0, dy=1.02)

    # --- B: sibling-species composition by site (100% stacked, horizontal, grouped by zone) ---
    species = [s for s in S.SPECIES_HEX if s in set(sp["sibling_species"])]
    sites_order, ylabels, yrows = [], [], []
    y = 0
    yticks, yticklabels = [], []
    for z in S.ZONE_ORDER:
        zs = sorted(sp[sp.zone == z]["site"].unique(),
                    key=lambda s: -sp[(sp.zone == z) & (sp.site == s)]["n"].sum())
        for s in zs:
            sub = sp[(sp.zone == z) & (sp.site == s)]
            tot = sub["n"].sum()
            left = 0.0
            for spc in species:
                v = sub[sub.sibling_species == spc]["n"].sum() / tot if tot else 0
                if v > 0:
                    axB.barh(y, v, left=left, color=S.SPECIES_HEX[spc], edgecolor="white",
                             linewidth=0.4, height=0.78)
                    # species identity is shown by colour + the italic legend below;
                    # no in-bar text (avoids clipping "coluzzii" -> "coluzzi").
                    left += v
            axB.text(1.01, y, f"{int(tot)}", va="center", fontsize=5, color="0.3")
            yticks.append(y)
            yticklabels.append(s)
            y += 1
        y += 0.6  # gap between zones
    axB.set_yticks(yticks)
    axB.set_yticklabels(yticklabels)
    axB.set_xlim(0, 1)
    axB.set_xlabel("Relative abundance")
    axB.set_title("Sibling-species composition", fontsize=7)
    axB.invert_yaxis()
    axB.legend(handles=[Patch(facecolor=S.SPECIES_HEX[s], label=S.italic(S.SPECIES_LABEL[s]))
               for s in species], loc="upper center", bbox_to_anchor=(0.5, -0.12),
               ncol=len(species), fontsize=5, columnspacing=1.0, handletextpad=0.4)
    S.panel_label(axB, "B", dx=-0.18, dy=1.02)
    fig.suptitle("Study system: Ghanaian bioclimatic zones, sampling sites and vector composition",
                 fontweight="bold", fontsize=8, y=1.02)
    S.save(fig, outdir, "figure_1_study_system")


# ============================================================ Figure 2 =========
def figure2(fd, outdir):
    rl = load(fd, "read_length_distribution.tsv")
    amp = load(fd, "amplicon_sizes.tsv")
    mr = load(fd, "marker_recovery.tsv")
    ac = load(fd, "assignment_confidence.tsv")

    fig = plt.figure(figsize=(7.2, 2.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.0, 1.0], wspace=0.42)
    axA, axB, axC = (fig.add_subplot(gs[0, i]) for i in range(3))

    # A: read-length by marker
    series = ["cyt_b", "co1_short", "co1_long", "unassigned"]
    scol = {"cyt_b": "#332288", "co1_short": "#117733", "co1_long": "#CC6677", "unassigned": "#BBBBBB"}
    for s in series:
        d = rl[rl.series == s].sort_values("length_bp")
        if d.empty:
            continue
        axA.fill_between(d.length_bp, d["count"], step="mid", alpha=0.35, color=scol[s])
        axA.step(d.length_bp, d["count"], where="mid", color=scol[s], lw=0.9,
                 label=S.MARKER_LABEL.get(s, "unassigned"))
    for _, r in amp.iterrows():
        axA.axvline(r["expected_bp"], ls="--", color=scol[r["marker"]], lw=0.6, alpha=0.8)
    axA.set_xlim(0, 900)
    axA.set_xlabel("Read length (bp)")
    axA.set_ylabel("Reads")
    axA.set_title("Read-length by marker", fontsize=7)
    axA.legend(fontsize=5, loc="upper right")
    S.panel_label(axA, "A")

    # B: marker recovery (grouped horizontal bars by marker)
    markers = ["cyt_b", "co1_short", "co1_long"]
    mcol = {"cyt_b": "#332288", "co1_short": "#117733", "co1_long": "#CC6677"}
    hosts = [h for h in S.HOST_ORDER if h in set(mr["host"])]
    yb = np.arange(len(hosts))
    bw = 0.26
    for i, mk in enumerate(markers):
        vals = [mr[(mr.host == h) & (mr.marker == mk)]["n_samples"].sum() for h in hosts]
        axB.barh(yb + (1 - i) * bw, vals, height=bw, color=mcol[mk], label=S.MARKER_LABEL[mk])
    axB.set_yticks(yb)
    axB.set_yticklabels([S.common(h) for h in hosts])
    axB.invert_yaxis()
    axB.set_xlabel("Samples recovered")
    axB.set_title("Marker complementarity", fontsize=7)
    axB.legend(fontsize=5, loc="lower right")
    S.panel_label(axB, "B")

    # C: assignment confidence — 2D density (all calls already pass >=97% id / >=80% cov,
    # so a scatter just piles into the top-right corner; hexbin shows where the mass sits).
    hb = axC.hexbin(ac.pident, ac.coverage, gridsize=24, cmap=S.SEQ_CMAP, mincnt=1,
                    linewidths=0.0)
    axC.set_xlim(97, 100.2)
    axC.set_ylim(78, 103)
    axC.axvline(97, ls="--", color="#C44E52", lw=0.6)
    axC.axhline(80, ls="--", color="#C44E52", lw=0.6)
    axC.set_xlabel("BLAST identity (%)")
    axC.set_ylabel("Query coverage (%)")
    axC.set_title("Assignment confidence", fontsize=7)
    cb = fig.colorbar(hb, ax=axC, pad=0.02, fraction=0.055, aspect=18)
    cb.set_label("Assignments", fontsize=5.5)
    cb.ax.tick_params(labelsize=5)
    axC.text(0.97, 0.04, f"n = {len(ac):,} assignments\nall ≥97% id, ≥80% cov",
             transform=axC.transAxes, ha="right", va="bottom", fontsize=4.6, color="0.35")
    S.panel_label(axC, "C")
    fig.suptitle("Multi-marker metabarcoding efficacy and assignment validation",
                 fontweight="bold", fontsize=8, y=1.04)
    S.save(fig, outdir, "figure_2_metabarcoding_efficacy")


# ============================================================ Figure 3 =========
def figure3(eco_overall, eco_bioclim, outdir):
    fig = plt.figure(figsize=(7.0, 3.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.45)
    axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])

    # A: overall host blood indices (descending, Wilson CI)
    ov = eco_rows(eco_overall, "overall")
    hbi = ov[ov.metric.str.startswith("host_blood_index::")].copy()
    hbi["host"] = hbi.metric.str.split("::").str[1]
    hbi = hbi.sort_values("value")
    y = np.arange(len(hbi))
    axA.barh(y, hbi.value, color=[S.host_color(h) for h in hbi.host], height=0.72,
             edgecolor="black", linewidth=0.3)
    axA.errorbar(hbi.value, y, xerr=[hbi.value - hbi.ci_low, hbi.ci_high - hbi.value],
                 fmt="none", ecolor="black", elinewidth=0.6, capsize=1.5)
    axA.set_yticks(y)
    axA.set_yticklabels([S.two_line(h) for h in hbi.host], fontsize=5.5)
    axA.set_xlim(0, 1)
    axA.set_xlabel("Host blood index")
    axA.set_title("Overall host use", fontsize=7)
    S.panel_label(axA, "A", dx=-0.42)

    # B: HBI / BBI / zoophily slopegraph across zones (south -> north)
    zt = eco_rows(eco_bioclim, "ecological_zone")
    zones = [z for z in S.ZONE_ORDER if z in set(zt.stratum)]
    series = [("human_blood_index", "HBI (human)", "#0072B2"),
              ("host_blood_index::Bos taurus", "BBI (cattle)", "#D55E00"),
              ("animal_blood_index_zoophily", "Zoophily", "#117733")]
    xs = np.arange(len(zones))
    for metric, lbl, col in series:
        ys = [zt[(zt.stratum == z) & (zt.metric == metric)]["value"].mean() for z in zones]
        axB.plot(xs, ys, "-o", color=col, ms=4, lw=1.2)
        axB.text(xs[-1] + 0.06, ys[-1], lbl, color=col, va="center", fontsize=5.5)
    axB.set_xticks(xs)
    axB.set_xticklabels([S.zone_label(z) + ("†" if z in S.ZONE_DESCRIPTIVE else "")
                         for z in zones], fontsize=5.5, rotation=45, ha="right")
    axB.set_xlim(-0.2, len(zones) - 1 + 1.05)
    axB.set_ylim(0, 1.0)
    axB.set_ylabel("Index value")
    axB.set_title("Anthropophily ↔ zoophily by zone", fontsize=7)
    S.panel_label(axB, "B", dx=-0.2)
    if any(z in S.ZONE_DESCRIPTIVE for z in zones):
        axB.text(1.0, -0.34, "† descriptive only (small n)", transform=axB.transAxes,
                 ha="right", va="top", fontsize=4.8, color="0.4", style="italic")
    fig.suptitle("Vertebrate host use across Ghanaian bioclimatic zones",
                 fontweight="bold", fontsize=8, y=1.02)
    S.save(fig, outdir, "figure_3_host_use_by_zone")


# ============================================================ Figure 4 =========
def figure4(fd, eco_bioclim, outdir):
    mc = load(fd, "mixed_combinations.tsv").head(9)
    mult = load(fd, "feeding_multiplicity.tsv")
    fig = plt.figure(figsize=(7.4, 3.1))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.6, 1.0], wspace=0.4)
    axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])

    # A: UpSet-style — membership dot-matrix (left) + count bars (right)
    hosts = []
    for combo in mc.combination:
        for h in combo.split(" + "):
            if h not in hosts:
                hosts.append(h)
    hosts = [h for h in ["Human", "Cattle", "Sheep", "Goat", "Pig", "Dog", "Donkey",
                         "Chicken", "Cat"] if h in hosts]
    n = len(mc)
    yk = np.arange(n)[::-1]                    # top combo at top
    xcol = np.arange(len(hosts))
    matrix_w = len(hosts)
    for yi, (_, row) in zip(yk, mc.iterrows()):
        members = row.combination.split(" + ")
        present_x = [xcol[i] for i, h in enumerate(hosts) if h in members]
        for xi, h in enumerate(hosts):
            axA.plot(xcol[xi], yi, "o", ms=5,
                     color=("#333333" if h in members else "#dddddd"), zorder=2)
        if present_x:
            axA.plot([min(present_x), max(present_x)], [yi, yi], color="#333333", lw=1.2, zorder=1)
    # count bars to the right of the matrix
    bar_x0 = matrix_w + 0.4
    maxc = mc["count"].max()
    barscale = (len(hosts) * 0.9) / maxc
    for yi, (_, row) in zip(yk, mc.iterrows()):
        col = "#1b7837" if row["type"] == "animal_animal" else "#5e3c99"
        axA.barh(yi, row["count"] * barscale, left=bar_x0, height=0.55, color=col)
        axA.text(bar_x0 + row["count"] * barscale + 0.1, yi, int(row["count"]), va="center",
                 fontsize=5)
    axA.set_xticks(list(xcol) + [bar_x0 + maxc * barscale / 2])
    axA.set_xticklabels(hosts + ["mixed meals"], fontsize=5, rotation=45, ha="right")
    axA.tick_params(axis="x", length=0)
    axA.set_yticks([])
    axA.set_xlim(-0.6, bar_x0 + maxc * barscale + 1.2)
    axA.set_ylim(-0.8, n - 0.2)
    axA.set_title("Mixed-host blood-meal combinations", fontsize=7)
    axA.legend(handles=[Patch(facecolor="#5e3c99", label="human + animal"),
                        Patch(facecolor="#1b7837", label="animal + animal")],
               loc="lower right", fontsize=5)
    for sp in ("left", "bottom"):
        axA.spines[sp].set_visible(False)
    S.panel_label(axA, "A", dx=-0.05, dy=1.10)

    # B: mixed-feeding rate by zone
    zt = eco_rows(eco_bioclim, "ecological_zone")
    zones = [z for z in S.ZONE_ORDER if z in set(zt.stratum)]
    mf = zt[zt.metric == "mixed_feeding_rate"].set_index("stratum")
    x = np.arange(len(zones))
    vals = [mf.loc[z, "value"] for z in zones]
    lo = [vals[i] - mf.loc[z, "ci_low"] for i, z in enumerate(zones)]
    hi = [mf.loc[z, "ci_high"] - vals[i] for i, z in enumerate(zones)]
    bars = axB.bar(x, vals, color=[S.zone_color(z) for z in zones], width=0.66,
                   edgecolor="black", linewidth=0.3)
    # hatch + de-emphasise descriptive (small-n) zones so their wide CI is not over-read
    for bi, z in enumerate(zones):
        if z in S.ZONE_DESCRIPTIVE:
            bars[bi].set_hatch("////")
            bars[bi].set_alpha(0.55)
    axB.errorbar(x, vals, yerr=[lo, hi], fmt="none", ecolor="black", elinewidth=0.6, capsize=2)
    axB.set_xticks(x)
    axB.set_xticklabels([S.zone_label(z) + ("†" if z in S.ZONE_DESCRIPTIVE else "")
                         for z in zones], fontsize=5.5, rotation=45, ha="right")
    axB.set_ylim(0, max(0.3, max(vals) * 1.3))
    axB.set_ylabel("Mixed-feeding rate")
    axB.set_title("Mixed feeding by zone", fontsize=7)
    S.panel_label(axB, "B", dx=-0.22)
    if any(z in S.ZONE_DESCRIPTIVE for z in zones):
        axB.text(1.0, -0.42, "† descriptive only (small n)", transform=axB.transAxes,
                 ha="right", va="top", fontsize=4.8, color="0.4", style="italic")
    mixn = int(mult[mult.n_hosts >= 2]["n_samples"].sum())
    totn = int(mult["n_samples"].sum())
    fig.suptitle(f"Mixed-host blood feeding resolved by metabarcoding "
                 f"({mixn}/{totn} = {mixn/totn*100:.0f}% of meals mixed)",
                 fontweight="bold", fontsize=8, y=1.03)
    S.save(fig, outdir, "figure_4_mixed_feeding")


# ============================================================ Figure 5 =========
def figure5(fd, outdir):
    """Host-use ecology (main): (A) zooprophylaxis by zone, (B) Levins' niche breadth by
    sibling species, (C) vector-host bipartite feeding network. Pianka overlap and
    Bray-Curtis turnover are moved to Appendix B (build_supplementary_figures.figureB1)."""
    idx = load(fd, "host_ecology_indices.tsv")
    vhm = load(fd, "vector_host_matrix.tsv").set_index("species")

    fig = plt.figure(figsize=(7.4, 2.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.5], wspace=0.5)
    axA, axB, axC = (fig.add_subplot(gs[0, i]) for i in range(3))

    # A: Zooprophylaxis Index (realized blood-meal diversion) by zone, bootstrap 95% CI.
    zpi = idx[idx["index"] == "zooprophylaxis_index"].dropna(subset=["value"])
    zorder = ["Overall"] + [z for z in S.ZONE_ORDER if z in set(zpi.stratum)]
    zpi = zpi.set_index("stratum").reindex([z for z in zorder if z in set(zpi.stratum)])
    xa = np.arange(len(zpi))
    cols = ["#666666" if s == "Overall" else S.zone_color(s) for s in zpi.index]
    axA.bar(xa, zpi.value, color=cols, width=0.66, edgecolor="black", linewidth=0.3)
    axA.errorbar(xa, zpi.value, yerr=[zpi.value - zpi.ci_low, zpi.ci_high - zpi.value],
                 fmt="none", ecolor="black", elinewidth=0.6, capsize=2)
    axA.set_xticks(xa)
    axA.set_xticklabels([("Overall" if s == "Overall" else S.zone_label(s)) for s in zpi.index],
                        fontsize=6, rotation=45, ha="right")
    axA.set_ylim(0, 1)
    axA.set_ylabel("Zooprophylaxis index")
    axA.set_title("Blood-meal diversion", fontsize=7.5)
    S.panel_label(axA, "A", dx=-0.30)

    # B: Levins' niche breadth (B_A; 0 specialist -> 1 generalist) per sibling species.
    nb = idx[idx["index"] == "niche_breadth_BA"].set_index("stratum")
    sp = [s for s in S.SPECIES_HEX if s in nb.index]
    xb = np.arange(len(sp))
    vals = nb.loc[sp, "value"].values
    axB.bar(xb, vals, color=[S.SPECIES_HEX[s] for s in sp], width=0.62,
            edgecolor="black", linewidth=0.3)
    axB.errorbar(xb, vals, yerr=[vals - nb.loc[sp, "ci_low"].values,
                                 nb.loc[sp, "ci_high"].values - vals],
                 fmt="none", ecolor="black", elinewidth=0.6, capsize=2)
    axB.set_xticks(xb)
    axB.set_xticklabels([S.italic(S.SPECIES_LABEL[s]) for s in sp], fontsize=6,
                        rotation=45, ha="right")
    axB.set_ylim(0, max(0.4, vals.max() * 1.35))
    axB.set_ylabel("Niche breadth (B$_A$)")
    axB.set_title("Host specialisation", fontsize=7.5)
    axB.text(0.5, 0.94, "0 = specialist · 1 = generalist", transform=axB.transAxes,
             ha="center", va="top", fontsize=4.6, color="0.4")
    S.panel_label(axB, "B", dx=-0.30)

    # C: vector-host bipartite feeding network (species left, hosts right; edge width &
    # colour ~ meals). connectance / H2' annotated from host_ecology_indices.
    hosts = [h for h in S.HOST_ORDER if h in vhm.columns]
    species = [s for s in S.SPECIES_HEX if s in vhm.index]
    yS = np.linspace(0.78, 0.22, len(species))           # 3 species
    yH = np.linspace(0.97, 0.03, len(hosts))             # 9 hosts
    xS, xH = 0.10, 0.90
    wmax = float(vhm.values.max())
    cmap = plt.get_cmap(S.SEQ_CMAP)
    for si, s in enumerate(species):
        for hi, h in enumerate(hosts):
            w = float(vhm.loc[s, h])
            if w <= 0:
                continue
            frac = w / wmax
            axC.plot([xS, xH], [yS[si], yH[hi]], "-", lw=0.4 + 2.6 * frac,
                     color=cmap(0.25 + 0.7 * frac), alpha=0.85, zorder=1,
                     solid_capstyle="round")
    # nodes: size ~ total meals through each node
    sp_tot = vhm.loc[species].sum(axis=1)
    h_tot = vhm[hosts].sum(axis=0)
    axC.scatter([xS] * len(species), yS, s=18 + 120 * (sp_tot / sp_tot.max()).values,
                color=[S.SPECIES_HEX[s] for s in species], edgecolor="black",
                linewidth=0.3, zorder=3)
    axC.scatter([xH] * len(hosts), yH, s=10 + 90 * (h_tot / h_tot.max()).values,
                color=[S.host_color(h) for h in hosts], edgecolor="black",
                linewidth=0.3, zorder=3)
    for si, s in enumerate(species):
        axC.text(xS - 0.04, yS[si], S.italic(S.SPECIES_LABEL[s]), ha="right", va="center",
                 fontsize=5.2)
    for hi, h in enumerate(hosts):
        axC.text(xH + 0.04, yH[hi], S.common(h), ha="left", va="center", fontsize=5.2)
    axC.set_xlim(-0.32, 1.30)
    axC.set_ylim(-0.04, 1.04)
    axC.axis("off")
    conn = idx.loc[idx["index"] == "network_connectance", "value"]
    h2 = idx.loc[idx["index"].isin(["network_H2prime_proxy", "network_H2prime"]), "value"]
    sub = []
    if len(conn):
        sub.append(f"connectance = {float(conn.iloc[0]):.2f}")
    if len(h2):
        sub.append(f"H$_2'$ proxy = {float(h2.iloc[0]):.2f}")
    axC.set_title("Vector–host feeding network", fontsize=7.5)
    if sub:
        axC.text(0.5, -0.02, "   ·   ".join(sub), transform=axC.transAxes, ha="center",
                 va="top", fontsize=5, color="0.35")
    S.panel_label(axC, "C", dx=-0.06)

    fig.suptitle("Vector host-use ecology: specialisation, zooprophylaxis and feeding structure",
                 fontweight="bold", fontsize=8, y=1.04)
    S.save(fig, outdir, "figure_5_host_use_ecology")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figure-data", required=True)
    ap.add_argument("--eco-overall", required=True)
    ap.add_argument("--eco-bioclim", required=True)
    ap.add_argument("--geo-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    S.apply_house_style()
    eco_o = pd.read_csv(a.eco_overall, sep="\t")
    eco_b = pd.read_csv(a.eco_bioclim, sep="\t")
    only = {s.strip() for s in a.only.split(",") if s.strip()} or None
    want = lambda f: only is None or f in only
    print(f"main figures -> {a.outdir}")
    if want("1"):
        figure1(a.figure_data, a.geo_dir, a.outdir)
    if want("2"):
        figure2(a.figure_data, a.outdir)
    if want("3"):
        figure3(eco_o, eco_b, a.outdir)
    if want("4"):
        figure4(a.figure_data, eco_b, a.outdir)
    if want("5"):
        figure5(a.figure_data, a.outdir)
    print("done.")


if __name__ == "__main__":
    main()
