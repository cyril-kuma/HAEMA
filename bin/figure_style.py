"""
HÆMA publication-figure global style (Objective 1).

Implements the agreed house style: sans-serif (Helvetica/Arial, DejaVu fallback), base 6 pt /
axis-title 7 pt / panel-label 8 pt bold; no top/right spines; pure-white canvas; no heavy
gridlines; direct labels over legends where possible; NO rotated axis text; colour-blind-safe
palettes with FIXED bioclimatic-zone hex codes; vector (PDF + EPS) export.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- fixed bioclimatic-zone palette + canonical labels/order (south -> north) ----------
ZONE_HEX = {"Coastal_Savannah": "#4477AA", "Forest": "#117733", "Northern_Savannah": "#DDAA33"}
ZONE_LABEL = {"Coastal_Savannah": "Coastal Savannah", "Forest": "Forest",
              "Northern_Savannah": "Northern Savannah"}
ZONE_ORDER = ["Coastal_Savannah", "Forest", "Northern_Savannah"]
# Zones with too few host-identified meals for inferential reporting (descriptive only).
ZONE_DESCRIPTIVE = {"Coastal_Savannah"}
# GADM level-1 region -> study bioclimatic zone (documented 3-zone collapse for Ghana)
REGION_ZONE = {
    "GreaterAccra": "Coastal_Savannah", "Central": "Coastal_Savannah",
    "Ashanti": "Forest", "Ahafo": "Forest", "Bono": "Forest", "BonoEast": "Forest",
    "Eastern": "Forest", "Western": "Forest", "WesternNorth": "Forest", "Oti": "Forest",
    "Volta": "Forest",
    "Northern": "Northern_Savannah", "NorthEast": "Northern_Savannah",
    "Savannah": "Northern_Savannah", "UpperEast": "Northern_Savannah",
    "UpperWest": "Northern_Savannah",
}

# ---- colour-blind-safe host key (Okabe-Ito + Tol), distinct roles from the zone palette ----
HOST_HEX = {
    "Homo sapiens": "#0072B2", "Bos taurus": "#D55E00", "Ovis aries": "#009E73",
    "Capra hircus": "#E69F00", "Canis lupus familiaris": "#CC79A7", "Sus scrofa": "#56B4E9",
    "Equus asinus": "#F0E442", "Gallus gallus": "#882255", "Felis catus": "#999933",
    "unassigned": "#BBBBBB",
}
HOST_COMMON = {"Homo sapiens": "Human", "Bos taurus": "Cattle", "Ovis aries": "Sheep",
               "Capra hircus": "Goat", "Canis lupus familiaris": "Dog", "Sus scrofa": "Pig",
               "Equus asinus": "Donkey", "Gallus gallus": "Chicken", "Felis catus": "Cat"}
HOST_ORDER = ["Homo sapiens", "Bos taurus", "Ovis aries", "Capra hircus",
              "Canis lupus familiaris", "Sus scrofa", "Equus asinus", "Gallus gallus",
              "Felis catus"]
# sibling-species palette (distinct from the zone palette, esp. Guinea #DDAA33, and host blues).
# Keys match the samplesheet / endpoint convention used elsewhere in the pipeline.
SPECIES_HEX = {
    "Anopheles_coluzzii": "#004488",
    "Anopheles_gambiae_s.s": "#BB5566",
    "Anopheles_arabiensis": "#44AA99",
}
SPECIES_LABEL = {
    "Anopheles_coluzzii": "An. coluzzii",
    "Anopheles_gambiae_s.s": "An. gambiae s.s.",
    "Anopheles_arabiensis": "An. arabiensis",
}

MARKER_LABEL = {"cyt_b": "cyt b", "co1_short": "short COI", "co1_long": "long COI"}

# single sequential, perceptually-ordered, colour-blind-safe colourmap for all
# density/heatmap panels (Fig 2C hexbin, Fig 5 network, Appendix overlap/turnover).
SEQ_CMAP = "YlGnBu"


def apply_house_style():
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "font.size": 6, "axes.titlesize": 7, "axes.labelsize": 7,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
        "axes.titleweight": "bold",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.6, "axes.grid": False,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "legend.frameon": False, "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "savefig.bbox": "tight", "lines.linewidth": 1.0,
    })


def common(sci):
    return HOST_COMMON.get(sci, sci)


def italic(sci):
    return r"$\it{" + str(sci).replace(" ", r"\ ") + r"}$"


def two_line(sci):
    c = HOST_COMMON.get(sci)
    return (f"{c}\n" + italic(sci)) if c else sci


def zone_label(z):
    return ZONE_LABEL.get(z, z)


def zone_color(z):
    return ZONE_HEX.get(z, "#888888")


def host_color(sci):
    return HOST_HEX.get(sci, "#555555")


def panel_label(ax, s, dx=-0.13, dy=1.06):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=8, fontweight="bold",
            va="top", ha="left")


def save(fig, outdir, name, preview=True):
    """Vector PDF + EPS into outdir; optional PNG preview into outdir/previews/."""
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, f"{name}.pdf"))
    try:
        fig.savefig(os.path.join(outdir, f"{name}.eps"))
    except Exception as e:  # EPS lacks alpha/transparency; non-fatal
        print(f"    (eps skipped for {name}: {e})")
    if preview:
        pv = os.path.join(outdir, "previews")
        os.makedirs(pv, exist_ok=True)
        fig.savefig(os.path.join(pv, f"{name}.png"), dpi=300)
    plt.close(fig)
    print(f"  wrote {name}.pdf (+eps, preview png)")
