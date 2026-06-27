# Manuscript figures

HÆMA generates a set of **publication-ready figures** directly from its endpoint tables, so the
panels in a manuscript trace back to the same files a reviewer can download. Figures are produced by
the `BUILD_FIGURES` step (on by default) and written to:

```text
results/<run>/07_figures/
├── figure_01_workflow.{pdf,svg,png}
├── figure_02_sequencing_qc.{pdf,svg,png}
├── figure_03_depth_denoising.{pdf,svg,png}
├── figure_04_host_assignment.{pdf,svg,png}
├── figure_05_host_composition.{pdf,svg,png}
├── figure_06_mixed_host.{pdf,svg,png}
├── figure_07_controls_contamination.{pdf,svg,png}
├── figure_08_ecology.{pdf,svg,png}
├── figure_09_ecological_indices.{pdf,svg,png}
├── figure_10_temporal.{pdf,svg,png}
├── figure_11_phyloseq_composition.{pdf,png}     # phyloseq-native (R) — PDF + PNG
├── figure_12_phyloseq_alpha_diversity.{pdf,png}
├── figure_13_phyloseq_ordination.{pdf,png}
├── figure_14_phyloseq_heatmap.{pdf,png}
├── figure_15_phyloseq_decontam.{pdf,png}
├── figure_captions.md            # captions for the Python figures (01–10)
├── figure_manifest.tsv
├── phyloseq_figure_captions.md   # captions for the phyloseq figures (11–15)
└── phyloseq_figure_manifest.tsv
```

Each figure is exported as **vector PDF and SVG** (editable text — `svg.fonttype=none`,
`pdf.fonttype=42`) plus a **300 dpi PNG** for quick viewing and slides. Colours use the
colour-blind-safe Okabe–Ito palette, with a fixed host-taxon colour map so a species is the same
colour in every figure.

> Every figure is built **only from real pipeline outputs** — nothing is simulated or hand-drawn
> except the workflow schematic (Figure 1), which depicts the pipeline architecture. Where a value is
> evidence rather than a validated estimate (host fractions, mixed-host thresholds) the figure and
> its caption say so.

## What each figure shows

| Figure | Shows | Built from |
|---|---|---|
| **1 · Workflow** | The 9-stage pipeline architecture (methods schematic). | pipeline architecture (no data) |
| **2 · Sequencing & QC** | Read-processing funnel (raw→Q20→marker-assigned→retained); per-sample mean read quality; consensus length per marker (the 3 amplicon size classes); reads assigned per marker by primer vs length. | `qc_summary.tsv`, `bloodmeal_master_endpoint.tsv` |
| **3 · Depth & denoising** | Per-sample retained depth (log, by control class); retained ASVs vs assigned hosts; UMAP/HDBSCAN vs greedy-fallback usage per marker. | `sample_level_summary.tsv`, `qc_summary.tsv` |
| **4 · Host assignment** | Host-taxon frequency across field samples; host × marker recovery heatmap; BLAST identity vs coverage for accepted features with the 97%/80% thresholds. | `host_call_table.tsv`, `bloodmeal_master_endpoint.tsv` |
| **5 · Host composition** | Per-sample stacked host read fractions, one panel per marker; bars compose to 1.0 with a grey *unresolved* remainder; mixed feeds flagged (▲). | `host_call_table.tsv` |
| **6 · Mixed-host feeding** | Feeding-type composition (single/mixed/none); each detected mixed blood meal as its co-host fractions; mixed-host control sensitivity annotation. | `host_call_table.tsv`, `rambo_model_summary.tsv` |
| **7 · Controls & contamination** | Positive-control host recovery (expected vs recovered, status); phyloseq/decontam summary and per-marker negative-control background. | `positive_control_check.tsv`, `contamination_model_summary.tsv`, `qc_background_thresholds.tsv` |
| **8 · Ecology** | Host detections stratified by ecological zone and by *Anopheles gambiae* s.l. sibling species (per-group n shown). | `bloodmeal_master_endpoint.tsv` |
| **9 · Ecological indices** | Human Blood Index & zoophily forest plot (Wilson 95% CI), feeding-type partition, and host diversity, stratified by zone and sibling species. | `ecological_indices.tsv` |
| **10 · Temporal** | Sampling timeline, HBI & mixed-feeding per collection campaign (Wilson CI), and feeding-type composition over time — **descriptive only** (campaigns confounded with site/batch). | `ecological_indices.tsv` |

### phyloseq-native figures (R, `PHYLOSEQ_FIGURES`)
Built directly from `bloodmeal_phyloseq.rds` with **phyloseq + ggplot2** (vector PDF + 300 dpi PNG;
SVG is only produced for the Python figures). ASVs are agglomerated to `host_assignment`; controls
are excluded from the ecological views. This step runs in the R container and **no-ops cleanly** if
the object is an R-fallback (non-phyloseq) `.rds`.

| Figure | Shows | Built from |
|---|---|---|
| **11 · Composition** | `plot_bar` relative host abundance per field sample, faceted by ecological zone. | `bloodmeal_phyloseq.rds` |
| **12 · Alpha diversity** | `plot_richness` Observed + Shannon by zone — "host breadth" over ≤5 hosts (low-resolution; hosts only). | `bloodmeal_phyloseq.rds` |
| **13 · Beta diversity** | `plot_ordination` PCoA (Bray–Curtis) by zone — **descriptive**, separation reflects majority host, not a gradient (skips if degenerate). | `bloodmeal_phyloseq.rds` |
| **14 · Heatmap** | `plot_heatmap` host × sample read abundance (robust low-diversity companion to ordination). | `bloodmeal_phyloseq.rds` |
| **15 · decontam** | Per-ASV prevalence in negative controls vs field samples, contaminant flagged — demonstrates the decontaminated object. | `bloodmeal_phyloseq.rds`, `decontam_results.tsv` |

Full draft captions (with the exact input files) are written to `07_figures/figure_captions.md` on
every run.

## Curated Objective 1 publication suite

The `PUBLICATION_FIGURES` step (enabled with `--enable_publication_figures true`; on by default
outside the `test` profile) writes a streamlined manuscript/thesis suite to
`--publication_figures_dir` (default `results/figures/`):

```text
results/figures/
├── main/              figure_1_study_system.{pdf,eps} ... figure_5_host_use_ecology.{pdf,eps}
├── supplementary/     figure_S1_rarefaction.{pdf,eps}, figure_S2_*.{pdf,eps}, figure_B1_*.{pdf,eps}
└── figure_data/       plot-ready TSVs, vector_host_matrix.tsv, host_ecology_indices.tsv
```

This suite adds the Ghana GADM bioclimatic-zone map and vector-host ecology panels (Levins niche
breadth, Pianka overlap, Bray-Curtis turnover, network connectance / H2prime proxy). It needs the
geo-enabled `haema-figures:0.4.0` image.

## How to (re)generate

**As part of the pipeline** (default). `BUILD_FIGURES` runs after the endpoint/report steps; it only
consumes tables that exist, so figures needing host calls (5, 6) or R outputs (7) are skipped if
`--enable_rambo_model` / `--enable_r_outputs` are off. Toggle and configure with:

```bash
--enable_figures true|false      # default true (off in -profile test)
--figure_formats pdf,svg,png     # any subset
--enable_publication_figures true|false
--publication_figures_dir results/figures
--figures_container haema-figures:0.4.0
```

**Standalone**, against an existing results directory (fast iteration, no full pipeline run):

```bash
docker run --rm -u $(id -u):$(id -g) -v "$PWD":"$PWD" -w "$PWD" haema-figures:0.4.0 \
  python pipeline/bin/build_figures.py \
    --endpoint-dir results/<run>/05_endpoint_files \
    --reports-dir  results/<run>/06_reports \
    --manifest     results/<run>/05_endpoint_files/run_manifest.json \
    --outdir       results/<run>/07_figures \
    --formats pdf,svg,png
```

The script ([`bin/build_figures.py`](../bin/build_figures.py)) finds each table by its canonical
filename and **degrades gracefully**: a missing table skips its figure, and the step only fails if
*no* figure can be produced.

## Container

Figures use the project-specific **`haema-figures:0.4.0`** image (matplotlib + seaborn + geopandas
on a digest-pinned `python:3.11-slim` base — the same base as `haema-python`). Build it once, like
the other two custom images:

```bash
docker build -t haema-figures:0.4.0 -f containers/haema-figures/Dockerfile .
```

For HPC/publication, push it to a registry and pass its immutable `@sha256:` digest via
`--figures_container`.

## Notes & caveats

- **Host fractions are evidence summaries, not validated quantitative estimates**, and mixed-host
  detection thresholds are **not yet benchmarked** — Figures 5 and 6 state this; see
  [`limitations.md`](limitations.md).
- The grey *unresolved* segment in Figures 5–6 is the read fraction that was noise-filtered during
  denoising or left taxonomically unassigned; bars therefore compose to 1.0 honestly.
- Figure 8 describes the demonstration cohort (host detections by zone/vector); it is **not** a
  formal association test.
- Figures inherit the provenance of the endpoint tables they are built from. Regenerate them after
  any re-run so they match the current samplesheet and parameters.
