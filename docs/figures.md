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
├── figure_captions.md      # draft captions, naming the input file(s) for each figure
└── figure_manifest.tsv     # figure → files → inputs → caption (machine-readable)
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

Full draft captions (with the exact input files) are written to `07_figures/figure_captions.md` on
every run.

## How to (re)generate

**As part of the pipeline** (default). `BUILD_FIGURES` runs after the endpoint/report steps; it only
consumes tables that exist, so figures needing host calls (5, 6) or R outputs (7) are skipped if
`--enable_rambo_model` / `--enable_r_outputs` are off. Toggle and configure with:

```bash
--enable_figures true|false      # default true (off in -profile test)
--figure_formats pdf,svg,png     # any subset
--figures_container haema-figures:0.3.0
```

**Standalone**, against an existing results directory (fast iteration, no full pipeline run):

```bash
docker run --rm -u $(id -u):$(id -g) -v "$PWD":"$PWD" -w "$PWD" haema-figures:0.3.0 \
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

Figures use the project-specific **`haema-figures:0.3.0`** image (matplotlib + seaborn + pandas on a
digest-pinned `python:3.11-slim` base — the same base as `haema-python`). Build it once, like the
other two custom images:

```bash
docker build -t haema-figures:0.3.0 -f containers/haema-figures/Dockerfile .
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
