# HÆMA architecture review (`bin/` + `modules/local/`)

PI-led functional review, 2026-07-02. Method: cross-referenced every `bin/` script against
`modules/local/*.nf`, `main.nf`, other `bin/` scripts (imports), `tests/` and the `Makefile`.

## Headline finding

The tracked codebase is **already well-organised**: nearly every script maps 1:1 to a single
Nextflow process, shared logic is factored into helper libraries, there are **no hard-coded
absolute paths in reusable code**, and there is **no dead/uncalled tracked script** — only
standalone CLI tools with legitimate provenance/verification roles. The one genuine
consolidation opportunity is the two parallel figure systems (below). Aggressive merging was
**deliberately not performed** because it would risk a working figure/reporting stack that
cannot be fully rendered-tested in this environment, and the prompt's own guidance is not to
merge computationally distinct, separately-cached stages for aesthetics.

## Workflow map (stage → module → script)

| Stage | Module | Script | Status |
|-------|--------|--------|--------|
| Ingest | `fastq_ingest` | `merge_fastqs.py` (+`fastq_utils`) | keep |
| Validate | `input_validation` | `validate_inputs.py` | keep |
| Preflight | `production_preflight`, `medaka_preflight`, `r_preflight` | `production_preflight.py` | keep |
| Demux (opt) | `advanced_demux` | `demux_fastq_by_header.py` | keep |
| Preprocess | `trim_filter_split` | `trim_filter_split.py` | keep |
| Denoise (opt) | `denoise_mixed_templates` | `denoise_mixed_templates.py` | keep |
| Consensus/ASV | `cluster_consensus`, `dereplicate_asvs` | `cluster_consensus.py`, `dereplicate_asvs.py` | keep |
| Polish (opt) | `medaka` | (medaka) | keep |
| Taxonomy | `taxonomy` subwf → `blast/*`, `taxonomy_assign`, `curated_reference_metadata`, `no_taxonomy` | `parse_blast_assignments.py`, `build_curated_taxid_map.py`, `no_taxonomy.py` | **redesigned** (Modes A–D) |
| Aggregate | `aggregate_results` | `aggregate_results.py` | keep |
| Host model | `rambo_model` | `rambo_mixed_model.py` | keep |
| Ecology | `ecological_indices` | `compute_ecological_indices.py` | keep |
| Concordance | `marker_concordance` | `compute_marker_concordance.py` | **new, wired** |
| Comparisons | `host_ecology_comparisons` | `compute_host_ecology_comparisons.py` | **new, wired** |
| Report | `report` | `build_report.py` | keep |
| R outputs | `r_outputs` | `build_r_outputs.R`, `build_phyloseq_figures.R` | keep |
| MultiQC | `multiqc` | (multiqc) | keep |
| Figures (auto) | `figures` | `build_figures.py` | keep (see consolidation) |
| Figures (pub) | `publication_figures` | `build_main_figures.py`, `build_supplementary_figures.py`, `figure_data_prep.py`, `compute_host_ecology_indices.py` (+`figure_style`) | keep (see consolidation) |
| Manifest | `run_manifest` | `build_run_manifest.py` | keep (provenance extended) |

## `bin/` decisions

| Script | Role | Decision | Reason |
|--------|------|----------|--------|
| `fastq_utils.py` | shared FASTQ helpers (imported by 6 scripts) | **keep (library)** | Distinct shared dependency; correctly factored |
| `figure_style.py` | shared figure styling (imported by figure scripts) | **keep (library)** | Colourblind-safe palette/style shared by figure scripts |
| `parse_blast_assignments.py` | conservative LCA + marker guard | **keep (core)** | Scientific boundary; unit-tested |
| `build_reference_db.py` | broad reference-DB builder (new) | **keep (CLI tool)** | Reproducible DB provenance; standalone by design |
| `build_reference_panel.py` | curated-panel builder | **keep (CLI tool)** | Documents how the curated panel was built; not pipeline-called by design |
| `verify_reference_assets.py` | asset checksum/sidecar check | **keep (tooling)** | Used by `validate_release.sh` + `Makefile` |
| all other scripts | 1:1 with a process | **keep** | Each is a distinct, separately-cached, testable stage |
| `build_thesis_figures.py`, `figure3B_options.py` | **untracked** user scratch | **leave untouched** | Not committed, not wired; the user's own analysis scripts — not repo debris to delete |

No tracked script is uncalled-and-purposeless; nothing merged or deleted.

## `modules/local/` decisions

All 30 modules represent genuine Nextflow process boundaries (distinct tool/container/resource
or scientific meaning) and are retained. The taxonomy subworkflow was redesigned (Modes A–D) and
two additive modules (`marker_concordance`, `host_ecology_comparisons`) were wired in. No modules
merged: the DAG boundaries aid caching/scheduling and reproducibility.

## Consolidation recommendation (not executed)

Two figure systems coexist: `figures`/`build_figures.py` (1066 lines, "automated" figures) and
`publication_figures`/`build_main_figures.py`+`build_supplementary_figures.py` (curated thesis
suite). They overlap (host spectrum, HBI/ABI, mixed feeding, QC). **Recommendation:** designate
the `publication_figures` suite as canonical for the thesis/manuscript, fold the unique
`build_figures.py` panels into it (or the shared `figure_data_prep.py`), and retire the
`figures` module. **Deferred**, not executed, because figure rendering (geopandas/matplotlib) can
only be validated by running the containerised figure steps, and the current dual system is
working; a blind merge risks regressions. Tracked as a release-time task.

## Stale-file audit

No stale tracked files were found: every `docs/*.md` is current (methods/usage/parameters/output
updated this cycle; implementation_report + task_table corrected). The only working-tree scratch
files (`build_thesis_figures.py`, `figure3B_options.py`, `run_full.sh` variants) are **untracked**
user artifacts and were left untouched. No deletions were required.

## No hard-coded paths

`grep` over `main.nf`, `nextflow.config`, `modules/`, `subworkflows/`, `bin/` found no absolute
machine paths. Runtime paths are supplied via `--input`, `--raw_data_dir`, `--blast_db`,
`--blast_db_mount`, `--outdir` (the `run_*.sh` launchers carry machine-specific paths by design,
as run helpers, not reusable code).
