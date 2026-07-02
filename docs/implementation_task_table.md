# HÆMA Literature-Guided Re-engineering — Implementation Task Table

Generated: 2026-06-30 · **Statuses corrected: 2026-07-02**
Source: docs/literature_guided_reengineering_report.md

Status legend: **DONE** (implemented, wired and verified to the extent possible offline) ·
**SCAFFOLD+DOC** (config/docs present, deeper implementation is documented future work) ·
**DEFERRED** (post-thesis) · **NOT STARTED**.

> Verification note: the pipeline cannot be run end-to-end in the dev environment (BLAST/Medaka
> live only inside containers; no local `nt`/BOLD databases). "DONE" items are verified by Python
> unit tests, functional script runs in a venv, `nextflow config`, and `nextflow run -preview -stub`
> DAG builds for every reference mode — not by a full containerised run. See docs/implementation_report.md §10.

## Priority 1 — before thesis submission

| # | Recommendation | Type | Status |
|---|---------------|------|--------|
| A1 | Update docs/methods.md (framing, per-marker identity, RAMBO threshold, UMAP+RAMBO layering) | Documentation | **DONE** |
| A2 | Update docs/limitations.md (BOLD gap, NUMT, detection window, reference dependency) | Documentation | **DONE** |
| A3 | Add per-marker NUMT risk to run_manifest.json | Code | **DONE** (wired via parameters.json + verified) |
| A4 | Add mixed-host positive control entry to test samplesheet | Testing | **NOT STARTED** (deferred; needs test-data curation) |

## Priority 2 — before viva

| # | Recommendation | Type | Status |
|---|---------------|------|--------|
| B1 | Fisher's exact + Holm HBI/mixed-feeding comparisons | Code + Nextflow module | **DONE** (`compute_host_ecology_comparisons.py` fixed to real schema + `HOST_ECOLOGY_COMPARISONS` module + tests) |
| B2 | Multi-marker concordance analysis | Code + Nextflow module | **DONE** (`compute_marker_concordance.py` fixed to real schema + `MARKER_CONCORDANCE` module + tests) |
| B3 | Concordance heatmap figure | Visualisation | **DONE** (figure S3 wired into `PUBLICATION_FIGURES` via optional `--concordance-table`) |
| B4 | Pin haemavec-figures container digest | Config | **NOT STARTED** (needs a registry push first; see report §12) |
| B5 | Taxonomy redesign: Modes A (curated) / B (broad) / C (remote) / D (BOLD) | Nextflow + Code | **DONE** (`reference_mode` canonical, wired; all 4 modes DAG-verified) |
| B6 | Per-marker identity interpretation | Code + Documentation | **DONE** (enforced as a conservative confidence downgrade + documented) |
| B7 | Reference database provenance reporting | Code | **DONE** (run_manifest `reference_database`: mode, per-DB SHA-256, fallback chain, thresholds) |
| B8 | QC matrix documentation | Documentation | **DONE** (in limitations/methods; report §10 QC matrix) |

## Priority 3 — post-thesis / supervised extension

| # | Recommendation | Type | Status |
|---|---------------|------|--------|
| C1 | Live BOLD API tertiary fallback | Code + Nextflow | **SCAFFOLD+DOC** (`bold_mode=api_query` reserved; local BOLD FASTA implemented instead for reproducibility) |
| C2 | Optional `digestion_class` metadata column | Config + Code | **DEFERRED** |
| C3 | Kruskal-Wallis test for host richness | Code | **DEFERRED** (host-richness table emitted; formal test not yet added) |

## Documentation-only

| # | Recommendation | Status |
|---|---------------|--------|
| D1 | Length-fallback marker-misassignment risk | **DONE** |
| D2 | UMAP denoising vs RAMBO layering | **DONE** |
| D3 | RAMBO 1% vs Logue 10% threshold | **DONE** |
| D4 | Per-marker identity interpretation (97% global, 98% BOLD COI) | **DONE** |

## Not in the original report but required by "make honest"

| Item | Status |
|------|--------|
| Restore test suite from assets/tests/ back to tests/ | **DONE** |
| Fix concordance/ecology scripts that read non-existent columns (would silently mis-report) | **DONE** (+ new unit tests) |
| Reconcile the redundant `reference_mode` vs `taxonomy_strategy` config surfaces | **DONE** (reference_mode canonical; taxonomy_strategy an alias) |
| Correct the over-claiming implementation report + stale task table | **DONE** (this file + implementation_report.md) |
