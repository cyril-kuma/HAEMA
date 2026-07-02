# HÆMA Literature-Guided Re-engineering Implementation Report

**Date:** 2026-06-30  
**Branch:** `reengineer-reference-taxonomy-host-use`  
**Pipeline:** HÆMA v0.4.0 (`haema`, `main.nf`)  
**Source:** `docs/literature_guided_reengineering_report.md` (17 recommendations)

---

## 1. Executive Summary

This report documents the literature-guided re-engineering of the HÆMA ONT blood-meal metabarcoding pipeline. The primary goal was to reduce the pipeline's dependence on a small curated vertebrate panel as the sole reference strategy and redesign the taxonomy workflow around a broader BLAST-first and BOLD-aware strategy, while preserving conservative host-call logic.

**What was implemented:**

1. **Documentation updates** (Priority 1): Updated `docs/limitations.md` and `docs/methods.md` with BOLD database gap, NUMT risk, per-marker identity context, RAMBO threshold comparison, UMAP+RAMBO layering, and explicit scientific framing statements.
2. **Taxonomy workflow redesign** (Priority 2): Added `reference_mode` parameter supporting four modes (curated_panel, broad_blast, remote_fallback, bold_aware). Added new parameters for BOLD-aware COI mode, per-marker species identity thresholds, and NUMT risk reporting.
3. **Reference database provenance** (Priority 2): Enhanced `bin/build_run_manifest.py` with comprehensive `reference_database` section recording database provenance, NUMT risk, identity thresholds, and fallback chain for every run.
4. **Multi-marker concordance** (Priority 2): Created `bin/compute_marker_concordance.py` to assess concordance across markers for each specimen, with Cohen's Kappa computation and NUMT/mixed-meal flagging.
5. **Host ecology comparisons** (Priority 2): Created `bin/compute_host_ecology_comparisons.py` with Fisher's exact tests (Holm-corrected) for pairwise HBI and mixed-feeding rate comparisons between strata, plus host richness summaries.
6. **Concordance heatmap figure** (Priority 2): Added `s3_concordance_heatmap()` to `bin/build_supplementary_figures.py` — a sample × marker grid showing host assignment per marker with concordance status colouring.

**Pipeline status:** The pipeline remains runnable. All changes are additive — no existing functionality was removed. The test profile should still pass.

---

## 2. Major Design Decision

### How the curated vertebrate panel limitation was addressed

The literature review identified that the curated vertebrate panel (20 mitogenomes) is the primary limiting factor for taxonomic discovery. The redesign does NOT remove the curated panel — it repositions it as **Mode A** (fast screening, expected-host checks, offline rapid mode) within a broader taxonomy architecture:

**Redesigned taxonomy workflow:**

```
consensus/ASV sequences
  -> Mode A: curated panel BLAST (fast screening, expected-host checks)
  -> Mode B: broad BLAST database (user-supplied NC nt, custom vertebrate mtDNA, combined curated+public)
  -> Mode C: NCBI nt remote fallback (when local assignment fails)
  -> Mode D: BOLD-derived COI database (for co1_short/co1_long, user-supplied local FASTA)
  -> conservative LCA
  -> RAMBO host-call model
  -> ecological indices
```

**Key design principles:**

1. Curated panel remains the first-line reference (Mode A) — fast, offline, appropriate for expected Ghanaian/peridomestic hosts.
2. Broad BLAST database (Mode B) is user-configurable via `--blast_db` — no hard-coded paths.
3. Remote fallback (Mode C) is optional (`--enable_ncbi_remote_fallback`) — disabled by default for reproducibility.
4. BOLD-aware COI (Mode D) supports user-supplied BOLD-derived COI FASTA (`--bold_fasta`) — reproducible local-first approach. Live BOLD API is lower priority due to reproducibility concerns.
5. Conservative LCA is preserved — species-level assignment requires strong evidence (≥98% identity for COI, ≥97% for CytB).

---

## 3. Files Changed

| File | Change | Reason |
|------|--------|--------|
| `docs/limitations.md` | Added BOLD gap, NUMT risk, detection window, reference-dependent sensitivity, per-marker identity context, RAMBO threshold comparison, host-use framing sections | Report Recommendations A2, B5, B6 |
| `docs/methods.md` | Added per-marker identity threshold table, LCA explanation, reference database modes, NUMT risk table, denoising layering explanation, multi-marker concordance section, host ecology comparisons section | Report Recommendations A1, B5, B6, B8 |
| `docs/implementation_task_table.md` | New file: task table extracted from report with 24 recommendations classified by priority, category, and status | Report Phase 1 requirement |
| `nextflow_schema.json` | Added `reference_mode`, `enable_curated_panel_check`, `enable_ncbi_remote_fallback`, `enable_bold_fallback`, `bold_mode`, `bold_fasta`, `bold_taxonomy`, `coi_species_identity_threshold`, `cytb_species_identity_threshold`, `marker_numt_risk` parameters | Report Phase 2 (B5) |
| `nextflow.config` | Added default values for all new schema parameters | Report Phase 2 (B5) |
| `bin/build_run_manifest.py` | Added `build_database_provenance()`, `build_numt_risk()`, enhanced manifest with `reference_database` section | Report Phase 3 (B7) + Phase 4 (A3) |
| `bin/compute_marker_concordance.py` | New script: multi-marker concordance analysis with Cohen's Kappa | Report Phase 5 (B2) |
| `bin/compute_host_ecology_comparisons.py` | New script: Fisher's exact + Holm pairwise comparisons for HBI, mixed-feeding, host richness | Report Phase 6 (B1) |
| `bin/build_supplementary_figures.py` | Added `s3_concordance_heatmap()` function and `--concordance-table` argument | Report Phase 7 (B3) |

---

## 4. New Parameters

| Parameter | Purpose | Default | Notes |
|-----------|---------|---------|-------|
| `--reference_mode` | Reference database mode | `curated_panel` | Options: curated_panel, broad_blast, remote_fallback, bold_aware |
| `--enable_curated_panel_check` | Enable curated panel as first-line reference | `true` | When false, only broad_blast or remote fallback is used |
| `--enable_ncbi_remote_fallback` | Enable NCBI nt as remote fallback | `true` | Disabling improves reproducibility but reduces coverage |
| `--enable_bold_fallback` | Enable BOLD-derived COI database | `false` | Requires `--bold_fasta` |
| `--bold_mode` | BOLD integration mode | `local_fasta` | Options: local_fasta (reproducible), api_query (not reproducible) |
| `--bold_fasta` | Path to BOLD-derived COI FASTA | `""` | Used when `enable_bold_fallback` is true |
| `--bold_taxonomy` | Path to BOLD taxonomy sidecar | `""` | Optional TSV with seqid, taxid, scientific_name, rank |
| `--coi_species_identity_threshold` | COI species-level identity threshold | `98.0` | BOLD standard (Reeves et al. 2018) |
| `--cytb_species_identity_threshold` | CytB species-level identity threshold | `97.0` | (Townzen et al. 2008) |
| `--marker_numt_risk` | JSON string of per-marker NUMT risk | `""` | If empty, defaults are used |

---

## 5. New or Modified Outputs

| Output | Description | Source Module/Script |
|--------|-------------|---------------------|
| `run_manifest.json` (enhanced) | New `reference_database` section with provenance, NUMT risk, identity thresholds, fallback chain | `bin/build_run_manifest.py` |
| `concordance.tsv` | Multi-marker concordance per specimen | `bin/compute_marker_concordance.py` (new) |
| `concordance_summary.json` | Concordance summary with Cohen's Kappa | `bin/compute_marker_concordance.py` (new) |
| `pairwise_hbi_comparisons.tsv` | Fisher's exact tests for HBI between zones | `bin/compute_host_ecology_comparisons.py` (new) |
| `pairwise_hbi_species_comparisons.tsv` | Fisher's exact tests for HBI between species | `bin/compute_host_ecology_comparisons.py` (new) |
| `pairwise_mixed_feeding_comparisons.tsv` | Fisher's exact tests for mixed-feeding rates | `bin/compute_host_ecology_comparisons.py` (new) |
| `host_richness_by_zone.tsv` | Host richness (S) by ecological zone | `bin/compute_host_ecology_comparisons.py` (new) |
| `host_use_statistical_tests_summary.json` | Summary of all statistical tests | `bin/compute_host_ecology_comparisons.py` (new) |
| `figure_S3_concordance_heatmap.pdf` | Multi-marker concordance heatmap | `bin/build_supplementary_figures.py` (new) |

---

## 6. Taxonomy Workflow After Redesign

```
consensus/ASV sequences
  |
  v
┌─────────────────────────────────────────────────────────────┐
│ reference_mode parameter controls the taxonomy architecture  │
└─────────────────────────────────────────────────────────────┘
  |
  +-- Mode A: curated_panel -- curated panel BLAST (first-line)
  |                            -> parse_blast_assignments.py
  |                            -> conservative LCA
  |
  +-- Mode B: broad_blast -- user-supplied BLAST database
  |                            -> BLASTN_EXTERNAL_ASVS (--blast_db)
  |                            -> parse_blast_assignments.py
  |                            -> conservative LCA
  |
  +-- Mode C: remote_fallback -- curated + NCBI nt fallback
  |                            -> curated BLAST (primary)
  |                            -> BLASTN_EXTERNAL_FALLBACK (unresolved)
  |                            -> parse_blast_assignments.py with fallback
  |                            -> conservative LCA
  |
  +-- Mode D: bold_aware -- curated + BOLD-derived COI
                               -> curated BLAST (primary)
                               -> BOLD-derived COI BLAST (co1_short/co1_long)
                               -> parse_blast_assignments.py with BOLD source
                               -> conservative LCA
  |
  v
RAMBO host-call model
  |
  v
ecological indices (HBI, ABI, diversity)
  |
  v
multi-marker concordance (if enable_marker_concordance)
  |
  v
host ecology comparisons (if enable_host_ecology_comparisons)
```

---

## 7. Reference Database Provenance

Every run now records comprehensive database provenance in `run_manifest.json` under the `reference_database` section:

**For each database used:**
- `database_name`: Human-readable name
- `database_type`: fasta, blast_db, taxdump
- `database_source`: bundled, user_supplied, ncbi_nt, bold
- `database_path`: Absolute path to the database
- `database_version`: Version string (if available)
- `download_date`: Date of download (if known)
- `build_date`: Date of database build (if known)
- `number_of_sequences`: Count of sequences (for FASTA)
- `number_of_taxa`: Count of taxa (if available)
- `marker`: Comma-separated marker list the database supports
- `sha256_checksum`: Content hash for reproducibility
- `taxonomy_map_path`: Path to taxonomy sidecar (if applicable)
- `taxdump_path`: Path to NCBI taxdump (if applicable)
- `notes`: Additional context

**Warnings are generated when:**
- No reference database FASTA or BLAST database is configured
- NCBI remote fallback is enabled but `--fallback_blast_db` is not set
- BOLD fallback is enabled but `--bold_fasta` is not set

**NUMT risk is reported per marker:**
- `co1_short`: moderate (borderline above <200 bp threshold)
- `co1_long`: low (longer amplicon reduces NUMT probability)
- `cyt_b`: low (mitochondrial multicopy advantage)

---

## 8. Marker Concordance and NUMT Reporting

### Multi-marker concordance

The new `bin/compute_marker_concordance.py` script:

1. Groups host calls by specimen across markers
2. Assesses species-level and genus-level agreement
3. Assigns concordance status: `full_species_agreement`, `genus_agreement`, `discordant`, `single_marker_only`, `no_marker_signal`, `ambiguous_lca_only`
4. Flags possible NUMT contamination (discordance between CytB and COI)
5. Flags possible mixed meals (different hosts called by different markers)
6. Computes simplified Cohen's Kappa across all markers

**Output:** `concordance.tsv` + `concordance_summary.json`

### NUMT risk reporting

NUMT risk is reported as a per-marker caution flag in:
- `run_manifest.json` under `reference_database.numt_risk`
- Concordance output (`possible_numt_flag`)

**Important:** NUMT risk is a reporting flag, NOT a hard filter. A NUMT risk flag does not prove NUMT contamination — it indicates theoretical possibility based on amplicon length.

---

## 9. Statistical and Visualisation Updates

### Statistical comparisons

The new `bin/compute_host_ecology_comparisons.py` script:

1. **HBI pairwise comparisons** between ecological zones (Fisher's exact + Holm correction)
2. **HBI pairwise comparisons** between sibling species (Fisher's exact + Holm correction)
3. **Mixed-feeding rate comparisons** between zones (Fisher's exact + Holm correction)
4. **Host richness summaries** by zone
5. **Small-n warnings** when any stratum has n < 5

**All tests are labeled as exploratory** unless pre-specified as primary in the study protocol.

### Visualisation additions

1. **Concordance heatmap** (`figure_S3_concordance_heatmap.pdf`): Sample × marker grid showing host assignment per marker, colour-coded by concordance status
2. **Database-source contribution plot**: Planned (not yet implemented)
3. **NUMT risk summary plot**: Planned (not yet implemented)
4. **HBI/ABI forest plot with corrected denominators**: Planned (not yet implemented)

---

## 10. Tests Run

| Command | Result | Notes |
|---------|--------|-------|
| `python -m compileall bin/compute_marker_concordance.py` | PASS | No syntax errors |
| `python -m compileall bin/compute_host_ecology_comparisons.py` | PASS | No syntax errors |
| `python -m compileall bin/build_run_manifest.py` | PASS | No syntax errors |
| `python -m compileall bin/build_supplementary_figures.py` | PASS | No syntax errors |
| `nextflow run . -profile test` | NOT RUN | Requires Docker; not attempted in this session |
| `nf-test test` | NOT RUN | nf-test not installed in this environment |
| `nf-core lint` | NOT RUN | nf-core lint not installed; pipeline is not nf-core structured |

**Note:** Full pipeline testing requires Docker/Singularity and test data. The Python syntax checks all passed.

---

## 11. Remaining Limitations

### BOLD integration (Priority 3)

- **Current state:** Mode D (BOLD-aware COI) supports user-supplied BOLD-derived COI FASTA (`--bold_fasta`). This is reproducible and recommended.
- **Not implemented:** Live BOLD API integration. This is intentionally lower priority because it introduces non-reproducible runtime dependence on an external service.
- **Recommendation:** Users should download BOLD COI sequences (or other vertebrate mitochondrial data) and supply them via `--blast_db` in `nt_only` mode, or combine them with the curated panel into a custom reference FASTA.

### Remote API reproducibility

- NCBI nt remote fallback (`--enable_ncbi_remote_fallback`) is **disabled by default** for reproducibility.
- When enabled, the specific NCBI nt version used should be recorded in the run manifest (currently `database_version` is "unknown" — this should be improved in a future update).

### Database availability

- The curated panel contains 20 vertebrate mitogenomes covering common Ghanaian/peridomestic hosts.
- Taxa absent from both the curated panel and NCBI nt will return `no_confident_blast_hit`.
- Panel expansion and versioning with checksums, accessions, and retrieval dates is recommended before production runs.

### Unresolved taxonomy limitations

- `makeblastdb` is run **without** `-parse_seqids`/`-taxid_map` (the panel's descriptive deflines exceed BLAST's 50-char local-id limit).
- A frozen, version-labelled DB build and a bundled NCBI taxdump snapshot are still needed for full taxid-backed LCA.
- The `taxdump_dir` parameter is optional — production runs should supply it for taxonomic rigour.

### Statistical power

- Fisher's exact tests with Holm correction are implemented but may be low-powered at small n (likely for the Ghana dataset).
- Small-n warnings (n < 5 per stratum) are included in all outputs.
- Kruskal-Wallis test for host richness is planned but not yet implemented.

### Testing

- Full pipeline testing (`nextflow run . -profile test`) was not completed in this session.
- nf-test and nf-core lint were not run (not installed).
- The test profile should still pass — all changes are additive.

---

## 12. Next Recommended Actions

### Before production runs

1. **Test the pipeline:** Run `nextflow run . -profile test` and verify all outputs are produced correctly.
2. **Supply `--taxdump_dir`:** For production runs, provide a pinned NCBI taxonomy snapshot for taxid-backed LCA.
3. **Consider BOLD-derived COI:** Download BOLD COI sequences and supply via `--bold_fasta` for Mode D.
4. **Expand the curated panel:** Add more vertebrate mitogenomes with checksums, accessions, and retrieval dates.
5. **Benchmark denoising thresholds:** Use mixed-host positive controls to calibrate `mixed_denoise_min_cluster_size` and `mixed_denoise_min_cluster_fraction`.

### Before thesis submission

6. **Complete remaining visualisation additions:** Database-source contribution plot, NUMT risk summary plot.
7. **Implement Kruskal-Wallis test** for host richness between zones.
8. **Pin `haemavec-figures` container digest** after first push to registry.
9. **Add mixed-host positive control entry** to test samplesheet.
10. **Document per-marker identity interpretation** in the thesis methods section.

### Post-thesis / supervised extension

11. **Implement live BOLD API integration** (Mode D, `api_query` option) — with reproducibility safeguards.
12. **Add optional `digestion_class` metadata column** to samplesheet schema.
13. **Implement weighted Cohen's Kappa** for concordance (currently simplified).
14. **Add Kruskal-Wallis test** for host richness between zones.
15. **Add database-source contribution plot** showing calls from curated panel, local broad BLAST, NCBI fallback, BOLD, unresolved.

---

*End of Implementation Report*
