# Changelog

All notable changes to the HÆMA pipeline are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project aims to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- **qc_summary OOM (exit 137).** The aggregated `qc_summary.tsv` no longer concatenates the per-read
  denoise `cluster_membership` and trim `read_decisions` tables (millions of rows on real data — it
  reached ~1 GB and OOM-killed `AGGREGATE_RESULTS`/`BUILD_R_OUTPUTS`). Only per-sample/per-marker
  summaries are aggregated; the per-read tables remain published under `02_trimmed_filtered/qc/` and
  `03_consensus_variants/mixed_denoising/qc/` (`main.nf`, `subworkflows/local/preprocess/main.nf`).

### Changed
- **MultiQC report is now real.** The MultiQC step exposes the per-sample and per-marker summaries as
  MultiQC custom-content tables instead of emitting an empty "did not find native tool logs"
  placeholder.
- Removed the non-functional `rambo_external` value from the `mixed_denoise_backend` schema enum and
  the dead branch in `bin/denoise_mixed_templates.py` (it never ran an external RAMBO; it silently
  greedy-fell-back), so `--help` no longer advertises an unimplemented backend.
- Remove the unimplemented `pod5` value from the `input_type` enum (`nextflow_schema.json` + the
  `main.nf` guard). HÆMA ingests already-basecalled, already-demultiplexed MinKNOW FASTQ
  (`fastq_pass/barcodeXX`); `--help` no longer advertises a POD5/Dorado basecalling path that does
  not exist, and `--input_type pod5` is now rejected by schema validation. Supported: `fastq`,
  `pooled_fastq`.
- Set `manifest.homePage` (`nextflow.config`) and `repository-code` (`CITATION.cff`) to the public
  repository URL.

### Docs
- Consolidated the docs set for public release: removed planning/internal docs
  (`IMPLEMENTATION_PLAN.md`, `REFERENCE_PANEL_REVIEW.md`, `release_checklist.md`) and merged the
  implemented-vs-staged tracking from `DELIBERATE_LIMITATIONS.md` into a feature-status table in
  `limitations.md`; updated the docs index.
- Corrected the mixed-host controls' run label (RUN01 → RUN09) in `limitations.md`/`benchmarking.md`
  and disclosed the demonstration run's actual 50% (3/6) mixed-host recovery, including the
  *Bos taurus* no-host-signal control, to match the pipeline's own HTML report.
- Scrubbed real run/sample identifiers from the example samplesheet and docs (use `EXAMPLE_RUN`,
  `SAMPLE01`, etc.).

## [0.2.0] - 2026-06-13

### Added
- Self-contained repository: curated reference panel, taxonomy sidecar, target manifest,
  and primer file are now vendored under `assets/` so a fresh clone runs `-profile test`
  with no external data.
- Release artifacts: `LICENSE` (MIT), `CITATION.cff`, `CONTRIBUTING.md`,
  `CODE_OF_CONDUCT.md`, GitHub issue/PR templates, and a CI workflow.
- **Automated control checking & mixed-host recovery benchmark**: single-host *and* lab-prepared
  mixed-host controls are compared to their declared `expected_host_scientific_name` (semicolon-
  separated for mixtures), producing `positive_control_check.tsv` with per-host recovery, missing,
  and unexpected hosts, plus a `mixed_host_recovery_rate` in `rambo_model_summary.tsv` and the
  report. Lab-prepared mixtures can now calibrate the denoising thresholds — see `docs/benchmarking.md`.
- **Scientific-caveat safeguards**: the HTML report carries a caveats banner and a
  implemented/validated/staged feature-status section; host fractions are flagged as evidence
  (`host_fractions_benchmarked: false`) in the report and `run_manifest.json`; the manifest records
  Medaka execution mode (`disabled` / `enabled_cpu` / `enabled_gpu_profile`).
- **Reference integrity**: `assets/references/CHECKSUMS.sha256` and a `verify_reference_assets.py`
  checker (checksums + required sidecar fields), wired into CI and the pre-release script.
- **Testing & release tooling**: `Makefile`, `tests/validate_release.sh`, and Python unit tests
  (validation, positive/mixed-control logic, taxdump LCA, curated-panel taxid assignment,
  endpoint columns).

### Clarified
- Taxid-backed LCA (`--taxonomy_assignment_method taxid_lca --taxdump_dir`) works against the
  curated panel via sidecar taxid backfill (proven by `tests/test_taxid_assignment.py`); panel
  re-keying to accession seqids is **not** required, correcting an earlier overstated limitation.
- `gpu` profile (optional Medaka GPU acceleration; CPU is the default fallback).
- `taxid_source` / `taxid_map_records` provenance in the BLAST-DB step.
- Documentation: usage tutorial, output interpretation, methods draft, limitations, naming
  convention, and a release checklist under `docs/`.

### Changed
- **Schema-based parameter validation and grouped `--help`** via the pinned `nf-schema@2.2.0`
  plugin: `nextflow run . --help [group]` lists parameters by category, every run validates flags
  against `nextflow_schema.json` (mistyped params are flagged), and prints a non-default parameter
  summary. The schema is reorganised into documented groups.
- Single source of truth for reference assets: removed duplicate copies of the reference panel,
  primers, and samplesheet template from the project root and within `assets/`; all bundled data
  now lives only under `assets/`.
- Every run logs a "Feature gates" line so default on/off states are visible at a glance.
- `haema-r` image build hardened to use a reliable CRAN mirror (`cloud.r-project.org`) with a
  longer timeout, fixing intermittent `SSL connect error` failures; both custom images build and
  self-test their imports.
- Input validation accepts decimal `latitude`/`longitude` and the legacy `collection_cordinates`
  column without spurious warnings; samplesheets standardised on the full MIEM/MIMARKS header.
- All container images are pinned by immutable digest (or, for `ontresearch/medaka`, an
  immutable upstream sha-tag); nothing uses `:latest`.
- Output, log, and work directories now default to `launchDir`-relative paths (portable).
- `--input` and `--raw_data_dir` are now required with no default (fail-fast validation);
  use `-profile test` for a zero-config first run.
- Mixed-template denoising defaults tuned for real ONT depth
  (`min_cluster_size` 3→20, `min_cluster_fraction` 0.01→0.05) to avoid over-splitting.
- Per-attempt resource scaling with `resourceLimits` ceilings and OOM-aware retries.

### Fixed
- `makeblastdb` no longer fails fatally on the real panel's long (>50 char) pipe-delimited
  deflines (`-parse_seqids`/`-taxid_map` dropped; taxonomy resolved via the curated sidecar).
- UMAP/HDBSCAN now actually runs under Nextflow's non-root containers: `NUMBA_CACHE_DIR`/`HOME`
  point at a writable work dir (previously the denoiser silently fell back to greedy).
- Default Python image is the full `python:3.11` (not `-slim`) because Nextflow requires `ps`
  (procps) for task metrics; custom images install `procps` explicitly.
- Removed dead parameters (`allow_remote_blast`, `expected_positive_hosts`) and wired the
  previously-inert `cleanup` parameter to Nextflow's `cleanup` directive.
- `nextflow_schema.json` synced with all parameters.
- Multi-dimensional audit pass: `run_manifest.json` now records sha256 content hashes of the
  resolved reference assets (provenance); a runtime WARNING fires when `*_lca` is requested without
  a `--taxdump_dir` (ambiguous hits then collapse by defline genus, not taxid lineage); the
  taxonomy-skipped path emits the same column schema as the BLAST path; removed a dead `start`
  branch in `fastq_utils`; renamed the merge stat `mean_q`→`mean_read_q` (unweighted per-read mean);
  fixed the stale `homePage` and `.git`/checklist items.

## [0.1.0] - 2026-06-09

### Added
- Initial DSL2 pipeline scaffold: input validation, marker-aware trim/QC/split,
  pre-consensus denoising, consensus/dereplication, curated local-first BLAST + LCA,
  control-aware contamination flags, RAMBO-style mixed-host evidence, phyloseq/decontam
  endpoints, custom + MultiQC reports, and `test`/`production` profiles.
