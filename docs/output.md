# Output files & how to interpret them

All outputs are written under `--outdir` (default `<launchDir>/results`). Heavy intermediates stay
as work-directory symlinks; only the files below are copied (`publishDir mode: 'copy'`).

```text
results/
├── 01_ingested/            <run>/<barcode>/<sample_uid>.merged.fastq.gz + qc/*.merge_stats.tsv
├── 02_trimmed_filtered/    <run>/<barcode>/<sample_uid>.<marker>.fastq.gz + qc/ summaries & decisions
├── 03_consensus_variants/  mixed_denoising/ (cluster FASTQs, membership, summaries),
│                           consensus/ASV FASTAs + counts, medaka/ (if enabled)
├── 04_taxonomy/            blast_db/, curated_reference/, raw_blast/, assignments/, evidence/
├── 05_endpoint_files/      core downstream tables + .rds objects + run_manifest.json
└── 06_reports/             bloodmeal_pipeline_report.html, multiqc_report.html, decontam_*, qc_*
pipeline_info/              execution_timeline.html, execution_report.html, execution_trace.tsv, pipeline_dag.html
```

## The files you will actually use (`05_endpoint_files/`)

### `bloodmeal_master_endpoint.tsv` — the primary deliverable
One row per **feature × sample**. Every sample, including controls, is preserved — rows are
**flagged, never deleted**. Useful columns:

- `sample_uid`, `run_id`, `sample_id`, `barcode_id`, `control_status` — identity & control type
- `marker`, `cluster_id`, `asv_id`, `count`, `fraction` — the feature and its abundance
- `host_assignment`, `taxon_rank`, `assignment_status`, `confidence`, `pident`, `coverage` — taxonomy
- `blast_source` (`curated_panel` / `nt`) and `lca_taxid` / `top_staxids` — provenance
- `contamination_flag`, `contamination_reason`, `present_in_negative_control` — QC flags

### `host_call_table.tsv` — host calls per sample/marker
The RAMBO-style abundance/evidence layer. One row per retained host per sample/marker:

- `host_assignment`, `host_reads`, `host_fraction`, `n_supporting_features`
- `mixed_status` — **`single_host`**, **`mixed_host`**, or `no_host_signal`
- `best_confidence`, `best_assignment_status`, `control_status`

> **Reading a mixed feed:** two `mixed_host` rows for the same sample/marker with different
> `host_assignment` (e.g. `Homo sapiens` + `Capra hircus`) means two hosts were detected. Cross-marker
> agreement (same pair in `cyt_b` and `co1_short`) is the strongest evidence.

### Other endpoint files
- `asv_count_table.tsv` — feature × sample count matrix.
- `sample_level_summary.tsv` / `marker_level_summary.tsv` — per-sample / per-marker rollups.
- `contamination_flags.tsv` — features flagged against negative-control background.
- `bloodmeal_phyloseq.rds`, `bloodmeal_ecology_data.rds`, `*_decontaminated.rds` — R objects
  (formal phyloseq/decontam with `haema-r`, documented fallback objects otherwise).
- `run_manifest.json` — parameters, container images, output paths, workflow session — your
  provenance record for the methods section.

### `06_reports/positive_control_check.tsv` — control / recovery check
One row per control (any sample with a declared `expected_host_scientific_name`), comparing the
**observed** host call(s) to the declared composition. Handles both single-host controls and
**lab-prepared mixed-host controls** (semicolon-separated expected hosts):

- `control_kind` — `single_host_control` or `mixed_host_control`
- `expected_hosts` / `observed_hosts`, `n_expected` / `n_recovered`
- `missing_hosts` (declared but not recovered → false negatives) and `unexpected_hosts`
  (recovered but not declared → false positives)
- `status` — single: `pass` / `pass_genus` / `fail_unexpected_host` / `fail_no_host_signal`;
  mixed: `mixed_pass_all` / `mixed_pass_with_extra` / `mixed_partial` / `mixed_fail`;
  or `indeterminate_no_expected_host`.

This is a **recovery/integrity** check, **not** a benchmark of host *fractions*.
`rambo_model_summary.tsv` rolls it up (`single_host_controls_*`, `mixed_host_controls_*`,
`mixed_host_recovery_rate`) and the HTML report shows a banner. For lab-prepared mixtures this is a
**sensitivity benchmark** — see [`benchmarking.md`](benchmarking.md) for how to calibrate the
denoising thresholds with it.

## QC you should check (`06_reports/` & earlier QC tables)
- `02_trimmed_filtered/qc/*.trim_filter_split_summary.tsv` — reads passing Q/length and assigned to
  each marker (by primer vs by length). A low `written_reads` / high `no_marker` means primers or
  length windows don't match your amplicons.
- `03_consensus_variants/mixed_denoising/qc/*.mixed_denoise_summary.tsv` — `backend_used`
  (`umap_hdbscan` = real method, `greedy_identity` = fallback), `n_clusters`, `fallback_reason`.
- `06_reports/qc_background_thresholds.tsv` — negative-control background limits (with `haema-r`).
- `bloodmeal_pipeline_report.html` — human-readable run summary (primary interpretation report).

## File-naming convention
Output files are keyed by a **compound sample UID**: `sample_uid = <run_id>__<barcode_id>__<sample_id>`,
using a **double underscore** (`__`) between the three fields because each field may itself contain
single underscores (e.g. `EXAMPLE_RUN_2025`). Splitting a name on `__` unambiguously
recovers run, barcode, and sample. The marker and step are then appended with `.` (e.g.
`EXAMPLE_RUN_2025__barcode47__SAMPLE01.cyt_b.cluster001.consensus.fasta`). This is a
deliberate, consistent scheme, not a typo.
