# Parameters

Every parameter is documented with type and default in
[`nextflow_schema.json`](../nextflow_schema.json). Override any of them on the command line
(`--param value`) or with a `-params-file params.yml`. This page summarises the ones you are most
likely to set.

> **Discover & validate from the CLI:** `nextflow run . --help` prints all parameters grouped by
> category; `nextflow run . --help taxonomy` (or any group) drills in. Parameters are validated
> against the schema on every run (via the `nf-schema` plugin), so a mistyped flag like `--inputt`
> is flagged with a clear warning instead of being silently ignored.

## Inputs (required / core)
| Parameter | Default | Description |
|---|---|---|
| `--input` | — (required) | Samplesheet CSV. |
| `--raw_data_dir` | — (required, FASTQ mode) | Root holding `<run_id>/.../fastq_pass/<barcode>`. |
| `--primers` | `assets/primers.csv` | Primer CSV (`Gene,Forward_Primer,Reverse_Primer,Size`). |
| `--reference_fasta` | bundled panel | Curated vertebrate reference FASTA. |
| `--curated_reference_metadata` | bundled sidecar | `seqid,taxid,scientific_name,rank,...` taxonomy sidecar. |
| `--outdir` / `--log_dir` | `<launchDir>/results` / `logs` | Output / log locations. |

## Why Some Defaults Are `false`
The default profile is intentionally easy to launch on a laptop or workstation. `false` usually means
"optional gate disabled until you ask for it", not "missing configuration": advanced pooled-FASTQ
demultiplexing, Medaka polishing, strict metadata/barcode checks, required curated taxids, required
NCBI taxdump, strict Bioconductor enforcement, and work-directory cleanup all need extra inputs,
containers, or operational choices. `-profile production` turns on the strict scientific checks and
full production feature set. The tiny `test` profile also disables taxonomy/R/MultiQC by default for
speed; pass `--skip_taxonomy false` when you want the demo to exercise curated BLAST.

## QC & marker splitting
| Parameter | Default | Description |
|---|---|---|
| `--min_mean_q` | 20 | Minimum per-read mean Phred. |
| `--min_read_length` | 100 | Minimum read length before splitting. |
| `--marker_windows` | derived from primer `Size` ± `marker_window_padding` | `marker:min:max` length windows. |
| `--primer_max_error_rate` | 0.15 | Allowed primer mismatch fraction. |

## Mixed-template denoising
| Parameter | Default | Description |
|---|---|---|
| `--enable_mixed_denoising` | true | Run pre-consensus denoising. |
| `--mixed_denoise_backend` | `umap_hdbscan` | `umap_hdbscan` / `greedy` / `none` (falls back to greedy when UMAP is unavailable). |
| `--mixed_denoise_min_cluster_size` | 20 | Min reads per retained host cluster (lower = more sensitive). |
| `--mixed_denoise_min_cluster_fraction` | 0.05 | Min read fraction per retained cluster. |
| `--mixed_denoise_min_reads_for_umap` | 50 | Below this, use the greedy fallback. |

## Consensus & polishing
| Parameter | Default | Description |
|---|---|---|
| `--consensus_method` | `cluster` | `cluster` (greedy consensus) or `dereplicate` (exact). |
| `--enable_medaka` | **true** | Medaka consensus polishing before taxonomy (on by default; needs the ONT Medaka image + a model present in it; the `test` profile disables it for speed). |
| `--medaka_model` | `r1041_e82_400bps_sup_v4.3.0` | Must exist in the Medaka image (checked at runtime). |

## Taxonomy
| Parameter | Default | Description |
|---|---|---|
| `--skip_taxonomy` | false | Emit unassigned rows instead of running BLAST. |
| `--reference_mode` | `curated_panel` | Canonical reference selector: `curated_panel` (A) / `broad_blast` (B) / `remote_fallback` (C) / `bold_aware` (D). See docs/methods.md. |
| `--taxonomy_strategy` | `curated_then_fallback` | Legacy selector, still honoured and mapped onto `--reference_mode` when the latter is default. `curated_only` / `curated_then_fallback` / `nt_only`. |
| `--blast_db` | "" | Broad local BLAST db prefix for Mode B (`broad_blast`). |
| `--enable_curated_panel_check` | true | In Mode B, query the curated panel first and use `--blast_db` only for unresolved features. |
| `--fallback_blast_db` | "" | Local external BLAST db prefix for unresolved features (legacy `curated_then_fallback`). |
| `--enable_ncbi_remote_fallback` | true | Allow the NCBI remote fallback in Mode C (`remote_fallback`). |
| `--remote_blast_db` | `core_nt` | NCBI-hosted db queried with `blastn -remote` in Mode C. Not reproducible by default. |
| `--bold_fasta` / `--bold_taxonomy` | "" | BOLD-derived COI FASTA (+ optional taxonomy sidecar) for Mode D (`bold_aware`). |
| `--bold_mode` | `local_fasta` | `local_fasta` (reproducible; the only implemented path) or `api_query` (documented future work). |
| `--coi_species_identity_threshold` / `--cytb_species_identity_threshold` | 98 / 97 | Marker species-level confidence guard (downgrades `high`→`medium` below threshold; never loosens `--min_blast_identity`). |
| `--blast_db_mount` | "" | Host dir mounted read-only so the db prefix is visible in containers. |
| `--min_blast_identity` / `--min_blast_coverage` | 97 / 80 | Global hit acceptance thresholds. |
| `--taxonomy_assignment_method` | `conservative_lca` | `conservative_lca` / `taxid_lca` / `top_hit`. |
| `--marker_numt_risk` | "" | JSON overriding per-marker NUMT risk labels (else defaults: co1_short=moderate, co1_long=low, cyt_b=low). |
| `--enable_marker_concordance` | true | Compute multi-marker concordance per specimen. |
| `--enable_host_ecology_comparisons` | true | Compute exploratory Fisher/Holm host-use comparisons between strata. |

## Contamination, host model & R outputs
| Parameter | Default | Description |
|---|---|---|
| `--enable_rambo_model` | true | RAMBO-style mixed-host evidence & host calls. |
| `--rambo_min_host_reads` | 3 | Min supporting reads per host per sample/marker for a host call. |
| `--rambo_min_host_fraction` | 0.02 | Min host read fraction per sample/marker. **In-silico calibrated** (docs/denoising_calibration): 0% false positives at 100% detection ≥5% minor host; lower toward 0.01 to chase <5% minors (more false positives). |
| `--non_host_genera` | "" | Genera excluded as non-host (vector self-hits); empty = built-in vector default. |
| `--enable_r_outputs` | true | phyloseq/decontam endpoints (needs `--enable_rambo_model`). |
| `--strict_bioconductor` | false | Fail if formal phyloseq/decontam are missing (needs `haema-r`). |
| `--decontam_threshold` | 0.5 | decontam / negative-control prevalence threshold. |

## Figures
| Parameter | Default | Description |
|---|---|---|
| `--enable_figures` | true | Render publication figures to `07_figures/` (needs `haema-figures`; off in `-profile test`). See [figures.md](figures.md). |
| `--figure_formats` | `pdf,svg,png` | Comma-separated output formats (any subset of `pdf`, `svg`, `png`). |
| `--enable_publication_figures` | true | Render the curated Objective 1 figure suite to `--publication_figures_dir` (off in `-profile test`). |
| `--publication_figures_dir` | `results/figures` | Output directory for the curated suite (`main/`, `supplementary/`, `figure_data/`). |
| `--figure_bioclim_column` | `bioclimatic_zone` | Metadata column used for bioclimatic-zone stratification in the curated suite. |

## Containers (pinned defaults; see docs/CONTAINER_STRATEGY.md)
`--python_container`, `--blast_container`, `--r_container`, `--medaka_container`,
`--multiqc_container`, `--figures_container`, `--advanced_demux_container`.

## Execution
`--cleanup` (purge work dirs on success; disables `-resume`), plus the profiles
`test`, `local`, `docker`, `singularity`, `apptainer`, `slurm`, `gpu`, `production`.
