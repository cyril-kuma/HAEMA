# Parameters

Every parameter is documented with type and default in
[`nextflow_schema.json`](../nextflow_schema.json). Override any of them on the command line
(`--param value`) or with a `-params-file params.yml`. This page summarises the ones you are most
likely to set.

## Inputs (required / core)
| Parameter | Default | Description |
|---|---|---|
| `--input` | — (required) | Samplesheet CSV. |
| `--raw_data_dir` | — (required, FASTQ mode) | Root holding `<run_id>/.../fastq_pass/<barcode>`. |
| `--primers` | `assets/primers.csv` | Primer CSV (`Gene,Forward_Primer,Reverse_Primer,Size`). |
| `--reference_fasta` | bundled panel | Curated vertebrate reference FASTA. |
| `--curated_reference_metadata` | bundled sidecar | `seqid,taxid,scientific_name,rank,...` taxonomy sidecar. |
| `--outdir` / `--log_dir` | `<launchDir>/results` / `logs` | Output / log locations. |

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
| `--enable_medaka` | false | Medaka polishing (needs a Medaka container + model). |
| `--medaka_model` | `r1041_e82_400bps_sup_v4.3.0` | Must exist in the Medaka image (checked at runtime). |

## Taxonomy
| Parameter | Default | Description |
|---|---|---|
| `--skip_taxonomy` | false | Emit unassigned rows instead of running BLAST. |
| `--taxonomy_strategy` | `curated_then_fallback` | `curated_only` / `curated_then_fallback` / `nt_only`. |
| `--fallback_blast_db` | "" | External BLAST db prefix (e.g. `/path/nt/nt`) for unresolved features. |
| `--blast_db_mount` | "" | Host dir mounted read-only so the db prefix is visible in containers. |
| `--min_blast_identity` / `--min_blast_coverage` | 97 / 80 | Hit acceptance thresholds. |
| `--taxonomy_assignment_method` | `conservative_lca` | `conservative_lca` / `taxid_lca` / `top_hit`. |

## Contamination, host model & R outputs
| Parameter | Default | Description |
|---|---|---|
| `--enable_rambo_model` | true | RAMBO-style mixed-host evidence & host calls. |
| `--enable_r_outputs` | true | phyloseq/decontam endpoints (needs `--enable_rambo_model`). |
| `--strict_bioconductor` | false | Fail if formal phyloseq/decontam are missing (needs `haema-r`). |
| `--decontam_threshold` | 0.5 | decontam / negative-control prevalence threshold. |

## Containers (pinned defaults; see docs/CONTAINER_STRATEGY.md)
`--python_container`, `--blast_container`, `--r_container`, `--medaka_container`,
`--multiqc_container`, `--advanced_demux_container`.

## Execution
`--cleanup` (purge work dirs on success; disables `-resume`), plus the profiles
`test`, `local`, `docker`, `singularity`, `apptainer`, `slurm`, `gpu`, `production`.
