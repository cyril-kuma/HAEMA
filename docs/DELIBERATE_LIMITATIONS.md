# Deliberate Limitations And Staged Enhancements

Version `0.1.0` is now a production-gated staged pipeline rather than only a scaffold. The workflow structure implements the blueprint path for MinKNOW FASTQ input, pre-consensus mixed-template separation, cluster-level consensus/Medaka polishing, local-first BLAST/LCA, formal or fallback decontamination, phyloseq-ready outputs, and a machine-readable run manifest. Remaining limitations are mostly external-dependency, missing-resource, and scientific-validation issues.

## MinKNOW FASTQ Input

The primary input remains MinKNOW `fastq_pass/barcodeXX` folders. These reads were already generated on a GridION and basecalled with a Dorado SUP model, so they are comparable to re-basecalling from `.pod5` with the same model and settings.

The missing `.pod5` files would mainly preserve future optionality: newer Dorado models, raw-signal auditing, alternative raw demultiplexing, and recovery of reads that never reached `fastq_pass`.

## Advanced Demultiplexing

Implemented now:

- `--input_type pooled_fastq`
- `--advanced_demux_tool header_tag`
- external demux wrapper support through `--advanced_demux_command_template`
- production default `--demux_strategy pre_demultiplexed_minknow_trusted_folder`

Still external:

- Barbell and Deepbinner are not bundled or validated locally. They can be run through the command-template wrapper if the user supplies the tool, container, barcode model, and a command that writes MinKNOW-like `fastq_pass/barcodeXX` output under `{output}`.
- POD5/Dorado basecalling and raw-signal demultiplexing remain future work.

## Mixed Blood Meals

Implemented now:

- pre-consensus `DENOISE_MIXED_TEMPLATES` on marker FASTQs
- UMAP/HDBSCAN k-mer clustering with deterministic greedy fallback
- cluster FASTQs, cluster membership tables, and denoising summaries
- cluster-aware consensus, dereplication, BLAST, Medaka, and endpoint provenance
- exact dereplication and greedy cluster consensus modes
- feature-level mixed-template preservation
- `RAMBO_MIXED_MODEL`, a practical RAMBO-style abundance/evidence layer that keeps multiple host signals per sample/marker when they pass configurable read and fraction thresholds
- explicit `mixed_host_evidence.tsv` and `host_call_table.tsv`
- automated single-host positive-control checking against `expected_host_scientific_name` (`positive_control_check.tsv`)
- explicit `host_fractions_benchmarked: false` flags in the report and run manifest so fractions are not read as validated estimates

Still needs scientific validation:

- true RAMBO-style probabilistic mixed-template modelling against mixed-host positive controls
- calibration of `mixed_denoise_min_cluster_size` / `min_cluster_fraction` against known mixtures (a mixed-host positive-control set); the checking framework is implemented but needs the wet-lab controls
- benchmarked thresholds for low-abundance host calls

## Medaka

Implemented now:

- optional `MEDAKA_POLISH` module
- configurable `--medaka_model`, `--medaka_container`, and `--medaka_extra_args`
- clear failure if `medaka_consensus` is unavailable when enabled
- `MEDAKA_MODEL_PREFLIGHT` checks that the requested model is available when Medaka is enabled
- production profile requires `r1041_e82_400bps_sup_v4.3.0` unless overridden
- production pins the official `ontresearch/medaka` image by immutable sha-tag (never `:latest`) and verifies the model at runtime, so no custom Medaka image is required
- the optional `haema-medaka` Dockerfile (for build-time model assertion / air-gapped freezing) defaults to the same pinned sha-tag and is overridable through `MEDAKA_BASE_IMAGE`

Still needs validation:

- mixed-template positive controls to ensure polishing does not merge host evidence
- confirmation that the pinned Medaka build lists `r1041_e82_400bps_sup_v4.3.0` (asserted automatically by `MEDAKA_MODEL_PREFLIGHT` on the first production run)

## Taxonomy And LCA

Implemented now:

- curated-panel BLAST
- optional `nt` or other local fallback BLAST
- curated reference metadata sidecar with `seqid,taxid,scientific_name,rank,source_accession,sequence_md5`
- accession-driven Ghana/West African reference target manifest
- `build_reference_panel.py` to exclude placeholder or checksum-incomplete targets
- `makeblastdb -taxid_map` when taxids are present
- exact single-taxid LCA assignment without taxdump
- true taxdump-backed LCA when `--taxdump_dir` contains `nodes.dmp` and `names.dmp`

Still needed for publication-scale curated-panel expansion:

- marker-specific extracted reference panels
- database versioning and checksums
- accessions and checksums for target-manifest rows marked `needs_accession`
- an NCBI taxdump snapshot in `references/taxdump` or another configured location

## Contamination And R Ecology Outputs

Implemented now:

- negative-control prevalence/background threshold reporting
- `decontam_results.tsv`
- `contamination_model_summary.tsv`
- `bloodmeal_phyloseq.rds`
- `bloodmeal_ecology_data.rds`
- `bloodmeal_ecology_data_decontaminated.rds`
- `host_calls_decontaminated.tsv`
- `asv_count_table_decontaminated.tsv`
- `R_BIOCONDUCTOR_PREFLIGHT` for strict package checks

If `phyloseq` and `decontam` are installed in the R runtime, the R stage uses formal package objects/results. If they are unavailable and `--strict_bioconductor false`, the pipeline writes documented fallback objects and tables. Use `--strict_bioconductor true` to require formal Bioconductor packages.

## Reporting

Implemented now:

- custom HÆMA HTML report
- optional MultiQC HTML report
- machine-readable `run_manifest.json`
- production preflight report and summary under `pipeline_info/production_preflight`

The custom HÆMA report remains the primary interpretation report because many steps are custom scripts rather than native MultiQC-supported tools.
