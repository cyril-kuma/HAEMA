# HÆMA ONT Blood-Meal Nextflow Implementation Plan

## Blueprint Requirements

The architecture blueprint requires a reproducible DSL2 pipeline for ONT/GridION blood-meal metabarcoding from a tri-marker multiplex PCR assay: `cyt_b`, long COI, and short COI. Required stages are metadata validation, FASTQ or future POD5 ingestion, demultiplexing/barcode validation, marker-aware primer trimming, ONT read quality and length filtering, marker splitting, mixed-template-aware feature generation, local-first taxonomy, control-aware contamination assessment, structured endpoint outputs, and clear reports.

## Confirmed Inputs

- Project root: `<project-root>`
- Pipeline directory: `<repo>`
- Raw data directory: `<project-root>/data/Runs`
- Current input type: already-basecalled GridION `.fastq.gz` from MinKNOW `fastq_pass/barcodeXX` folders
- Samplesheet: `<project-root>/metadata/LUC_MOSBMA_samplesheet.csv`
- Primers: `<project-root>/metadata/primers.csv`
- Curated FASTA: `<project-root>/references/vertebrate_dna_ref_panel.fasta`
- Local NCBI nt prefix: `<path/to/blast_db>/nt/nt`
- Runtime layout: `results/`, `logs/`, and `work/` live under the project root, outside `pipeline/`

## Still Needed From The User

- Barcode kit/library prep, for example `SQK-RBK114-96`
- Exact Dorado SUP model used in MinKNOW, for manifest completeness
- Positive-control expected host identities
- Final acceptable amplicon windows if primer CSV `Size ± padding` is not sufficient
- Curated reference sidecar for the real panel: `seqid,taxid,scientific_name,rank`
- NCBI taxdump path if full ancestor-based LCA is required
- SLURM account/partition/queue and Apptainer/Singularity cache settings for HPC
- Whether coordinate values such as `8.785-1.535` should be normalized to comma-delimited latitude/longitude

## Implemented Decisions

- Docker is default for local execution.
- MinKNOW barcode folders are trusted as sample assignment; filename/folder mismatches are reported unless strict mode is enabled.
- Marker windows are read from the primer CSV unless `--marker_windows` is provided.
- Metadata, barcode, marker, run, and control state are retained in endpoint tables.
- Curated BLAST is primary; local `nt` is a fallback.
- Rows are flagged, not silently removed.
- Mixed host evidence is preserved at feature level and summarized in an optional RAMBO-style host-call table.

## Pipeline Structure

```text
pipeline/
├── main.nf
├── nextflow.config
├── nextflow_schema.json
├── bin/
├── modules/local/
├── subworkflows/local/
├── assets/
├── conf/
├── docs/
└── test_data/
```

## Results Structure

```text
results/
├── 00_demultiplexing/
├── 01_ingested/
├── 02_trimmed_filtered/
├── 03_consensus_variants/
├── 04_taxonomy/
│   ├── blast/
│   ├── curated_reference/
│   ├── assignments/
│   └── evidence/
├── 05_endpoint_files/
└── 06_reports/

logs/
└── pipeline_info/

work/
└── Nextflow task work directories
```

## Modules And Subworkflows

- `INPUT_VALIDATION`: samplesheet, primers, FASTQ inventory, barcode mismatch reporting
- `ADVANCED_DEMUX`: pooled FASTQ header-tag demux and external command-template demux wrapper
- `PREPROCESS_READS`: FASTQ merge, trim/filter/split into markers
- `DEREPLICATE_ASVS`: exact sequence ASV-like mode
- `CLUSTER_CONSENSUS`: greedy cluster/consensus mode
- `MEDAKA_POLISH`: optional consensus polishing
- `TAXONOMY_LOCAL_BLAST`: curated-only, curated-then-fallback, or external DB taxonomy
- `BUILD_CURATED_TAXID_MAP`: curated reference metadata validation and BLAST taxid map creation
- `ASSIGN_TAXONOMY`: conservative top-hit, exact-taxid LCA, or taxdump-backed LCA
- `AGGREGATE_RESULTS`: primary endpoints, QC, summaries, control flags
- `RAMBO_MIXED_MODEL`: practical mixed-host abundance/evidence layer
- `BUILD_R_OUTPUTS`: formal/fallback phyloseq and decontam outputs
- `MULTIQC_REPORT`: optional MultiQC HTML
- `BUILD_REPORT`: custom HÆMA report
- `BUILD_RUN_MANIFEST`: machine-readable run manifest

## Key Configuration Parameters

- Input: `input`, `primers`, `raw_data_dir`, `input_type`, `pooled_fastq`, `require_fastqs`
- Demux: `enable_advanced_demux`, `advanced_demux_tool`, `advanced_demux_command_template`, `demux_run_id`
- Preprocess: `min_mean_q`, `min_read_length`, `primer_search_window`, `primer_max_error_rate`, `marker_windows`
- Feature generation: `consensus_method`, `min_asv_reads`, `min_asv_fraction`, `min_cluster_identity`
- Medaka: `enable_medaka`, `medaka_model`, `medaka_container`, `medaka_extra_args`
- Taxonomy: `reference_fasta`, `curated_reference_metadata`, `taxonomy_strategy`, `fallback_blast_db`, `blast_db_mount`, `taxdump_dir`
- LCA strictness: `require_curated_taxids`, `require_taxdump_for_lca`, `taxonomy_assignment_method`
- Mixed host model: `enable_rambo_model`, `rambo_min_host_reads`, `rambo_min_host_fraction`, `rambo_include_contaminants`
- R outputs: `enable_r_outputs`, `enable_phyloseq`, `enable_decontam`, `strict_bioconductor`, `decontam_threshold`
- Reporting: `enable_multiqc`, `outdir`, `log_dir`

## Tested Commands

```bash
python3 -m py_compile pipeline/bin/*.py
python3 -m json.tool pipeline/nextflow_schema.json >/tmp/haema_schema_check.json
nextflow config pipeline -profile test
```

```bash
NXF_ANSI_LOG=false nextflow -log logs/nextflow-r-output-test.log run pipeline \
  -profile test,local \
  --skip_taxonomy false \
  --enable_r_outputs true \
  --enable_phyloseq true \
  --enable_decontam true \
  --strict_bioconductor false \
  --enable_multiqc false \
  --outdir results/test_r_outputs \
  --log_dir logs/test_r_outputs
```

```bash
NXF_ANSI_LOG=false nextflow -log logs/nextflow-advanced-demux-test.log run pipeline \
  -profile test,local \
  --input_type pooled_fastq \
  --enable_advanced_demux true \
  --pooled_fastq 'pipeline/assets/test_data/pooled/pooled.fastq' \
  --demux_run_id RUN_TEST \
  --skip_taxonomy true \
  --enable_r_outputs false \
  --enable_multiqc false \
  --outdir results/test_advanced_demux \
  --log_dir logs/test_advanced_demux
```

```bash
NXF_ANSI_LOG=false nextflow -log logs/nextflow-medaka-stub-test.log run pipeline \
  -profile test,local \
  -stub-run \
  --enable_medaka true \
  --medaka_container python:3.11 \
  --skip_taxonomy false \
  --enable_r_outputs false \
  --enable_multiqc false \
  --outdir results/test_medaka_stub \
  --log_dir logs/test_medaka_stub
```

## Remaining Validation Work

- Build the real curated reference metadata sidecar and taxdump-backed LCA test.
- Validate positive-control expected hosts in the project samplesheet.
- Benchmark RAMBO-style thresholds against known mixed-host controls.
- Pin a Medaka image/model after checking the exact installed model names.
- Add an HPC profile once SLURM and Apptainer details are known.
