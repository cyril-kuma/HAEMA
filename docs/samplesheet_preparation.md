# Preparing a samplesheet

HÆMA is driven by a single CSV samplesheet (`--input`) plus a raw-data directory
(`--raw_data_dir`). The pipeline discovers `*.fastq.gz` under
`raw_data_dir/<run_id>/<minknow_run_folder>/fastq_pass/<barcode_id>/` and matches each
barcode to a samplesheet row by `(run_id, barcode_id)`.

## Required and recommended columns

The schema (validated by `bin/validate_inputs.py`) expects, per row:

| Column | Required | Purpose |
|--------|----------|---------|
| `run_id` | yes | Top-level run folder name under `--raw_data_dir` |
| `minknow_run_folder` | yes | MinKNOW output folder holding `fastq_pass/` |
| `barcode_id` | yes | e.g. `barcode07`; the `fastq_pass/` subfolder |
| `sample_id` | yes | Human-readable sample label |
| `specimen_id` | yes | **Stable specimen key for cross-objective integration — keep unique** |
| `sample_type` | yes | `sample` / `positive_control` / `negative_control` |
| `expected_host_scientific_name` | controls | Expected host(s); `;`-separate a mixed control (e.g. `Homo sapiens;Bos taurus`) |
| `sibling_species`, `collection_region`, `bioclimatic_zone`, `collection_date` | recommended | Stratification (zone / species / season) for ecological indices and comparisons |
| `barcode_kit`, `flowcell`, `basecalling_model` | recommended | ONT run metadata / provenance |

Controls are **excluded from ecological indices** and used for contamination flagging;
`expected_host_scientific_name` drives the RAMBO positive-control check.

## Building a subset (test/validation) samplesheet

For a fast end-to-end validation, subset to the rows that exercise the science:

- **all** `positive_control` and `negative_control` rows (these are your assertions),
  including any mixed-fed positive controls (`A;B` expected hosts);
- a few field `sample` rows spanning each `bioclimatic_zone` / `sibling_species`.

Only the listed barcodes are ingested, so a subset drawn from large runs stays cheap.
Always confirm the raw data exists for each chosen row before running, e.g.:

```bash
# keep only rows whose fastq_pass/<barcode_id> directory has reads
awk -F, 'NR==1{print;next}
  {p="'"$RAW"'/"$2"/*/fastq_pass/"$4; cmd="ls "p"/*.fastq.gz 2>/dev/null | head -1";
   cmd|getline f; close(cmd); if(f!="") print}' samplesheet.csv > subset.csv
```

This repository's real run used `metadata/subset_validation_samplesheet.csv` (27 rows:
13 positive incl. 3 mixed-fed, 8 negative, 6 field across 3 zones) for validation and
`metadata/full_run_samplesheet.csv` (all rows with data on disk) for the full run. Both
were generated from `metadata/LUC_MOSBMA_RUN_samplesheet.csv` by keeping only rows whose
`(run_id, barcode_id)` resolve to real reads.

## Building the full samplesheet

The full samplesheet is the master samplesheet filtered to rows whose raw data is present.
Rows without data on disk are dropped (and reported) rather than failing the run.

## Reproducibility notes

- Paths are **not** hard-coded in the pipeline; pass `--input` and `--raw_data_dir` at
  runtime. Keep the samplesheet and raw data on the same machine or a shared mount.
- Another user preparing their own experiment only needs to reproduce the column schema
  above and point `--raw_data_dir` at their MinKNOW output tree. See
  `assets/samplesheets/example_samplesheet.csv` and
  `assets/samplesheets/mixed_host_control_example.csv` for templates.
