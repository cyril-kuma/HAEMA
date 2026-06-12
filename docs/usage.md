# Usage tutorial

A step-by-step walkthrough from install to a real run. For parameter details see
[`parameters.md`](parameters.md) and [`nextflow_schema.json`](../nextflow_schema.json); for output
interpretation see [`output.md`](output.md).

## 1. Install

```bash
curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/
nextflow info          # confirm Nextflow >= 24.10.0 and Java 17-24
docker --version       # or: singularity --version
git clone <repository-url> haema && cd haema
```

## 2. First run (bundled demo)

```bash
nextflow run . -profile test,docker --skip_taxonomy false --outdir results/test
```

Expect `Execution complete` in ~1 minute and a populated `results/test/05_endpoint_files/`.
This proves Nextflow, your container engine, and the pipeline logic all work.

## 3. Prepare your inputs

**Samplesheet** — copy and edit the template:

```bash
cp assets/samplesheets/example_samplesheet.csv my_samplesheet.csv
```

- `run_id` must equal the top-level folder under your data root.
- `barcode_id` (`barcodeNN`) must match the MinKNOW `fastq_pass/barcodeNN` folder.
- `sample_type` is one of `sample`, `positive_control`, `negative_control`. **Include at least one
  negative control** — background subtraction depends on it.
- For each `positive_control`, set `expected_host_scientific_name` (e.g. `Bos taurus`) so the
  pipeline can automatically check that the control recovered its known host
  (`06_reports/positive_control_check.tsv`; pass/fail surfaced in the report).
- Fill ecological metadata as completely as you can (coordinates, zone, species, batches).
  `-profile production` requires the full MIEM/MIMARKS field set.

**Data layout** the validator expects:

```text
Runs/
└── <run_id>/<minknow_subfolder>/fastq_pass/<barcodeNN>/*.fastq.gz
```

**Primers / references** — the bundled `assets/primers.csv` and `assets/references/` work out of
the box; override with `--primers` / `--reference_fasta` for your own.

## 4. Enable real UMAP/HDBSCAN denoising (recommended for real data)

The default `python:3.11` image cannot run UMAP, so denoising falls back to greedy clustering.
Build the scientific image once and pass it in:

```bash
docker build -t haema-python:0.3.0 -f containers/haema-python/Dockerfile .
```

## 5. Run on your data

```bash
nextflow run . -profile local \
  --input        "$PWD/my_samplesheet.csv" \
  --raw_data_dir "$PWD/Runs" \
  --taxonomy_strategy curated_only \
  --python_container haema-python:0.3.0 \
  --enable_rambo_model true \
  --outdir results/myrun --log_dir logs/myrun -resume
```

Use **absolute paths**. Add `--fallback_blast_db /path/to/nt/nt --blast_db_mount /path/to/blast_db`
to enable the `nt` fallback for hosts outside the curated panel.

## 6. Inspect results

```bash
column -t -s$'\t' results/myrun/05_endpoint_files/host_call_table.tsv | less -S
open results/myrun/06_reports/bloodmeal_pipeline_report.html      # macOS; use xdg-open on Linux
```

Look for `mixed_status = mixed_host` rows (a mosquito that fed on >1 host), check the negative
control has no host calls, and confirm the positive control matches its known host.

## 7. Tune for your depth

- Shallow samples / want more sensitivity to minor hosts → lower
  `--mixed_denoise_min_cluster_fraction` (e.g. 0.02) and `--mixed_denoise_min_cluster_size`.
- Over-splitting / too many clusters and slow taxonomy → raise both.
- No host expected outside the panel → keep `--taxonomy_strategy curated_only` (faster, frozen).

## 8. Validate & test (developers / before a release)

```bash
make lint                      # python compile, schema, config (all profiles), refs, unit tests
make test                      # bundled test-profile run (needs Docker)
bash tests/validate_release.sh --run   # full pre-release validation incl. a test run
python3 bin/verify_reference_assets.py --assets-dir assets/references   # checksum + sidecar check
```

## 9. HPC

```bash
export NXF_SINGULARITY_CACHEDIR=/shared/singularity
nextflow run . -profile slurm \
  --input ... --raw_data_dir ... --outdir ... -resume
```

Set your partition/account via `process.queue` / `process.clusterOptions` (a `-c custom.config`)
and adjust `process.resourceLimits` to your node sizes.
