# Benchmarking mixed-host recovery with lab-prepared controls

Lab-prepared mixed blood meals (samples assembled from **known** host DNA in the wet lab) are the
gold-standard controls for calibrating this pipeline's mixed-template denoising. In the bundled
Ghana dataset these are the `mf`-noted RUN01 samples (barcode34, 36, 47, 56) — they are constructed
mixtures, not field unknowns, so their true host composition is known to the people who made them.

HÆMA can now use them as an automated **recovery benchmark**: did the pipeline recover every
declared host, and did it invent any it shouldn't have?

## 1. Declare the known composition

For each lab-prepared control, set `expected_host_scientific_name` in the samplesheet to the
semicolon-separated list of the hosts that were mixed in (use accepted scientific names):

```csv
run_id,barcode_id,sample_id,sample_type,control_type,expected_host_scientific_name,notes
LUC_MOSBMA_RUN01_20250721,barcode47,NS098,positive_control,mixed_host_positive_control,Homo sapiens;Capra hircus,lab-prepared mixture
```

> **Fill these from your lab records — do not copy the pipeline's own calls** (that would make the
> benchmark circular). A single-host control uses one name; a mixture uses `;`-separated names.
> Marking `sample_type=positive_control` (or leaving `sample` with a declared expected host) both work.

A ready-to-edit example is in
[`assets/samplesheets/mixed_host_control_example.csv`](../assets/samplesheets/mixed_host_control_example.csv).

## 2. Run and read the recovery report

Run as usual (with real UMAP via `--python_container haema-python:0.3.0`). The host-call model then
writes `06_reports/positive_control_check.tsv`:

| column | meaning |
|---|---|
| `control_kind` | `single_host_control` or `mixed_host_control` |
| `expected_hosts` / `observed_hosts` | declared vs recovered |
| `n_expected` / `n_recovered` | how many declared hosts were recovered |
| `missing_hosts` | declared but **not** recovered (false negatives → sensitivity loss) |
| `unexpected_hosts` | recovered but **not** declared (false positives / contamination) |
| `status` | `mixed_pass_all`, `mixed_pass_with_extra`, `mixed_partial`, `mixed_fail` (or single-host `pass`/`fail`) |

`06_reports/rambo_model_summary.tsv` rolls these up: `mixed_expected_hosts`,
`mixed_expected_hosts_recovered`, and `mixed_host_recovery_rate` (sensitivity). The HTML report shows
the same as a banner.

## 3. Calibrate the denoising thresholds

Sweep the thresholds and pick the setting that maximises recovery without introducing unexpected
hosts:

```bash
for cs in 5 10 20 30; do for cf in 0.02 0.05 0.10; do
  nextflow run . -profile local --input controls.csv --raw_data_dir Runs \
    --python_container haema-python:0.3.0 \
    --mixed_denoise_min_cluster_size $cs --mixed_denoise_min_cluster_fraction $cf \
    --outdir results/sweep_cs${cs}_cf${cf} --log_dir logs/sweep_cs${cs}_cf${cf}
  echo "cs=$cs cf=$cf:"; grep mixed_host_recovery_rate results/sweep_cs${cs}_cf${cf}/06_reports/rambo_model_summary.tsv
done; done
```

Choose the most conservative thresholds (largest `min_cluster_size`/`min_cluster_fraction`) that
still recover all declared hosts with no/low `unexpected_hosts`. Record the chosen values and the
recovery rate in your methods.

## 4. Honest scope

- This benchmarks **detection/recovery (presence-absence)** of declared hosts. It does **not** yet
  validate quantitative host *fractions* — that needs controls prepared at **known mixing ratios**.
  Until then `host_fractions_benchmarked` stays `false` in the report and manifest.
- Recovery depends on sequencing depth; benchmark at depths comparable to your real samples.
- Genus-level matches count as recovered (`pass_genus`) when the panel lacks the exact species.
