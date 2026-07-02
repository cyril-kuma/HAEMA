# In-silico denoising calibration results

The "bridge" for the absent wet-lab known-ratio controls (see
[`../denoising_redesign_plan.md`](../denoising_redesign_plan.md) §6). Synthetic known-ratio
mixed-host controls were built by pooling reads from *pure* single-host sources at defined ratios,
then classified with reference-guided read classification (ONT-adapted Logue 2016
classify-then-count), to measure the detection floor and calibrate the minor-host fraction
threshold **without any wet-lab work**.

## How to regenerate

```bash
# 1. Build synthetic mixtures (here: simulated from the curated panel + primers; --pool accepts
#    REAL pure single-host FASTQs once the full run yields cattle/goat/sheep-only samples).
python bin/make_insilico_mixtures.py --simulate \
  --reference assets/references/vertebrate_dna_ref_panel.fasta --primers assets/primers.csv \
  --hosts Homo_sapiens,Bos_taurus,Capra_hircus,Gallus_gallus \
  --reads-per-pool 3000 --depth 2000 --ratios 99:1,95:5,90:10,75:25,50:50 --replicates 3 \
  --out-dir mix --seed 42

# 2. minimap2 (ONT preset), one pass over all mixtures (read ids prefixed <mixture_id>|)
minimap2 -x map-ont --secondary=no assets/references/vertebrate_dna_ref_panel.fasta combined.fastq > combined.paf

# 3. Aggregate to detection/threshold curves
python bin/benchmark_mixed_detection.py --paf combined.paf --manifest mix/mixture_manifest.tsv \
  --output-prefix benchmark
```

## Experiment

180 synthetic mixtures = 4 hosts (Human, Cattle, Goat, Chicken — a bird included to probe
cross-taxon behaviour) × every ordered pair × 5 ratios (99:1 … 50:50) × 3 replicates, 2000 reads
each, classified against the curated panel with minimap2 `map-ont`. Reads simulated with an ONT-like
error model (~3.3 % sub+indel).

## Headline result — calibrated threshold

| `min_host_fraction` | Minor detection @ true 5 % | @ 10 % | @ 25 % | False-positive rate (any rung) |
|-----------------------|---------------------------|--------|--------|-------------------------------|
| 0.005 | 100 % | 100 % | 100 % | **25–28 %** (unacceptable) |
| 0.01 *(RAMBO current)* | 100 % | 100 % | 100 % | **11–17 %** (too permissive) |
| **0.02 (recommended)** | **100 %** | **100 %** | **100 %** | **0 %** |
| 0.05 | 50 % (misses half) | 100 % | 100 % | 0 % |
| 0.10 *(Logue 2016)* | 0 % | 56 % | 100 % | 0 % |

**Conclusions:**
- **Limit of detection ≈ 5 % minor host** at this depth/error (100 % detection for ≥ 5 %; 1 % is
  below the reliable floor).
- **`min_host_fraction = 0.02` is the operating point**: 0 % false positives, 100 % detection ≥ 5 %.
  This is more conservative than the current 1 % (which sits in the noise floor: 11–17 % false
  positives) and more sensitive than Logue's 10 % (which misses half of true-5 % and true-10 %
  minors).
- **Observed minor fractions track the true ratios** (e.g. true 25 %→~0.21, 10 %→~0.09, 5 %→~0.04),
  monotonic with a mild negative bias → fractions are qualitatively informative but **must not** be
  read as exact proportions (`host_fractions_benchmarked` stays `false`).

Full data: [`benchmark_sweep.tsv`](benchmark_sweep.tsv) (detection + false-positive by rung ×
threshold), [`benchmark_per_mixture.tsv`](benchmark_per_mixture.tsv) (observed-vs-true per mixture),
[`mixture_manifest.tsv`](mixture_manifest.tsv) (truth).

## Relationship to the pipeline, the samplesheet, and mixed-feeding controls

- **`make_insilico_mixtures.py` is a standalone developer/validation tool**, not part of the
  Nextflow workflow. It is run manually (see "How to regenerate"), **not** via `nextflow run` and
  **not** driven by the samplesheet.
- **It does not change the samplesheet format.** The production samplesheet schema
  (`docs/samplesheet_preparation.md`) is unchanged.
- **It does not remove the need for mixed-feeding controls** — it changes what they are needed *for*:
  - The tool is a computational substitute for known-*ratio* controls **for the single purpose of
    picking the detection threshold** (which previously blocked us).
  - **Real mixed-feeding positive controls are still required and still run through the pipeline**
    (`sample_type = positive_control`, `expected_host = A;B`): they validate end-to-end
    detection/recovery on *real* data via `positive_control_check.tsv` — which simulation cannot do.
  - **Quantitative** validation (read fraction ↔ blood proportion) still requires wet-lab
    known-ratio controls (`../mixed_host_control_protocol.md`); `host_fractions_benchmarked` stays
    `false` until then.

## Honest scope (critical appraisal)

- **What this calibrates:** the *bioinformatic* detection floor and the minor-host fraction
  threshold for the read-classification denoising mode — the exact "denoising thresholds not
  benchmarked" gap.
- **What it does NOT:** extraction/PCR-amplification bias, mtDNA-copy differences (nucleated avian
  RBCs), digestion, or the read-fraction↔blood-volume relationship. Those still need the wet-lab
  known-ratio controls in [`../mixed_host_control_protocol.md`](../mixed_host_control_protocol.md).
- **Simulation caveat:** reads are simulated with an approximate ONT error model, not real R10.4.1;
  re-run with `--pool` on **real** pure single-host read pools (available once the full run yields
  cattle/goat/sheep-only samples) to confirm before locking the threshold. The broad RefSeq-mito
  database (18k taxa, many near-neighbours) may increase cross-mapping, which would make the 2 %
  false-positive suppression even more valuable.
