# In Silico Mixed-Host Validation and Pipeline Integration Report

**Date:** 2026-07-02 · **Branch:** `pi-led-reengineering`

## 1. Executive summary

`make_insilico_mixtures.py` was inspected, refined, run, and integrated as a **standalone
developer/regression tool** (not a production module). Deterministic synthetic fixtures now live in
`tests/data/insilico_mixtures/`; a fast, container-free **regression test** protects the mixed-host
detection logic; and explicit validation-status metadata (`host_fractions_benchmarked`,
`wetlab_known_ratio_validation`) is now emitted in the RAMBO summary, the run manifest, and the
mixture metadata. `run_full.sh` was reduced to a machine-specific wrapper with production defaults
moved into the `workstation` profile/config. No quantitative validation is claimed;
`host_fractions_benchmarked` remains **false**.

## 2. Status of `make_insilico_mixtures.py`

**It had already been run — this session** (commit `96fa0bb`), producing 180 synthetic mixtures and
the calibration in `docs/denoising_calibration/`, which set `rambo_min_host_fraction = 0.02`. It is
**not called by any Nextflow module** (verified). This task refined it (validation metadata + seed
recording), re-ran a small deterministic subset as committed fixtures, and added regression tests.

## 3. Script review

| Aspect | Finding |
|--------|---------|
| Purpose | Build known-ratio synthetic mixed-host read sets for computational detection testing |
| Inputs | `--simulate` (reference FASTA + primers) **or** `--pool NAME=fastq` (real pure pools) |
| Outputs | per-mixture FASTQ, `mixture_manifest.tsv` (truth), `run_metadata.json` (caveats + seed) |
| Determinism | Yes — `random.Random(--seed)`; seed recorded in metadata |
| Read provenance | Read IDs `MAJ_<host>_/MIN_<host>_` encode source host |
| Ratios / hosts | Configurable; single-host (`100:0`) + two-/multi-host supported |
| gzip | Output `.gz` supported; fixtures stored gzipped |
| Limitation | Simulated reads use an approximate ONT error model + insert-only amplicons (no primers), so they feed the read-classification harness, **not** the full primer-splitting pipeline as-is |

## 4. Revised validation plan

**A — validated now (in silico):** mixed-host detection logic, threshold behaviour, known input
read-ratio recovery *as read fraction only*, vector exclusion, output schema, and regression
stability. Result: LOD ≈ 5% minor host; `min_host_fraction = 0.02` (0% false positives).

**B — NOT validatable in silico:** blood-volume quantification, extraction/PCR/primer bias,
digestion-dependent detectability, and the true read-fraction↔blood-proportion relationship.

**C — future wet-lab plan:** single-host + two-host known-ratio mixtures across 95:5…5:95, marker-
specific, replicated, extraction-to-sequencing, with positive/negative controls, and explicit
criteria to flip `host_fractions_benchmarked` — fully specified in
[`mixed_host_control_protocol.md`](mixed_host_control_protocol.md).

## 5. Mixtures generated (committed fixtures)

| Mixture ID | Source hosts | Ratio | Depth | Output | Seed |
|-----------|--------------|-------|-------|--------|------|
| Homo_sapiens__Bos_taurus__100-0__r0 | Human (pure) | 100:0 | 300 | `.fastq.gz` | 42 |
| Homo_sapiens__Bos_taurus__90-10__r0 | Human+Cattle | 90:10 | 300 | `.fastq.gz` | 42 |
| Homo_sapiens__Bos_taurus__50-50__r0 | Human+Cattle | 50:50 | 300 | `.fastq.gz` | 42 |
| Bos_taurus__Homo_sapiens__100-0__r0 | Cattle (pure) | 100:0 | 300 | `.fastq.gz` | 42 |
| Bos_taurus__Homo_sapiens__90-10__r0 | Cattle+Human | 90:10 | 300 | `.fastq.gz` | 42 |
| Bos_taurus__Homo_sapiens__50-50__r0 | Cattle+Human | 50:50 | 300 | `.fastq.gz` | 42 |

Manifest + `run_metadata.json` (seed, hosts, ratios, caveats) alongside. Regenerate:
`docs/denoising_calibration/README.md`.

## 6. Integration decision

**Standalone developer/regression utility** (prompt options A + C), **not** a production module —
synthetic mixture generation is never part of routine runs. Consumed by a pytest regression test at
the detection-logic level (fast, deterministic, no containers/Nextflow). A full-pipeline test-profile
path (primer-inclusive simulation + `nextflow run -profile test`) is documented as a future option.

## 7. Files changed

| File | Change | Reason |
|------|--------|--------|
| `bin/make_insilico_mixtures.py` | `run_metadata.json` with seed + validation-status caveats | carry caveats with data |
| `bin/classify_reads.py` | extract testable `denoise_and_call()` | unit-testable detection logic |
| `bin/rambo_mixed_model.py` | add `wetlab_known_ratio_validation=false` | prevent overclaiming |
| `bin/build_run_manifest.py` | add `validation_status` block | explicit status in manifest |
| `nextflow.config` | `workstation` profile now sets `haema-*` containers | defaults not hidden in run_full.sh |
| `run_full.sh` | reduced to machine-specific wrapper | core defaults live in profile/config |
| `docs/usage.md` | out-of-the-box + broad-run commands | runnable without hidden flags |
| `tests/test_insilico_mixed_detection.py` | new (8 tests) | regression protection |
| `tests/data/insilico_mixtures/` | committed fixtures | reproducible test data |

## 8. New/updated outputs

| Output | Purpose | Produced by |
|--------|---------|-------------|
| `run_metadata.json` (mixtures) | seed + validation caveats | `make_insilico_mixtures.py` |
| RAMBO summary `wetlab_known_ratio_validation` | explicit non-quantitation | `rambo_mixed_model.py` |
| manifest `validation_status` | explicit status block | `build_run_manifest.py` |

## 9. Test assertions added

Single-host → **not** mixed; two-host (10%) → **mixed**; minor at 1% (below 2% floor) → **dropped**;
minor at ~2% → **detected**; **vector (Anopheles) excluded**; global single-read filter; PAF
identity/coverage filtering; fixture metadata asserts `host_fractions_benchmarked == false`.

## 10. `run_full.sh` review

It was a machine-specific wrapper that **also hid production defaults** (containers, curated-check,
Medaka). Moved: `haema-python/r/figures` → `workstation` profile. Kept in the wrapper (correctly,
machine-specific): `--input`, `--raw_data_dir`, `--reference_mode broad_blast`, `--blast_db`,
`--blast_db_mount`, `--outdir`. `broad_blast`+db stays a per-run choice (can't be a safe global
default — it needs a database). Effective params unchanged → `-resume` cache preserved.

## 11. Reference/database handling

Bundled: curated vertebrate panel + taxonomy sidecar + primers (small, in `assets/`). User-supplied:
the broad BLAST database (too large to commit) — built reproducibly by `bin/build_reference_db.py`
(RefSeq mitochondrion) with a provenance JSON + SHA-256, pointed to via `--blast_db`. Documented in
`usage.md`, `reproducibility.md`, `CONTAINER_STRATEGY.md`.

## 12. Tests run

| Command | Result | Notes |
|---------|--------|-------|
| `compileall bin/` | PASS | |
| `make_insilico_mixtures.py --help` / `classify_reads.py --help` | PASS | |
| `pytest tests/` | **PASS (24)** | +8 new regression tests |
| `nextflow config -profile workstation,docker` / `-profile test` | PASS | workstation resolves haema-* |
| `nextflow run . -profile test` | **DEFERRED** | full run is using `pipeline/work`; running it now would contend the cache. Run after the full run finishes (or with a separate `-w`). |

## 13. Remaining limitations

- **No wet-lab known-ratio controls**; `host_fractions_benchmarked = false`,
  `wetlab_known_ratio_validation = false`.
- In-silico mixtures validate **computational detection only** — not blood-volume quantification or
  extraction/PCR/primer/degradation bias.
- Threshold calibrated on **simulated** reads; re-confirm with `--pool` on real single-host pools.
- Synthetic reads are insert-only (no primers) → they feed the classification harness, not the full
  primer-splitting pipeline, without a primer-inclusive simulation refinement.

## 14. Next actions

1. Finish the full run; `-resume` to apply the 0.02 threshold to final host calls.
2. Re-run the calibration harness on **real** single-host pools once available.
3. Run `nextflow run . -profile test` after the full run frees `pipeline/work`.
4. (Optional) primer-inclusive simulation + a full test-profile mixed-host assertion.
5. Wet-lab known-ratio validation per the protocol to ever flip the benchmarked flag.
