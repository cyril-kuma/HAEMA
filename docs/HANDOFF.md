# HÆMA — session handoff / resume guide

**Last updated:** 2026-07-02 · **Active branch:** `pi-led-reengineering` (19 commits ahead of `main`)
**Repo:** `pipeline/` is the git repo (remote `origin` → github.com/cyril-kuma/HAEMA.git).
The outer project dir (`01_Bloodmeal_Host_Metabarcoding/`) is **not** a git repo.

This file lets another session/device resume without re-deriving context. Read it first.

---

## 1. Where things stand

**Done and committed on `pi-led-reengineering`:**
- **Taxonomy redesign** — `--reference_mode` Modes A–D (curated / broad_blast / remote / bold_aware);
  validated on real data in `broad_blast` mode.
- **Broad reference DB built** — NCBI RefSeq mitochondrion r235 (18,198 seqs / 18,026 taxa) at
  `data/reference_db/refseq_mito_r235/` (**outside git**, ~570 MB). Built by `bin/build_reference_db.py`.
- **Provenance + NUMT** in `run_manifest.json`; **vector self-hit exclusion** (`--non_host_genera`);
  **marker concordance** + **host-ecology comparisons** (Fisher+Holm, Kruskal-Wallis) wired.
- **Deferred items** A4 (real mixed controls validated), B4 (container digest-pin procedure), C3
  (Kruskal-Wallis) resolved; C2 (digestion metadata) documented-deferred.
- **Denoising calibration bridge** — in-silico known-ratio mixtures → **`rambo_min_host_fraction`
  default raised 0.01 → 0.02** (calibrated: 0% false positives, 100% detection ≥5% minor host).
- **Docs** current: methods, limitations, parameters, output, usage, reproducibility,
  CONTAINER_STRATEGY, architecture_review, samplesheet_preparation, mixed_host_control_protocol,
  denoising_redesign_plan, denoising_calibration/, pi_reengineering_report.

**In flight:** the **full 714-sample run** (`results/full_run_v2/`) — see §3.

**Recovery snapshot:** a pre-work tar of the tree is in the session scratchpad (not the repo).

## 2. Key paths, params, tools

| Thing | Value |
|-------|-------|
| Reference DB (broad) | `data/reference_db/refseq_mito_r235/blastdb/refseq_mito_r235` (prefix) |
| Full samplesheet | `metadata/full_run_samplesheet.csv` (714 rows with data) |
| Subset samplesheet | `metadata/subset_validation_samplesheet.csv` (27: 13 pos/8 neg/6 field) |
| Full-run launcher | `pipeline/run_full.sh` (broad_blast + RefSeq mito, Medaka on) |
| Subset launcher | `pipeline/run_subset_validation.sh` (Medaka off, fast) |
| Nextflow engine | `export NXF_VER=24.10.5` (24.10 LTS; newer strict parser breaks this pipeline) |
| Profile | `-profile workstation,docker` (16 CPU / 56 GB; GPU Medaka via docker profile) |
| Custom images | `haema-python:0.3.0`, `haema-r:0.3.0`, `haema-figures:0.4.0` (local build; not pushed) |
| Calibrated host-call threshold | `rambo_min_host_fraction = 0.02` (default, out of the box) |

**Standalone dev/validation tools (NOT in the pipeline, run manually):**
`bin/build_reference_db.py`, `bin/make_insilico_mixtures.py`, `bin/classify_reads.py`,
`bin/benchmark_mixed_detection.py`. They do **not** use or change the samplesheet.

## 3. Resume the full run

The run is detached (`nextflow` JVM). Check and resume:

```bash
cd pipeline
pgrep -af nextflow-24                 # is it still running?
tail -f ../logs/full_run_v2/nextflow.log
# if stopped or a task failed, just re-launch — it -resumes from cache:
bash run_full.sh
```

**Apply the calibrated 0.02 threshold to the full run:** the in-flight run was launched with the old
`0.01` (Nextflow fixes params at launch). After it completes, run `bash run_full.sh` once more —
`-resume` re-runs **only** `RAMBO_MIXED_MODEL` + downstream (ecology, concordance, comparisons,
report, figures) with `0.02`, reusing all cached upstream (cheap, minutes).

**On failure:** read `.command.err`/`.command.sh` in the failed task's `work/` dir (path printed in
the log), fix, `bash run_full.sh` to resume. Do not delete `work/`.

**Final outputs:** `results/full_run_v2/05_endpoint_files/` — `bloodmeal_master_endpoint.tsv`,
`host_call_table.tsv`, `marker_concordance.tsv`, `host_ecology_comparisons/`, `run_manifest.json`.

## 4. Re-run the in-silico denoising calibration (optional)

Standalone; see `docs/denoising_calibration/README.md`. Once the full run yields **real** cattle/
goat/sheep-only samples, re-run with `make_insilico_mixtures.py --pool NAME=reads.fastq …` (real
pools) instead of `--simulate`, to confirm the 0.02 threshold on real reads before locking it.

## 5. Outstanding / next steps (priority order)

1. **Finish the full run**, then `-resume` to apply 0.02 (§3).
2. **Confirm 0.02 on real single-host pools** (§4) before treating it as final.
3. **(Optional, larger) wire `read_classification`/`hybrid` mode** into Nextflow per
   `denoising_redesign_plan.md` §5 (module + subworkflow, default-off) — only in-silico-validated so
   far; do not make default until real-data validated.
4. **Container reproducibility (B4):** push `haema-*` images to a registry, pin by `@sha256:`.
5. **Wet-lab quantitative validation** (`mixed_host_control_protocol.md`) to ever flip
   `host_fractions_benchmarked` to `true`.
6. **Figure-system consolidation** (`architecture_review.md`) — deferred.
7. Merge `pi-led-reengineering` → `main` when the full run + review are complete.

## 6. Known limitations (see `docs/limitations.md` for the full, current list)

- Denoising LOD ≈ 5% minor host; threshold calibrated on **simulated** reads (confirm on real).
- Host fractions are support evidence only, **not** blood volume (`host_fractions_benchmarked:false`).
- Broad DB is all-taxa RefSeq mito → vector/microbial matches excluded by genus list, not by a
  vertebrate-restricted DB.
- Curated panel not yet checksum-versioned; no bundled taxdump snapshot; `read_classification` mode
  not yet wired; live BOLD API intentionally unimplemented; `digestion_class` metadata deferred.
