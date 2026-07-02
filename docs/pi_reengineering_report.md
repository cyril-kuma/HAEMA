# PI-Led HÆMA Pipeline Re-engineering Report

**Date:** 2026-07-02 · **Branch:** `pi-led-reengineering` · **Pipeline:** HÆMA (Nextflow DSL2)

## 1. Executive summary

Acting as PI, I took the reviewed pipeline and (a) moved host identification off a
curated-panel-first strategy onto a **broad, reproducible BLAST database**, (b) resolved the
outstanding deferred items, (c) built real test + full samplesheets from the actual data,
(d) **validated the whole pipeline end-to-end on real data**, and (e) **launched the full
714-sample run**, which is progressing.

The end-to-end subset validation (27 samples: 13 positive incl. 3 mixed-fed, 8 negative, 6 field)
**completed green** in `broad_blast` mode against a purpose-built NCBI RefSeq mitochondrion (release
235) database. All single-host positive controls recovered their expected host at species level; the
mixed-fed controls recovered multiple hosts; the broad database resolved 24 features the curated
panel could not. This validation also surfaced — and I then handled — a genuine scientific
consequence of the broad database (vector self-amplification; §8).

## 2. Documentation synthesis (Phase 1)

`docs/` was reviewed and is current and useful (no stale debris): methods, usage, parameters,
output, limitations, ecological_indices, figures, reproducibility, CONTAINER_STRATEGY, benchmarking,
README, plus the re-engineering artifacts. Actionable themes and their resolution:

| Theme | Action taken |
|-------|--------------|
| Curated panel too limiting | Repositioned as optional pre-check/validation; broad_blast (RefSeq mito) is now the recommended real-run mode (`methods.md`, `usage.md`) |
| Reference provenance | `run_manifest.json.reference_database` records mode, per-DB SHA-256, fallback chain, thresholds, NUMT risk |
| Samplesheet guidance | New `docs/samplesheet_preparation.md` |
| Container pinning (B4) | Build→push→pin procedure + recorded local image IDs in `CONTAINER_STRATEGY.md` |
| Architecture clarity | New `docs/architecture_review.md` (bin/module review + stale audit) |

## 3. Live BOLD API decision (Phase 3)

**Decision: local, versioned database — live BOLD API deliberately not used by default.** The prior
review deferred live BOLD for **reproducibility** (a network call returns different results over
time), API instability/rate-limits, and the lack of taxid-compatible output for LCA. I confirmed
this reasoning and implemented the reproducible alternative: `broad_blast`/`bold_aware` consume a
**fixed, checksummed local database** built by `bin/build_reference_db.py`. For the real run I built
that database from **NCBI RefSeq mitochondrion release 235** (18,198 sequences / 18,026 taxa, COI +
CytB), which is broad enough to catch unexpected vertebrate hosts while remaining offline and
version-pinnable. `bold_aware` accepts a user-supplied BOLD-derived COI FASTA by the same mechanism;
`--bold_mode api_query` is reserved and documented as future work, never a default.

## 4. Deferred items audit (Phase 8)

| Item | Meaning | Original reason | Final decision | Action taken |
|------|---------|-----------------|----------------|--------------|
| A4 | Mixed-host positive control | No bundled control | **Resolved** | Real mixed-fed controls (`Homo sapiens;Capra hircus`, `;Bos taurus`) exist in the samplesheet and are in the validation subset; run confirms `mixed_pass_with_extra` / `mixed_partial`. Illustrative template also present. |
| B4 | Figures container digest pin | Image not registry-pushed | **Documented** | Hard-pinning a local image ID is non-reproducible; added exact build→push→pin procedure + recorded local image IDs (`CONTAINER_STRATEGY.md`). |
| C2 | `digestion_class` metadata column | Not in samplesheet | **Deferred** | No digestion/Sella field in the real data; documented as optional future enhancement. |
| C3 | Kruskal-Wallis host richness | Not implemented | **Implemented** | Added per-specimen host-richness Kruskal-Wallis across zones to `compute_host_ecology_comparisons.py` (exploratory; degenerate-case safe). |

## 5. Pipeline architecture review (Phases 2, 4, 11)

See `docs/architecture_review.md`. Findings: clean ~1:1 script↔module mapping; shared helpers
(`fastq_utils`, `figure_style`) correctly factored; **no hard-coded paths in reusable code**; **no
dead tracked scripts**. No modules merged (each is a genuine, separately-cached process boundary).
The one real consolidation opportunity — two overlapping figure systems (`figures` vs
`publication_figures`) — is documented as a deferred, non-blind release task (figure rendering can
only be validated by running the containerised steps). Untracked user scratch scripts were left
untouched.

## 6. `bin/` review — see `docs/architecture_review.md` (full table).

Highlights: `parse_blast_assignments.py` (conservative LCA + new marker guard) kept as core;
`build_reference_db.py` added; `compute_marker_concordance.py` / `compute_host_ecology_comparisons.py`
fixed to the real RAMBO schema and wired; figure scripts hardened. Nothing deleted.

## 7. Taxonomy / reference redesign (Phases 5–7)

`--reference_mode` is the canonical selector (Modes A curated_panel / B broad_blast / C
remote_fallback / D bold_aware); legacy `--taxonomy_strategy` maps onto it. Conservative assignment
(identity/coverage/e-value/top-bitscore-delta, LCA, per-marker species guard, confidence labels,
fallback + source tracking) is preserved and recorded per feature. The real run uses **Mode B with
the curated panel as pre-check + the RefSeq mito broad database** — curated resolves expected
peridomestic hosts (taxid-aware), the broad DB resolves the rest. Provenance (mode, per-DB SHA-256,
fallback chain, thresholds, NUMT risk) is written to `run_manifest.json`.

## 8. Key scientific finding — vector self-amplification

The broad database revealed the **vector's own mtDNA co-amplifying** (*Anopheles gambiae* returned
for 24 subset features) — impossible with the vertebrate-only curated panel. Blood-meal hosts are
vertebrates, so I added a **configurable, genus-aware non-host (vector) exclusion**
(`--non_host_genera`; default *Anopheles/Culex/Aedes/…* + ticks/sandflies/tsetse) applied to the
ecological indices, concordance, and comparisons. Vector calls are dropped from host-use summaries
but stay visible in the raw endpoint tables and the positive-control check for transparency/QC.
Documented in `limitations.md`; a vertebrate-restricted DB is noted as a cleaner long-term option.

## 9. Reference provenance (Phase 7)

`run_manifest.json.reference_database` records: resolved `mode`; per-database provenance
(name/type/source/path/SHA-256/sequence count) for the curated panel and broad DB; the `fallback_chain`;
per-marker identity thresholds; and per-marker NUMT risk. The broad DB is additionally described by
its own `*.provenance.json` from `build_reference_db.py` (source, release 235, build date, SHA-256,
counts). Verified in the subset `run_manifest.json` (`mode: broad_blast`; curated_panel + broad_blast_db).

## 10. Test samplesheet and assertions (Phase 9)

`metadata/subset_validation_samplesheet.csv` (27 rows). Assertions verified by the run:
pipeline completes; all endpoint files produced; **single-host controls → species match (`pass`)**;
**mixed controls → multiple hosts recovered**; negative controls produce no valid host call; broad
DB source recorded and contributing (24 fallback features); concordance/LCA/confidence populated;
figures + report generated; vector excluded from indices (0 *Anopheles* in concordance).

## 11. Full samplesheet (Phase 9)

`metadata/full_run_samplesheet.csv` (714 rows) = the master samplesheet filtered to rows whose raw
data is present on disk (2 controls dropped for missing data, reported). Reproducible by another
user via `docs/samplesheet_preparation.md`. Not machine-specific beyond `--input`/`--raw_data_dir`.

## 12. Figures/reporting refactor (Phase 10)

No blind merge (see §5). Concrete hardening: per-figure try/except resilience in
`build_main_figures.py` and `build_supplementary_figures.py` (a degenerate stratum no longer aborts
the whole figure step); empty-legend guard. Consolidation of the two figure systems recommended for
release.

## 13. Stale files cleaned (Phase 11)

None deleted. Audit found no stale tracked files; only untracked user scratch
(`build_thesis_figures.py`, `figure3B_options.py`, `run_full.sh` variants) which were left in place.

## 14. Tests run

| Command | Result | Notes |
|---------|--------|-------|
| `compileall bin/` | PASS | all scripts compile |
| `pytest tests/` | PASS (16) | incl. concordance + host-ecology suites |
| vector-exclusion + C3 functional checks | PASS | `is_host`, `per_mosquito_hosts`, `kruskal_richness` |
| `nextflow config -profile test` | PASS | new params resolve |
| `nextflow -preview -stub` (all modes) | PASS | valid DAG per reference mode |
| **Subset validation** `-profile workstation,docker`, real data, broad_blast + RefSeq mito | **PASS (exit 0)** | controls validated; vector excluded; §10 |

## 15. Full run status

**LAUNCHED and progressing** (`run_full.sh`; `results/full_run_v2/`, `logs/full_run_v2/`).
Preflight + Medaka-model preflight ✔; **all 714 samples validated**; preprocessing underway
(reusing the 27 cached subset barcodes via `-resume`). Configuration: `broad_blast` (curated
pre-check + RefSeq mito r235), Medaka on, `workstation,docker`. Expected multi-hour runtime
(CPU/GPU Medaka over ~700 barcodes). Monitoring: the run halts and reports on any task failure;
re-run `bash pipeline/run_full.sh` to `-resume` after fixing any blocking task. As of this report it
had produced no errors; final completion is pending.

## 16. Remaining limitations

- Full 714-sample run not yet complete at time of writing (launched, progressing).
- Broad DB is all-taxa RefSeq mito; vector/microbial matches are excluded by genus list, not by a
  vertebrate-restricted database (cleaner long-term option; §8).
- `broad_blast_db` provenance records the DB path but not a file SHA-256 (a BLAST DB is a multi-file
  prefix); the FASTA build has a SHA-256 in its provenance JSON.
- Remote mode (C) remains non-reproducible by design; live BOLD API not implemented.
- Two figure systems not yet consolidated; figure resilience mitigates fragility.
- C2 (`digestion_class`) deferred (no data).

## 17. Recommended next steps

1. Let the full run finish; inspect `results/full_run_v2/05_endpoint_files/` (endpoint, host calls,
   concordance, ecology comparisons, `run_manifest.json`).
2. Review low-frequency/unexpected host calls surfaced by the broad DB (QC).
3. Consider a Chordata-restricted broad database build for the final analysis.
4. Push the three `haema-*` images to a registry and pin by digest (B4).
5. Consolidate the two figure systems into the `publication_figures` suite.
6. Supply `--taxdump_dir` (pinned NCBI snapshot) if taxid-backed LCA on the broad DB is wanted.

---
*End of PI-Led Re-engineering Report*
