# HÆMA Literature-Guided Re-engineering — Implementation Report

**Date:** 2026-07-02
**Branch:** `reengineer-reference-taxonomy-host-use`
**Pipeline:** HÆMA (`haema`, `main.nf`, Nextflow DSL2)
**Source review:** `docs/literature_guided_reengineering_report.md`

> This report supersedes an earlier draft of the same name that described the taxonomy
> Modes A–D as implemented when they were only config/documentation scaffolding. It records
> what is *actually* wired and how it was verified.

---

## 1. Executive summary

The re-engineering removes the curated vertebrate panel's status as the *sole* limiting reference
strategy and rebuilds the taxonomy workflow around a canonical four-mode reference architecture
(`--reference_mode`), while preserving the conservative, LCA-based host-call logic and every
scientific-framing safeguard. The previously orphaned analysis scripts and provenance code are now
fixed to the real data schema, wired into the workflow as gated steps, and verified.

**Pipeline status: runnable.** A full `nextflow run -profile test,docker -stub-run` completes
(exit 0) with the new steps executing, and every reference mode produces a valid DAG under
`-preview`. What has **not** been done is a full *containerised, real-data* run (no BLAST/Medaka
binaries or `nt`/BOLD databases are available in the dev environment) — so numerical outputs from
live BLAST are unverified. See §10–11.

The work was committed in reviewable increments; the starting state (a partial, over-claiming prior
session) was captured as an honest checkpoint first, and a full working-tree snapshot was archived
before any change.

---

## 2. Major design decision — reference architecture

**Problem.** The curated panel (a small vertebrate mitogenome set) was the dominant route to host
identification, which the literature (Channumsin 2021; Kipp 2023; Santos 2019) identifies as the
primary limiting factor. A prior session added a *parallel, inert* `reference_mode` config surface
that no `.nf` consumed, duplicating the existing, working `taxonomy_strategy` engine.

**Resolution.** `--reference_mode` is now the **canonical** selector and is actually consumed by the
taxonomy subworkflow. The legacy `--taxonomy_strategy` is retained as a **backward-compatible alias**
that maps onto a reference mode when `--reference_mode` is left at its default, so existing configs
and tests are unchanged. The curated panel is **kept** — repositioned as an optional first-line
reference (Modes A/B pre-check and the primary in C/D), not the only route. Each mode uses at most
one fallback database, which keeps `parse_blast_assignments.py`'s single-fallback contract intact:

| Mode | `--reference_mode` | Behaviour |
|------|--------------------|-----------|
| A | `curated_panel` | Curated panel only (fast, offline). |
| B | `broad_blast` | User `--blast_db`; curated pre-check on by default, else broad-only. |
| C | `remote_fallback` | Curated primary → NCBI `nt` via `blastn -remote` (gated; not reproducible). |
| D | `bold_aware` | Curated primary → reproducible local BOLD-derived COI FASTA (`--bold_fasta`). |

BOLD is implemented as a **reproducible local FASTA** built into a BLAST DB at runtime; a live BOLD
API (`--bold_mode api_query`) is deliberately not implemented so routine runs stay reproducible.

---

## 3. Files changed

| File | Change | Reason |
|------|--------|--------|
| `subworkflows/local/taxonomy/main.nf` | Rewritten to resolve + branch on `reference_mode` (Modes A–D); single-fallback join | Report §4.3, task B5 |
| `modules/local/blast/external/main.nf` | New `remote` input → conditional `blastn -remote` | Mode C |
| `modules/local/blast/makeblastdb/main.nf` | Generalised with `db_tag` (name + publishDir) | Second (BOLD) DB without collision |
| `modules/local/blast/blastn/main.nf` | New `BLASTN_BOLD_ASVS` (distinct `.bold.blast.tsv`) | Mode D runtime BOLD DB |
| `modules/local/taxonomy_assign/main.nf` | `ASSIGN_TAXONOMY_WITH_FALLBACK` takes fallback-source label; all assign procs pass `--marker` + thresholds | Provenance + per-marker guard |
| `bin/parse_blast_assignments.py` | `--marker`/`--coi-species-identity`/`--cytb-species-identity`; conservative high→medium downgrade | Report §4.3–4.4 |
| `bin/compute_marker_concordance.py` | Fixed to real RAMBO schema (`host_assignment`/`best_confidence`), rank-1 per marker, control exclusion, always-header | Task B2 (was silently empty) |
| `bin/compute_host_ecology_comparisons.py` | Rewrote data layer: per-`sample_uid` host-sets, strata from master endpoint, season HBI | Task B1 (was double-counting) |
| `modules/local/marker_concordance/main.nf` | **New** module | Wire B2 |
| `modules/local/host_ecology_comparisons/main.nf` | **New** module | Wire B1 |
| `modules/local/publication_figures/main.nf` | Optional concordance input → figure S3 heatmap | Task B3 |
| `modules/local/run_manifest/main.nf` | Serialise new taxonomy params into parameters.json | Wire provenance/NUMT (A3/B7) |
| `main.nf` | Invoke `MARKER_CONCORDANCE`, `HOST_ECOLOGY_COMPARISONS`; pass concordance to figures | Wire B1/B2/B3 |
| `nextflow.config`, `nextflow_schema.json` | `remote_blast_db`, `enable_marker_concordance`, `enable_host_ecology_comparisons` | New params |
| `tests/test_marker_concordance.py`, `tests/test_host_ecology_comparisons.py` | **New** unit tests (11 cases) | Verify B1/B2 |
| `docs/methods.md`, `usage.md`, `parameters.md`, `output.md` | Document modes, params, outputs | Report §9 |
| `tests/*` | Restored from an accidental relocation into `assets/tests/` | Make-honest |

---

## 4. New parameters

| Parameter | Purpose | Default | Notes |
|-----------|---------|---------|-------|
| `--reference_mode` | Canonical reference selector (A–D) | `curated_panel` | Maps legacy `--taxonomy_strategy` when default |
| `--remote_blast_db` | NCBI remote DB name (Mode C) | `core_nt` | `blastn -remote`; not reproducible |
| `--enable_marker_concordance` | Run multi-marker concordance | `true` | Needs RAMBO host calls |
| `--enable_host_ecology_comparisons` | Run exploratory Fisher/Holm comparisons | `true` | Needs RAMBO host calls |

Params added by the prior session and now **actually consumed**: `--enable_curated_panel_check`,
`--enable_ncbi_remote_fallback`, `--enable_bold_fallback`, `--bold_mode`, `--bold_fasta`,
`--bold_taxonomy`, `--coi_species_identity_threshold`, `--cytb_species_identity_threshold`,
`--marker_numt_risk`.

---

## 5. New or modified outputs

| Output | Description | Source |
|--------|-------------|--------|
| `05_endpoint_files/marker_concordance.tsv` (+ `06_reports/…_summary.json`) | Per-specimen species/genus agreement + NUMT/mixed caution flags | `MARKER_CONCORDANCE` |
| `05_endpoint_files/host_ecology_comparisons/*.tsv` (+ summary JSON) | Fisher/Holm HBI, HBI-by-species, mixed-feeding, host richness | `HOST_ECOLOGY_COMPARISONS` |
| `run_manifest.json` → `reference_database` | Mode, per-DB SHA-256 provenance, fallback chain, per-marker thresholds + NUMT risk | `build_run_manifest.py` |
| `supplementary/figure_S3_concordance_heatmap` | Specimen × marker concordance heatmap | `PUBLICATION_FIGURES` (optional) |

---

## 6. Taxonomy workflow after redesign

```
consensus/ASV sequences
  └─ reference_mode selector (legacy taxonomy_strategy mapped in)
       ├─ A curated_panel  : curated BLAST ─────────────► conservative LCA
       ├─ B broad_blast    : [curated pre-check] + broad --blast_db (fallback) ► LCA
       ├─ C remote_fallback: curated + NCBI nt (blastn -remote, gated) ────────► LCA
       └─ D bold_aware     : curated + local BOLD COI FASTA (runtime DB) ──────► LCA
                                   │ (marker species-identity guard: COI 98 / CytB 97)
                                   ▼
                         AGGREGATE_RESULTS ─► RAMBO host-call model
                                   ├─► ECOLOGICAL_INDICES (Wilson CIs; controls excluded)
                                   ├─► MARKER_CONCORDANCE (NUMT/mixed caution flags)
                                   └─► HOST_ECOLOGY_COMPARISONS (Fisher + Holm; exploratory)
```

---

## 7. Reference database provenance

`run_manifest.json.reference_database` now records, per run: the resolved `mode`; a `provenance`
block per database used (curated_panel, broad/fallback, BOLD, taxdump) with name, type, source,
path, **SHA-256 checksum**, sequence count where applicable, and taxonomy-map/taxdump paths; the
`fallback_chain`; the per-marker `identity_thresholds`; and the per-marker `numt_risk`. Missing
provenance raises warnings rather than being silently treated as complete (e.g. remote fallback
enabled without a fallback DB). Verified end-to-end against a fixture manifest.

---

## 8. Marker concordance and NUMT reporting

`MARKER_CONCORDANCE` reduces each specimen's RAMBO host calls to one dominant (rank-1) host per
marker, excludes controls, and classifies agreement (`full_species_agreement`, `genus_agreement`,
`discordant`, `single_marker_only`, `no_marker_signal`, `ambiguous_lca_only`), with
`possible_numt_flag` raised on CytB↔COI disagreement and `possible_mixed_meal_flag` on species-level
disagreement. **NUMT risk is a caution flag, not proof** and never a hard filter; per-marker labels
(co1_short = moderate, co1_long = low, cyt_b = low) are configurable via `--marker_numt_risk` and
reported in the manifest and concordance outputs.

---

## 9. Statistical and visualisation updates

`HOST_ECOLOGY_COMPARISONS` computes **exploratory** pairwise Fisher's exact tests (Holm-corrected
within family) for HBI between zones, HBI between seasons, HBI between sibling species, and
mixed-feeding rate between zones, plus a descriptive host-richness table — all on host-identified
field specimens with a mosquito counted **once** (host-set per `sample_uid`), controls excluded,
consistent with `compute_ecological_indices.py`. scipy gives the exact test; a chi-square
approximation (n ≥ 10) is the documented fallback. Small-n strata (< 5) are flagged; outputs are
labelled exploratory. Visualisation: figure S3 concordance heatmap is wired (degrades cleanly when
concordance is off). Database-source contribution and NUMT-summary plots remain future work (§11).

---

## 10. Tests run

| Command | Result | Notes |
|---------|--------|-------|
| `python -m compileall bin/` | **PASS** | All bin scripts compile |
| `pytest tests/` | **PASS** (16) | 8 restored + 2 new suites (concordance, host-ecology) |
| `tests/test_lca.py`, `tests/test_taxid_assignment.py` | **PASS** | LCA unchanged by marker guard |
| parse_blast marker-guard functional check | **PASS** | COI 97.5%→medium, 98.5%→high |
| manifest provenance/NUMT functional check | **PASS** | mode, per-DB SHA-256, NUMT labels populate |
| `nextflow config -profile test` | **PASS** | Config + new params resolve |
| `nextflow run -preview -stub` (A/B/broad+curated/C/D) | **PASS** | Valid DAG + correct topology per mode |
| `nextflow run -profile test,docker -stub-run` | **PASS (exit 0)** | Completes; `MARKER_CONCORDANCE` + `HOST_ECOLOGY_COMPARISONS` execute |
| Full containerised real-data run | **NOT RUN** | No BLAST/Medaka binaries or nt/BOLD DBs in dev env |
| `nf-test`, `nf-core lint` | **NOT RUN** | Not installed; pipeline is not nf-core-structured |

---

## 11. Remaining limitations

- **No live-BLAST verification.** Numerical taxonomy outputs depend on BLAST against real databases,
  which could not run here. Modes are DAG- and stub-verified only; run `-profile test` with real
  containers before trusting numbers.
- **Mode fallback chains are single-tier.** curated→broad, curated→nt-remote, curated→BOLD. A
  3-tier curated→nt→BOLD chain is not implemented (documented as future work).
- **BOLD** is local-FASTA only (reproducible). Live API (`--bold_mode api_query`) is intentionally
  unimplemented. BOLD fallback runs on all markers but is effective only for COI.
- **Remote (Mode C)** is not reproducible by default; `database_version` for `nt` is recorded as
  "unknown". Reserve for exploratory look-ups.
- **Cohen's kappa** in concordance is a simplified measure (documented); weighted kappa is future work.
- **Deferred:** mixed-host positive-control test-samplesheet entry (A4), haemavec-figures digest pin
  (B4), Kruskal-Wallis host-richness test (C3), `digestion_class` metadata (C2), database-source and
  NUMT-summary plots.
- The test profile still defaults `skip_taxonomy = true`; pass `--skip_taxonomy false` to exercise
  the taxonomy path in CI.

---

## 12. Next recommended actions

1. Run `nextflow run . -profile test,docker --skip_taxonomy false` (real containers) and confirm
   taxonomy + concordance + ecology outputs on the bundled data.
2. Do one real-data Mode A run, then a Mode B run against a local `nt`/vertebrate-mtDNA DB, and
   diff host calls to quantify the coverage gain.
3. Add a 2-host positive-control row (e.g. *Homo sapiens* + *Bos taurus*) to the test samplesheet (A4).
4. Prepare a versioned, checksummed BOLD-derived COI FASTA and validate Mode D on it.
5. Supply `--taxdump_dir` (pinned NCBI snapshot) for taxid-backed LCA in production.
6. After first registry push, pin the `haema-figures` image digest (B4).
7. Add the database-source contribution plot and Kruskal-Wallis host-richness test.

---

*End of Implementation Report*
