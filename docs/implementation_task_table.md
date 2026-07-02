# HÆMA Literature-Guided Re-engineering — Implementation Task Table

Generated: 2026-06-30
Source: docs/literature_guided_reengineering_report.md

## Priority 1 — High impact, low effort (before thesis submission)

| # | Recommendation | Priority | Files Likely Affected | Implementation Type | Status |
|---|---------------|----------|----------------------|---------------------|--------|
| A1 | Update docs/methods.md: framing, per-marker identity context, RAMBO threshold comparison, UMAP+RAMBO layering | High | docs/methods.md | Documentation | NOT STARTED |
| A2 | Update docs/limitations.md: BOLD gap, NUMT risk, detection window, reference dependency | High | docs/limitations.md | Documentation | NOT STARTED |
| A3 | Add marker_numt_risk lookup to run_manifest.json output | High | bin/build_run_manifest.py | Code | NOT STARTED |
| A4 | Add mixed-host positive control entry to test samplesheet | High | assets/tests/ (or test data) | Testing | NOT STARTED |

## Priority 2 — Medium impact, moderate effort (before viva)

| # | Recommendation | Priority | Files Likely Affected | Implementation Type | Status |
|---|---------------|----------|----------------------|---------------------|--------|
| B1 | Add Fisher's exact + Holm correction for HBI pairwise stratum comparisons | Medium | bin/compute_host_ecology_comparisons.py (new), pipeline/modules/local/ (new), nextflow_schema.json, nextflow.config | Code + Nextflow module | NOT STARTED |
| B2 | Add multi-marker concordance analysis | Medium | bin/compute_marker_concordance.py (new), pipeline/modules/local/ (new), nextflow_schema.json, nextflow.config | Code + Nextflow module | NOT STARTED |
| B3 | Add concordance heatmap figure | Medium | bin/build_main_figures.py or bin/build_supplementary_figures.py | Visualisation | NOT STARTED |
| B4 | Pin haemavec-figures container digest in nextflow.config | Medium | nextflow.config | Config | NOT STARTED |
| B5 | Redesign taxonomy workflow: support Mode A (curated panel), Mode B (broad BLAST), Mode C (remote fallback), Mode D (BOLD-aware) | High | pipeline/subworkflows/local/taxonomy/main.nf, pipeline/modules/local/blast/, bin/parse_blast_assignments.py, nextflow_schema.json, nextflow.config | Nextflow module + Code | NOT STARTED |
| B6 | Add per-marker identity interpretation documentation | Medium | docs/methods.md, docs/output.md | Documentation | NOT STARTED |
| B7 | Add reference database provenance reporting | Medium | bin/build_run_manifest.py, bin/build_report.py, pipeline/modules/local/run_manifest/ | Code | NOT STARTED |
| B8 | Add QC matrix documentation | Medium | docs/methods.md | Documentation | NOT STARTED |

## Priority 3 — Enhancement (post-thesis or supervised extension)

| # | Recommendation | Priority | Files Likely Affected | Implementation Type | Status |
|---|---------------|----------|----------------------|---------------------|--------|
| C1 | Add BOLD API as tertiary COI fallback | Low | bin/bold_query.py (new), pipeline/subworkflows/local/taxonomy/main.nf, nextflow_schema.json | Code + Nextflow module | NOT STARTED |
| C2 | Add optional digestion_class metadata column | Low | nextflow_schema.json, validate_inputs.py, compute_ecological_indices.py | Config + Code | NOT STARTED |
| C3 | Add formal Kruskal-Wallis test for host richness | Low | bin/compute_host_ecology_comparisons.py | Code | NOT STARTED |

## Documentation-Only Recommendations

| # | Recommendation | Priority | Files Likely Affected | Implementation Type | Status |
|---|---------------|----------|----------------------|---------------------|--------|
| D1 | Document length-fallback marker misassignment risk | Medium | docs/methods.md | Documentation | NOT STARTED |
| D2 | Document UMAP denoising vs RAMBO layering | Medium | docs/methods.md | Documentation | NOT STARTED |
| D3 | Document RAMBO 1% vs Logue 10% threshold | Medium | docs/methods.md, rambo output summary | Documentation | NOT STARTED |
| D4 | Document per-marker identity interpretation (97% global, 98% BOLD COI standard) | Medium | docs/methods.md, run_manifest.json | Documentation | NOT STARTED |

## Summary by Category

| Category | Count | Items |
|----------|-------|-------|
| Documentation | 7 | A1, A2, D1, D2, D3, D4, B6, B8 |
| Code (Python) | 5 | A3, B1, B2, B7, C2, C3 |
| Nextflow module | 4 | B1, B2, B5, C1 |
| Config/schema | 4 | B4, B5, B7, C1, C2 |
| Visualisation | 1 | B3 |
| Testing | 1 | A4 |
| Reference/database | 1 | B5 (BOLD) |
| Release readiness | 1 | B4 |
| **Total** | **24** | |

## Implementation Order (refined)

### Phase 1: Documentation (no code changes)
1. A2 — Update docs/limitations.md
2. A1 — Update docs/methods.md
3. D1-D4 — Additional documentation items
4. B6 — Per-marker identity interpretation

### Phase 2: Core taxonomy redesign
5. B5 — Redesign taxonomy workflow (subworkflows/local/taxonomy/main.nf)
6. B7 — Add reference database provenance reporting
7. A3 — Add NUMT risk to run_manifest

### Phase 3: New analytical modules
8. B2 — Multi-marker concordance analysis
9. B1 — Host ecology comparisons
10. B3 — Concordance heatmap figure

### Phase 4: Testing and release readiness
11. A4 — Mixed-host positive control test
12. B4 — Pin container digest
13. Run tests

### Phase 5: Post-thesis enhancements
14. C1 — BOLD API integration
15. C2 — Digestion class metadata
16. C3 — Kruskal-Wallis test
