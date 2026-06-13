# Release checklist

Track readiness for each milestone. Checked = done in this repository; unchecked = remaining.

## Engineering / repository
- [x] Self-contained repo (references + primers vendored under `assets/`)
- [x] `launchDir`-relative outputs; no hardcoded absolute paths in code/config
- [x] Required inputs validated with fail-fast messages
- [x] All containers pinned by digest / immutable tag (no `:latest`)
- [x] `nextflow_schema.json` covers all parameters
- [x] `LICENSE`, `CITATION.cff`, `CHANGELOG.md`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`
- [x] GitHub issue/PR templates and a CI workflow
- [x] `docs/`: usage, parameters, output, methods, limitations, reproducibility, container strategy
- [x] `Makefile` + `tests/validate_release.sh` one-command build/lint/release entry points
- [x] Git repository initialised with commit history
- [ ] Set the real repository URL in `CITATION.cff` (repository-code) **and** `nextflow.config` (manifest.homePage) — both currently the `USER/haema` placeholder
- [ ] Fill remaining `CITATION.cff` placeholders (ORCID, affiliation, Zenodo DOI)
- [ ] Push `haema-python` / `haema-r` to a public registry and pin digests (`make push REGISTRY=...`)
- [ ] Tag a release (`v0.2.0`) and archive for a DOI (e.g. Zenodo)

## Testing
- [x] `nextflow config` parses for all 8 profiles
- [x] `-profile test` real run completes (curated BLAST taxonomy)
- [x] Whole-DAG stub run (all features) completes
- [x] Real-data demonstration run completes with mixed-host detection
- [x] Python unit tests: input validation + positive-control logic (`tests/`)
- [x] Reference-asset verification (checksums + sidecar fields); endpoint-column assertions
- [x] One-command pre-release validation (`bash tests/validate_release.sh` / `make lint`)
- [x] taxdump/LCA + curated-panel taxid-assignment unit tests (synthetic fixtures)
- [ ] CI green on GitHub-hosted runners (workflow + tests added; verify on first push)

## Scientific
- [x] Implemented vs staged features documented
- [x] Limitations and assumptions documented
- [x] Reference panel checksummed (`CHECKSUMS.sha256`) and field-validated
- [x] Single- and **mixed**-host controls auto-checked vs `expected_host_scientific_name`
- [x] Mixed-host recovery benchmark framework (`docs/benchmarking.md`) for lab-prepared controls
- [x] Host fractions explicitly flagged as evidence (not benchmarked) in report + manifest
- [x] taxdump-backed LCA works for the curated panel via sidecar taxids (tested)
- [ ] Denoising thresholds benchmarked against the lab-prepared mixed controls (declare composition + sweep)
- [ ] Curated panel expanded with more West African sylvatic taxa
- [ ] Quantitative host-fraction validation (needs controls at known mixing ratios)

## Readiness summary
| Use | Ready? |
|---|---|
| Basic testing | ✅ yes |
| Internal use | ✅ yes |
| External beta testing | ✅ yes (with documented caveats) |
| Public GitHub release | 🟡 after CITATION/registry/tag items above |
| Manuscript / thesis methods | 🟡 methods text yes; quantitative host-fraction claims need control benchmarking |
