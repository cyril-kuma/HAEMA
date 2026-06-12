# Contributing to HÆMA

Thanks for your interest in improving HÆMA. This is a research pipeline under active
development; contributions, bug reports, and scientific feedback are all welcome.

## Ways to contribute

- **Report a bug** — open an issue using the bug-report template. Include the exact command,
  the Nextflow version (`nextflow info`), the profile used, and the relevant lines from the
  `.command.err`/`.command.log` of the failing task work directory.
- **Request a feature** — open an issue using the feature-request template.
- **Submit a change** — fork, branch, and open a pull request against `main`.

## Development setup

```bash
git clone <your-fork-url> haema && cd haema
# Requirements: Nextflow >= 24.10.0, Java 17-24, and Docker (or Singularity/Apptainer).
nextflow run . -profile test,docker --outdir results/test    # smoke test (~1 min)
```

## Pull-request checklist

- [ ] `nextflow config .` parses for the `test`, `local`, and `production` profiles.
- [ ] `nextflow run . -profile test,docker --outdir results/test` completes (exit 0).
- [ ] New/changed processes have a `stub:` block and a resource `label`.
- [ ] New parameters are added to `nextflow_schema.json` with a description and default.
- [ ] Container images are pinned by digest or an immutable tag — never `:latest`.
- [ ] No hardcoded absolute paths; no large or private data committed.
- [ ] `CHANGELOG.md` updated under `## [Unreleased]`.
- [ ] Docs updated if behaviour, inputs, or outputs changed.

## Coding conventions

- **Nextflow**: DSL2; lowercase channel factories (`channel.of`), explicit closure parameters,
  `def` for local variables, named `emit:` outputs, resource `label`s (`process_low/medium/high`).
- **Python** (`bin/`): standard-library only for the default-image steps; scientific
  dependencies (numpy/sklearn/umap/hdbscan) belong in `haema-python` and must degrade gracefully.
- **Style**: keep new code consistent with the surrounding file (naming, comment density, idiom).

## Scientific changes

Changes that affect denoising, taxonomy, contamination, or host-calling logic should include a
short rationale in the PR and, where possible, evidence on the bundled test data or a documented
positive control. Do not overclaim: distinguish implemented behaviour from planned enhancements.

## Code of conduct

By participating you agree to abide by the [Code of Conduct](CODE_OF_CONDUCT.md).
