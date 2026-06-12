<!-- Thanks for contributing to HÆMA! Please complete the checklist below. -->

## Summary

<!-- What does this PR change, and why? -->

## Type of change

- [ ] Bug fix
- [ ] New feature
- [ ] Documentation
- [ ] Scientific logic (denoising / taxonomy / contamination / host-calling)

## Checklist

- [ ] `nextflow config .` parses for `test`, `local`, and `production` profiles
- [ ] `nextflow run . -profile test,docker --outdir results/test` completes (exit 0)
- [ ] New/changed processes have a `stub:` block and a resource `label`
- [ ] New parameters are in `nextflow_schema.json` (description + default)
- [ ] Containers pinned by digest / immutable tag (no `:latest`)
- [ ] No hardcoded absolute paths; no large or private data committed
- [ ] `CHANGELOG.md` updated under `## [Unreleased]`
- [ ] Docs updated if behaviour / inputs / outputs changed

## Scientific note (if applicable)

<!-- Rationale and any evidence on test data / positive controls. -->
