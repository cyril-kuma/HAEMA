# HÆMA Pipeline Containers

These Dockerfiles define the project-specific runtime images used by the `production` profile.

Development and test profiles use **public, digest-pinned** images by default, so anyone can run
the scaffold on a clean machine with nothing to build. Production mode switches to these custom
images only where a public base cannot provide the scientific dependency stack:

- `haema-python` — UMAP/HDBSCAN mixed-template denoising (`umap-learn`, `hdbscan`, `scikit-learn`,
  `numpy`, `scipy`, `pandas`, `biopython`).
- `haema-r` — formal `phyloseq` + `decontam` (Bioconductor 3.20) ecology endpoints.
- `haema-figures` — matplotlib/seaborn/geopandas publication figures, including the Ghana
  bioclimatic-zone map.

Medaka is **not** a custom image: the production profile pins the official upstream
`ontresearch/medaka` image (immutable sha-tag) and verifies the required model at runtime via
`MEDAKA_MODEL_PREFLIGHT`. See `../docs/CONTAINER_STRATEGY.md` for the full rationale.

All custom Dockerfiles pin their base by digest, pin every package version, install `procps` (Nextflow
needs `ps` for task metrics), and fail the build if the stack cannot be imported — so a successful
build is reproducible.

## Build (from the `pipeline/` directory)

```bash
docker build -t haema-python:0.3.0 -f containers/haema-python/Dockerfile .
docker build -t haema-r:0.3.0     -f containers/haema-r/Dockerfile .
docker build -t haema-figures:0.4.0 -f containers/haema-figures/Dockerfile .
```

## Optional hardened Medaka image

Only needed if you want the SUP/R10 model asserted at **build** time (the pipeline already asserts
it at **run** time) or baked into a frozen image for an air-gapped HPC:

```bash
docker build -t haema-medaka:0.3.0 -f containers/haema-medaka/Dockerfile .
```

The base defaults to a pinned immutable upstream sha-tag (not `:latest`). Override it explicitly:

```bash
docker build \
  --build-arg MEDAKA_BASE_IMAGE='ontresearch/medaka:<pinned-sha-tag-or-digest>' \
  -t haema-medaka:0.3.0 \
  -f containers/haema-medaka/Dockerfile .
```

Verify the model is present:

```bash
docker run --rm haema-medaka:0.3.0 medaka tools list_models | grep -F r1041_e82_400bps_sup_v4.3.0
```

## Registry / HPC

For production/HPC release, push `haema-python:0.3.0`, `haema-r:0.3.0`, and
`haema-figures:0.4.0` to a registry and set `--python_container` / `--r_container` /
`--figures_container` to the resulting immutable `@sha256:` digests. For
Apptainer/Singularity, set a shared `NXF_SINGULARITY_CACHEDIR` (the profiles read it) so the
converted SIFs are built once and reused across nodes.
