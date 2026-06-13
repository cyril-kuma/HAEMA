# Container Strategy

The HÆMA pipeline uses a **hybrid, digest-pinned** container strategy.

- **Public, digest-pinned images** are used wherever a tool is standard and available as a
  maintained community image (BLAST, MultiQC, Medaka), and for the stdlib-only Python and
  fallback-R steps. These run on a clean Docker machine with nothing to build.
- **Custom project images** are used only where the public base image cannot satisfy a
  scientific dependency stack: UMAP/HDBSCAN denoising (Python) and formal phyloseq/decontam
  (R/Bioconductor). These are required only by the `production` profile.

Every image is pinned to an immutable reference (an `@sha256:` digest, or for `ontresearch/medaka`
an immutable upstream `sha…` build tag) so that the same bytes are pulled on every machine and at
every future date. `:latest` and floating minor tags are not used anywhere.

> **On version numbers:** the custom images (`haema-python:0.3.0`, `haema-r:0.3.0`) are versioned
> **independently** of the pipeline (`0.2.0`). The image tag tracks changes to that image's
> dependency stack (it was bumped 0.2.0→0.3.0 when the bases were digest-pinned and `procps`/repo
> hardening were added), not the pipeline release. For publication, push the images to a registry
> and pin them by `@sha256:` digest, after which the tag is only a human label.

## Current Images

| Scope | Image | Type | Why |
| --- | --- | --- | --- |
| Default Python modules | `python:3.11@sha256:0c3b4dde…` | public, digest-pinned | Runs the stdlib-only Python scripts and the documented greedy denoising fallback. The full (not -slim) image is used because Nextflow needs `ps`/procps in the container to collect task metrics. |
| Production Python modules | `haema-python:0.3.0` | custom (slim + procps) | Adds pinned `biopython`, `pandas`, `numpy`, `scipy`, `scikit-learn`, `umap-learn`, `hdbscan` for true UMAP/HDBSCAN mixed-template denoising. Slim base + explicit `procps`. |
| BLAST | `quay.io/biocontainers/blast:2.16.0--hc155240_3@sha256:44d9d07f…` | public, digest-pinned | Standard BioContainer; works under Docker and Apptainer/Singularity. |
| Advanced demux wrapper | `python:3.11-slim-bookworm@sha256:e2d3af73…` | public, digest-pinned | The built-in `header_tag` demux is a stdlib Python wrapper. External Barbell/Deepbinner commands must provide their own validated runtime. |
| Default R outputs | `rocker/r-ver:4.4.2@sha256:df267491…` | public, digest-pinned | Sufficient for fallback RDS/table outputs when `--strict_bioconductor false` (base R only). |
| Production R outputs | `haema-r:0.3.0` | custom | Adds `phyloseq`, `decontam`, `ape` (Bioconductor 3.20, version-pinned) for formal ecological endpoints and decontamination. |
| Medaka (default and production) | `ontresearch/medaka:shaf3943918…` | public, immutable sha-tag | Official ONT image pinned to an immutable build. The required SUP/R10 model is verified at runtime by `MEDAKA_MODEL_PREFLIGHT`. |
| MultiQC | `multiqc/multiqc:v1.25@sha256:ad08efae…` | public, digest-pinned | Standard reporting image. |

## Why a hybrid strategy (and not all-public or all-custom)

- **All-public is not possible for the science.** No single maintained public image bundles
  `umap-learn` + `hdbscan` + `biopython` for the denoising step, or `phyloseq` + `decontam`
  for the ecology step. Splitting those across per-tool BioContainers would also break the
  shared `bin/fastq_utils.py` helper the Python steps import. A custom image is the correct,
  reproducible choice for exactly those two stacks.
- **All-custom is wasteful and less portable.** BLAST, MultiQC, and Medaka have excellent,
  maintained upstream images. Re-wrapping them adds local-build friction and a maintenance
  burden for no scientific gain. Medaka in particular only needed a *model check*, which the
  pipeline already performs at runtime — so production uses the official pinned image directly.

The result is two custom images (`haema-python`, `haema-r`) instead of three, and a fully
public, build-free path for development, testing, and any non-`production` run.

## Building the custom images (production only)

From the `pipeline/` directory (or use the `Makefile` targets `make images` / `make push`):

```bash
docker build -t haema-python:0.3.0 -f containers/haema-python/Dockerfile .
docker build -t haema-r:0.3.0     -f containers/haema-r/Dockerfile .

# Push to your registry and get the digests to pin (replace REGISTRY with your namespace):
make push REGISTRY=ghcr.io/USER        # then set --python_container / --r_container to @sha256:...
```

Both Dockerfiles pin their base image by digest, pin every package version, and fail the build
if the scientific stack cannot be imported, so a successful build is reproducible.

The optional `haema-medaka:0.3.0` image is **not required** — it exists only to additionally
assert the SUP/R10 model at *build* time (the pipeline already asserts it at *run* time):

```bash
docker build -t haema-medaka:0.3.0 -f containers/haema-medaka/Dockerfile .   # optional
```

## Long-term reproducibility and HPC

- For publication/HPC release, push `haema-python:0.3.0` and `haema-r:0.3.0` to a registry and
  set `--python_container` / `--r_container` to the resulting immutable `@sha256:` digests.
- On HPC, run with `-profile apptainer` or `-profile singularity` and set a shared
  `NXF_SINGULARITY_CACHEDIR` (the profiles read it) so converted SIFs are pulled once and reused.
- All public images already carry digests in `nextflow.config`, so no extra pinning is needed
  for them.

## Known Risks / Limitations

- `haema-python:0.3.0` and `haema-r:0.3.0` are local tags until built or pushed to a registry;
  a fresh machine running `-profile production` must build them first (or be pointed at registry
  digests). The development/test path needs no build.
- The pinned Medaka image is multi-GB; pre-pull it (or pre-stage the SIF) before offline runs.
- Barbell/Deepbinner remain external. Their commands and images must be pinned by the user when
  real pooled/raw inputs are introduced via `--advanced_demux_command_template`.
- Image digests were resolved on 2026-06-11. They are immutable; re-resolve only if you
  deliberately upgrade a tool version.
