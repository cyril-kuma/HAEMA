# Reproducibility statement

HÆMA is designed so that the **same pipeline revision + the same pinned containers + `-resume`**
produce identical results on any machine with Nextflow and a container engine.

## What is pinned
- **Nextflow engine** — pin with `export NXF_VER=24.10.0` (the version this release was tested on).
- **Pipeline revision** — clone a tagged release, or `nextflow run <repo> -r <tag>`.
- **Plugins** — `nf-schema@2.2.0` is version-pinned in `nextflow.config` (parameter validation +
  `--help`). It is downloaded once and cached; for offline/air-gapped runs pre-stage it with
  `nextflow plugin install nf-schema@2.2.0` or point `NXF_PLUGINS_DIR` at a shared cache.
- **Containers** — every public image is referenced by an immutable `@sha256:` digest (and the
  Medaka image by an immutable upstream sha-tag) in `nextflow.config`; nothing uses `:latest`.
  The two custom images (`haema-python`, `haema-r`) pin their base by digest and every package
  version, and self-test their imports at build time. For publication, push them to a registry and
  set `--python_container` / `--r_container` to the resulting digests.
- **Reference data** — the curated panel, taxonomy sidecar, and target manifest are vendored under
  `assets/references/` and versioned with the repo.
- **Parameters & provenance** — every run writes `05_endpoint_files/run_manifest.json` recording
  parameters, container images, profile, command line, and session ID. Nextflow also emits an
  execution trace, timeline, report, and DAG under `pipeline_info/`.

## Environments
- **Local / Docker:** `-profile local` (or `docker`). Digest-pinned images pulled once and cached.
- **HPC / Singularity-Apptainer:** `-profile slurm` (or `singularity` / `apptainer`). Set a shared
  `NXF_SINGULARITY_CACHEDIR` so SIFs are converted once and reused; the profiles read it.
- **GPU:** `-profile gpu` adds `--gpus all` for Medaka only; requires the NVIDIA Container Toolkit.
  CPU is the automatic fallback.
- **Offline / air-gapped:** pre-pull images (or pre-build SIFs), set `NXF_OFFLINE=true`, and use a
  local registry or `NXF_SINGULARITY_CACHEDIR`.

## No machine-specific state
- Output, log, and work directories default to `launchDir`-relative paths.
- `--input`, `--raw_data_dir`, and external BLAST databases are the only paths you supply; pass
  them as absolute paths so they resolve inside containers.
- No tool is assumed to be installed on the host (use `-profile host` only for debugging).

## Known reproducibility caveats
- UMAP/HDBSCAN uses `random_state=42` but parallel numba kernels can introduce minor
  run-to-run variation at the read-cluster level; host-level calls are stable.
- The curated panel is not yet checksum-versioned end to end (see
  [`limitations.md`](limitations.md)).
