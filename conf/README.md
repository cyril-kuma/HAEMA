# Configuration

Profiles currently live in `nextflow.config`:

- `test` — tiny bundled demo data for a zero-config first run
- `local` — local workstation via Docker (laptop-safe resource ceiling)
- `docker` / `singularity` / `apptainer` — select the container engine
- `host` — no container engine (host-installed tools; debugging only, not reproducible)
- `slurm` — HPC via SLURM + Singularity
- `gpu` — optional NVIDIA GPU acceleration for Medaka (CPU fallback otherwise)
- `production` — strict metadata + custom images + full feature set

Institution-specific SLURM queue, account, cache, and storage settings can be added here as separate `*.config` files and layered with `-c`.
