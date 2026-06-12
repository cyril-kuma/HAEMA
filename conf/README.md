# Configuration

Profiles currently live in `nextflow.config`:

- `local`
- `docker`
- `singularity`
- `apptainer`
- `slurm`
- `test`
- `production`

Institution-specific SLURM queue, account, cache, and storage settings can be added here as separate `*.config` files and layered with `-c`.
