#!/usr/bin/env bash
# HÆMA FULL real-data run (machine-specific convenience wrapper).
#
# This is only a thin launcher: production DEFAULTS live in nextflow.config + the `workstation`
# profile (containers, Medaka on, curated pre-check, denoising, RAMBO thresholds). This wrapper
# sets only MACHINE-SPECIFIC values (input, data dir, the broad BLAST database path, output).
#   * Quick start (no broad DB; curated panel + bundled reference, works out of the box):
#       nextflow run . -profile workstation,docker --input SHEET.csv --raw_data_dir DIR --outdir results
#   * Broad run (this script): adds --reference_mode broad_blast + a local RefSeq-mito database.
# Re-run verbatim to -resume after any interruption (work dir preserved). Effective parameters are
# identical to the previous explicit form, so -resume caching is preserved.
set -euo pipefail
export NXF_VER=24.10.5          # 24.10 LTS: the host's newer strict parser breaks this pipeline

PROJ="/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding"
DB="$PROJ/data/reference_db/refseq_mito_r235/blastdb"
cd "$PROJ/pipeline"

nextflow -log "$PROJ/logs/full_run_v2/nextflow.log" run . \
  -profile workstation,docker -resume \
  -w "$PROJ/pipeline/work" \
  --input        "$PROJ/metadata/full_run_samplesheet.csv" \
  --raw_data_dir "$PROJ/data" \
  --reference_mode broad_blast \
  --blast_db       "$DB/refseq_mito_r235" \
  --blast_db_mount "$DB" \
  --outdir  "$PROJ/results/full_run_v2" \
  --log_dir "$PROJ/logs/full_run_v2"
# Containers, --enable_curated_panel_check, --enable_medaka and RAMBO thresholds now come from the
# workstation profile / config defaults — no longer duplicated here.
