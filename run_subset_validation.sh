#!/usr/bin/env bash
# HÆMA subset validation run (controls + per-zone field samples) on real data,
# broad_blast reference mode against the RefSeq mitochondrion r235 database.
# Medaka disabled for speed; validates taxonomy/DB/concordance/ecology wiring end-to-end.
set -euo pipefail
export NXF_VER=24.10.5
cd "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/pipeline"
nextflow -log "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/logs/subset_validation/nextflow.log" run . \
  -profile workstation,docker -resume \
  -w "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/pipeline/work" \
  --input        "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/metadata/subset_validation_samplesheet.csv" \
  --raw_data_dir "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/data" \
  --python_container  haema-python:0.3.0 \
  --r_container       haema-r:0.3.0 \
  --figures_container haema-figures:0.4.0 \
  --reference_mode broad_blast \
  --blast_db       "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/data/reference_db/refseq_mito_r235/blastdb/refseq_mito_r235" \
  --blast_db_mount "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/data/reference_db/refseq_mito_r235/blastdb" \
  --enable_curated_panel_check true \
  --enable_medaka false \
  --outdir  "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/results/subset_validation" \
  --log_dir "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/logs/subset_validation"
