#!/usr/bin/env bash
# HÆMA FULL real-data run (all field samples + controls with data on disk),
# broad_blast reference mode: curated panel pre-check + RefSeq mitochondrion r235
# broad database (catches unexpected vertebrate hosts). Medaka polishing ON.
# Re-run verbatim to -resume after any interruption (work dir preserved).
set -euo pipefail
export NXF_VER=24.10.5          # 24.10 LTS: the host's newer strict parser breaks this pipeline
cd "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/pipeline"
nextflow -log "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/logs/full_run_v2/nextflow.log" run . \
  -profile workstation,docker -resume \
  -w "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/pipeline/work" \
  --input        "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/metadata/full_run_samplesheet.csv" \
  --raw_data_dir "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/data" \
  --python_container  haema-python:0.3.0 \
  --r_container       haema-r:0.3.0 \
  --figures_container haema-figures:0.4.0 \
  --reference_mode broad_blast \
  --blast_db       "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/data/reference_db/refseq_mito_r235/blastdb/refseq_mito_r235" \
  --blast_db_mount "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/data/reference_db/refseq_mito_r235/blastdb" \
  --enable_curated_panel_check true \
  --enable_medaka true \
  --outdir  "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/results/full_run_v2" \
  --log_dir "/home/waccbip/Documents/bioinformatics/01_Bloodmeal_Host_Metabarcoding/logs/full_run_v2"
