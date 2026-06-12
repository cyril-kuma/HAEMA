# Assets

Small, versioned files shipped with the pipeline so a fresh clone is self-contained.

| Path | What it is |
|---|---|
| `primers.csv` | Tri-marker assay primers (`cyt_b`, `co1_short`, `co1_long`). Default for `--primers`. |
| `references/vertebrate_dna_ref_panel.fasta` | Curated Ghana/West-African vertebrate mitogenome panel (~20 records). Default for `--reference_fasta`. |
| `references/vertebrate_dna_ref_panel.taxonomy.tsv` | Taxonomy sidecar (`seqid,taxid,scientific_name,rank,...`). Default for `--curated_reference_metadata`. |
| `references/ghana_west_africa_reference_targets.tsv` | Accession-driven target manifest used by production preflight. |
| `references/CHECKSUMS.sha256` | SHA256 integrity manifest for the reference files (checked by `bin/verify_reference_assets.py`). |
| `samplesheets/example_samplesheet.csv` | MIEM/MIMARKS samplesheet template — copy and edit for your run. |
| `samplesheets/mixed_host_control_example.csv` | Illustrative single + lab-prepared mixed-host control rows (see `docs/benchmarking.md`). |
| `empty_curated_reference_metadata.tsv` | Fallback used when no taxonomy sidecar is supplied. |
| `NO_FILE` | Placeholder for optional process inputs. |
| `test_data/` | Tiny synthetic fixtures used by `-profile test`. |

These assets are the **single source of truth** for the pipeline's bundled reference data — there
are no duplicate copies elsewhere in the tree. Override any default by passing the corresponding
`--` parameter. Large databases (e.g. NCBI `nt`) are **not** vendored — supply them with
`--fallback_blast_db` / `--blast_db_mount`. The curated panel still needs publication-scale
expansion and checksum versioning (see `docs/limitations.md`).
