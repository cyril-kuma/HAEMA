# Reference Panel Review

## Decision

Use the curated Ghana/West African vertebrate panel as the primary database, with local NCBI `nt` as fallback.

This gives the pipeline a reproducible and inspectable first-pass taxonomy layer while still allowing unexpected hosts to be recovered. Querying `nt` as the primary database would maximize breadth, but it is slower, harder to freeze for manuscript reproducibility, and more likely to assign reads to nearest available non-local relatives when West African taxa are underrepresented.

## Current Panel Contents

`references/vertebrate_dna_ref_panel.fasta` currently contains:

- 20 records
- 332,209 total bases
- complete mitochondrial genomes, 16,215-17,229 bp each
- domestic or expected hosts including human, cattle, zebu cattle, goat, sheep, pig, dog, cat, chicken, duck, guinea fowl, donkey, horse, mouse, rat, and house sparrow
- several African/sylvatic-relevant records including Nile grass rat, Natal multimammate rat, common warthog, and guinea fowl

Minor sequence-quality notes:

- `Homo_sapiens|NC_012920.1` contains one `N`
- `Passer_domesticus|NC_025611.1` contains one `N` and one `Y`

These are not blockers for BLAST, but they should be recorded in database provenance.

## Improvements Needed

Before publication-scale interpretation, improve the panel by adding:

1. A taxonomy sidecar with `seqid,taxid,scientific_name,rank,common_name,accession,source,retrieval_date`. **(Done — present and field-validated.)**
2. Accession provenance and checksums for every sequence. **(Partly done — file-level SHA256 in `assets/references/CHECKSUMS.sha256`; per-sequence `sequence_md5` in the sidecar.)**
3. Marker-specific extracted `cyt_b`, `co1_short`, and `co1_long` FASTAs. *(Remaining — `bin/build_reference_panel.py` is the starting tool.)*
4. Multiple representative haplotypes for common domestic and peridomestic hosts. *(Remaining — needs curated sequence retrieval.)*
5. Additional Ghana/West African sylvatic mammals, birds, reptiles, and other plausible blood-meal hosts. *(Remaining — needs curated sequence retrieval.)*
6. A frozen BLAST database build with version labels and checksums. *(Partly done — checksums present; an explicit DB version label remains.)*

### Repository-side tooling now available
- `assets/references/CHECKSUMS.sha256` — file integrity for the bundled panel, sidecar, and target manifest.
- `bin/verify_reference_assets.py` — verifies checksums and required sidecar fields (`seqid, taxid, scientific_name, rank, source_accession`); run by CI and `tests/validate_release.sh`.
- `bin/build_reference_panel.py` — filters a source FASTA + sidecar by the target manifest, excluding placeholder/checksum-incomplete targets, and emits a `needs-accession` report.

**Do not fabricate accessions, taxids, retrieval dates, or checksums.** Real NCBI retrieval is the
required next step for panel expansion; the tooling above validates whatever real records are added.

## Current Pipeline Behavior

The pipeline supports:

- `--taxonomy_strategy curated_then_fallback`
- `--reference_fasta`
- `--curated_reference_metadata`
- `--require_curated_taxids`
- `--fallback_blast_db <path/to/blast_db>/nt/nt`
- `--blast_db_mount <path/to/blast_db>`
- endpoint columns `blast_source`, `fallback_used`, `primary_assignment_status`, `assignment_method`, `lca_taxid`, and `top_staxids`

If curated metadata are provided, `BUILD_CURATED_TAXID_MAP` validates them and `MAKE_BLAST_DB` uses `makeblastdb -taxid_map` when possible. Without an NCBI taxdump, the pipeline can treat one exact top taxid as an exact LCA. With `--taxdump_dir`, it can perform ancestor-based LCA across multiple top taxids.

## Recommended Policy

Use curated-panel assignments as the primary thesis/manuscript layer once the panel has complete taxid metadata and marker-specific references. Use `nt` fallback calls as exploratory unless they agree with the curated panel, controls, expected biogeography, and read-level evidence.
