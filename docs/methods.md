# Methods (draft for manuscript / thesis)

A methods-section-ready description of the HÆMA pipeline. Fill the bracketed values from your
`run_manifest.json` and update tool versions to match the images you used. This describes the
**implemented** workflow; see [`limitations.md`](limitations.md) for caveats and staged features.

## Bioinformatic analysis

Oxford Nanopore R10.4.1 reads, basecalled in MinKNOW with the Dorado super-accuracy model
(`[basecalling_model]`) and demultiplexed into `fastq_pass/barcodeXX` folders, were analysed with
the HÆMA pipeline (v`[version]`), a containerised Nextflow DSL2 workflow
(Nextflow `[NXF_VER]`; Di Tommaso et al., 2017). Each sample was linked to its run, barcode, and
MIEM/MIMARKS-compliant ecological metadata, and to its assay control class (biological sample,
extraction blank, PCR negative, or single-host positive control).

Per sample, FASTQ chunks were merged and reads were processed against the tri-marker assay
(cytochrome *b*, short COI, long COI). Terminal primers were detected and trimmed allowing an
error rate of `[primer_max_error_rate]` (IUPAC-aware), reads below mean Phred `[min_mean_q]` or
`[min_read_length]` bp were discarded, and surviving reads were assigned to a marker by primer
evidence and amplicon length window. Reads assigned by length fallback (when primer evidence is
absent) carry a higher marker misassignment risk; these are flagged with
`assignment_method = "length_fallback"` in the master endpoint table for downstream filtering.

### Mixed-template denoising

To resolve mixed blood meals, marker reads were denoised before consensus generation: 5-mer
spectra were embedded with UMAP (McInnes et al., 2018) and clustered with HDBSCAN
(McInnes et al., 2017), retaining clusters with ≥`[mixed_denoise_min_cluster_size]` reads and
≥`[mixed_denoise_min_cluster_fraction]` of the sample/marker reads; groups below
`[mixed_denoise_min_reads_for_umap]` reads used a deterministic greedy identity-clustering
fallback. Per retained cluster, a consensus/representative sequence was generated
(`[consensus_method]`; optional Medaka polishing, Oxford Nanopore Technologies, with model
`[medaka_model]`).

**Note on denoising layering:** UMAP+HDBSCAN denoising operates on *reads* (pre-BLAST) to partition
mixed templates by sequence composition. The RAMBO evidence model (below) operates on *host-assigned
reads* (post-BLAST) to accumulate evidence per host. These are complementary, not redundant,
detection layers. The concept is supported by Channumsin et al. (2021), who confirmed multiple
host templates co-occur in individual mosquito abdomens via TOPO-TA cloning, and by Kipp et al.
(2023), who showed that sufficient read depth enables assembly of multiple distinct mitogenome
sequences — the gold standard for mixed-template detection that UMAP+HDBSCAN approximates for
lower-depth ONT data.

### Taxonomic assignment

Consensus sequences were assigned taxonomy by BLASTn (Camacho et al., 2009; v`[blast_version]`)
against a curated Ghana/West-African vertebrate mitogenome panel, accepting hits with
≥`[min_blast_identity]`% identity and ≥`[min_blast_coverage]`% query coverage and resolving ties
with a conservative lowest-common-ancestor rule (`[taxonomy_assignment_method]`); taxids were taken
from the curated taxonomy sidecar. Features unresolved by the curated panel were optionally queried
against a local NCBI `nt` database.

**Per-marker identity threshold context:** The pipeline uses a global 97% minimum identity
threshold for all markers. This is a conservative choice:

| Marker | BOLD/standard species-level threshold | Pipeline global threshold | Assessment |
|--------|--------------------------------------|--------------------------|------------|
| COI (BOLD) | ≥98% (Reeves et al. 2018; Santos et al. 2019) | 97% | Conservative; high-confidence calls (pident ≥98%) can be filtered post-hoc |
| CytB | ≥97–98% (Townzen et al. 2008; Hadj-Henni et al. 2015) | 97% | Appropriate |
| 16S rRNA | ≥90% (Logue et al. 2016) | 97% | More conservative than standard |

These per-marker thresholds are now enforced automatically as a conservative confidence guard:
a species-level assignment that would otherwise be labelled `high` is downgraded to `medium` when
its top identity falls below the marker's species threshold (`--coi_species_identity_threshold`,
`--cytb_species_identity_threshold`). This never loosens the global pass filter
(`--min_blast_identity`, default 97%); it only prevents species-level overcalling. High-confidence
assignments can additionally be filtered post-hoc using the `pident` field in
`bloodmeal_master_endpoint.tsv`.

**Lowest Common Ancestor (LCA) assignment:** Where top BLAST hits map to multiple species at equal
or near-equal bitscore (within `top_bitscore_delta` = 2.0), the pipeline collapses the assignment
to the Lowest Common Ancestor. When `--taxdump_dir` is supplied, taxid-based lineage escalation
from NCBI taxonomy is used. Without `--taxdump_dir`, single-taxon hits receive their exact sidecar
taxid, but ambiguous multi-taxon hits are collapsed by reference-defline genus string (less
reliable for polyphyletic genus names). Production runs should supply `--taxdump_dir` for
taxonomic rigour.

**Reference database modes:** The reference strategy is selected with `--reference_mode`, which is
the canonical control. The legacy `--taxonomy_strategy` (`curated_only`/`curated_then_fallback`/
`nt_only`) remains honoured and is mapped onto a reference mode when `--reference_mode` is left at
its default, so existing configurations keep working. Each mode uses at most one fallback database:

- **Mode A — `curated_panel`:** Curated vertebrate panel only. Fast and offline; limited by panel
  completeness but appropriate for expected Ghanaian/peridomestic hosts.
- **Mode B — `broad_blast`:** User-supplied broad local BLAST database (`--blast_db`, e.g. an NCBI
  nt subset, a custom vertebrate mitochondrial database, or a combined curated + public reference).
  With `--enable_curated_panel_check true` (default) the curated panel is queried first and the
  broad database resolves only the unresolved features; set it false to query the broad database
  alone.
- **Mode C — `remote_fallback`:** Curated panel first, then NCBI `nt` via `blastn -remote`
  (`--remote_blast_db`, default `core_nt`) for features the curated panel cannot resolve. Gated by
  `--enable_ncbi_remote_fallback`. Remote queries are **not reproducible by default** (the remote
  database is versionless from the client side) and should be reserved for exploratory look-ups.
- **Mode D — `bold_aware`:** Curated panel first, then a reproducible, project-local BOLD-derived
  COI FASTA database (`--bold_fasta`, optionally `--bold_taxonomy`) built into a BLAST database at
  runtime. It is effective for the COI markers (`co1_short`, `co1_long`); CytB queries against a COI
  database simply return no hits and are unaffected. Live BOLD API querying (`--bold_mode
  api_query`) is intentionally **not** implemented — a network-dependent lookup would make routine
  runs non-reproducible; supply a downloaded BOLD FASTA instead.

The reference database actually used, its SHA-256 checksum, the fallback chain, and the per-marker
identity thresholds and NUMT risk are all recorded in `run_manifest.json` under `reference_database`.

**NUMT risk:** Nuclear mitochondrial pseudogenes (NUMTs) can be co-amplified alongside genuine
mitochondrial DNA. Risk by marker:

| Marker | Target size | NUMT risk | Rationale |
|--------|-------------|-----------|-----------|
| `co1_short` | 234 bp | **Moderate** | Borderline above the <200 bp risk threshold (Reeves et al. 2018) |
| `co1_long` | 359 bp | **Low** | Longer amplicon reduces NUMT co-amplification probability |
| `cyt_b` | 359 bp | **Low** | Mitochondrial multicopy advantage confirmed (Hadj-Henni et al. 2015) |

NUMT risk is reported as a per-marker caution flag, not a hard filter. A NUMT risk flag does
**not** prove NUMT contamination; it indicates theoretical possibility. Discordance between CytB
and COI host calls for the same sample may warrant NUMT investigation.

### Host evidence model (RAMBO)

Host assemblages per sample/marker were summarised with an abundance/evidence model that preserves
multiple hosts above configurable read and fraction thresholds, classifying each as single- or
mixed-host. The RAMBO evidence model requires ≥3 supporting reads AND ≥1% of marker reads per host
per sample.

**Threshold comparison:** The 1% fraction threshold is more permissive than the >10% threshold
used by Logue et al. (2016) for minor host detection. This increased sensitivity is appropriate
for a discovery screen but means the pipeline may report mixed meals more frequently than gel-based
or HRM methods. Stricter thresholds can be applied post-hoc by filtering the `host_fraction`
column in `rambo_host_call_table.tsv`.

**Important:** Host fractions are **supporting evidence only** and do **not** represent proportional
ingested blood volumes (Logue et al. 2016). The run manifest records `host_fractions_benchmarked:
false` to explicitly flag this limitation.

### Ecological indices

From the per-mosquito host calls (hosts unioned across markers, controls excluded), standard
vector–host ecological indices were computed: the Human Blood Index (proportion of host-identified
blood meals containing human; Garrett-Jones, 1964), the complementary animal/zoophily index, a
human-only/mixed/animal-only feeding-type partition, the mixed-feeding rate (≥2 host taxa),
host-specific blood indices, and host-community diversity (richness, Shannon, Gini–Simpson, Pielou
evenness), each with Wilson score 95% confidence intervals and reported overall and stratified by
ecological zone, vector sibling species, and collection period/season. Because sampling was
opportunistic and each collection campaign coincided with a distinct site and sequencing batch, the
temporal/seasonal strata are reported descriptively and are confounded with location and batch; no
temporal trend or seasonal-effect test was performed.

**Scientific framing:** Blood-meal host identifications describe **realised vertebrate blood-meal
host use** among blood-fed mosquitoes. In the absence of host availability (abundance) data, **no
inference about host preference, host selection, or diversion from host availability is warranted**
(Gueye et al. 2023; Altahir et al. 2022). All figures, tables, and reports must use "host use" or
"blood-meal composition" language. Labels such as "host preference profile" are prohibited without
host-availability data.

**Detection window:** Following Kent & Norris (2005), host DNA is reliably detectable within
24–30 hours post-feeding. Heavily digested blood meals (Sella scale 5–6; Santos et al. 2019) may
produce false-negative identifications even when the host is present in the curated panel.

**Reference-dependent sensitivity:** Following Channumsin (2021) and Kipp (2023), identifications
are limited by curated panel and NCBI `nt` coverage. Taxa absent from both databases will return
`no_confident_blast_hit`.

### Multi-marker concordance

For samples with host calls from two or more markers, concordance analysis assesses whether host
assignments agree at species or genus level. Discordant calls between CytB and COI for the same
sample may indicate NUMT co-amplification, reference database gaps, or genuine within-sample
complexity. Concordance status values include: `full_species_agreement`, `genus_agreement`,
`discordant`, `single_marker_only`, `no_marker_signal`, `ambiguous_lca_only`.

**Status:** Planned (Recommendation 10). Not yet implemented.

### Host ecology statistical comparisons

Formal statistical comparisons between HBI strata (Fisher's exact test with Holm-Bonferroni
correction) are planned for:

- Pairwise HBI comparisons between ecological zones
- Pairwise HBI comparisons between sibling species (if sample sizes permit)
- Mixed-feeding rate comparisons between zones
- Host richness comparisons between zones (Kruskal-Wallis, if meaningful)

**Status:** Planned (Recommendation 9). All tests should be framed as exploratory unless
pre-specified as primary. Small-n warnings are applied when any stratum has n < 5.

## Reproducibility

All steps ran in version-pinned containers (digest-pinned). The pipeline, parameters, container
images, and execution traces are recorded in `run_manifest.json` and `pipeline_info/`. See
[`reproducibility.md`](reproducibility.md).

## Key references (verify and format to your target journal)
- Di Tommaso P, et al. Nextflow enables reproducible computational workflows. *Nat Biotechnol*. 2017.
- McInnes L, Healy J, Melville J. UMAP. *arXiv*:1802.03426. 2018.
- McInnes L, Healy J, Astels S. hdbscan. *J Open Source Softw*. 2017.
- Camacho C, et al. BLAST+. *BMC Bioinformatics*. 2009.
- Davis NM, et al. decontam. *Microbiome*. 2018.
- McMurdie PJ, Holmes S. phyloseq. *PLoS ONE*. 2013.
- Ewels P, et al. MultiQC. *Bioinformatics*. 2016.
- Garrett-Jones C. The human blood index of malaria vectors… *Bull World Health Organ*. 1964;30:241–261.
- Orsborne J, et al. Using the human blood index to investigate host biting plasticity… *Malar J*. 2018;17:479.

> Verify every citation and the exact tool versions against the images you ran before submission.
