# Limitations & scientific caveats

HÆMA is research software. This page is the honest, consolidated statement of what is and is not
yet validated, so results are not over-interpreted. It covers both **scientific defensibility**
(the caveats below) and **feature status** (implemented vs staged — see the
[table at the end](#feature-status-implemented-vs-staged)).

> **Last updated:** 2026-06-30 (literature-guided re-engineering, branch `reengineer-reference-taxonomy-host-use`)

## Implemented and exercised on data

- Metadata/control validation, marker-aware primer trimming, Q20/length filtering, and
  length-based marker splitting (cyt b / short COI / long COI). Verified on real R10.4.1 data
  (~97% of reads passed Q20; primers detected).
- Pre-consensus **mixed-template denoising** with UMAP + HDBSCAN (k-mer spectra), with a
  deterministic greedy fallback for low-read groups.
- Curated **local-first BLAST** with a conservative LCA and an optional external `nt` fallback.
- Control-aware contamination flagging and a RAMBO-style mixed-host evidence/host-call layer.
- phyloseq/decontam endpoints (formal with `haema-r`; documented fallback otherwise).
- End-to-end real-data demonstration (a Ghana R10.4.1 set of 9 sequencing runs, 51 samples and
  controls): mixed feeding detected in multiple samples with biologically plausible hosts (e.g.
  Human + Goat corroborated across two markers). Single-host positive controls recovered their
  declared host (3/3), and negative controls showed only low-level *Ovis aries* background
  (count = 3) that the contamination layer flagged. **Mixed-host recovery was, however, only 50%
  (3 of 6 declared hosts):** one Human + *Bos taurus* control recovered no host signal and another
  recovered only Human — even though *Bos taurus* is in the curated panel, i.e. a
  denoising/detection-sensitivity gap, not a panel-coverage gap. The HTML report and
  `positive_control_check.tsv` surface this explicitly; see caveat 1.

## Caveats you must respect when interpreting results

1. **Denoising thresholds are not yet benchmarked.** `mixed_denoise_min_cluster_size` (default 20)
   and `mixed_denoise_min_cluster_fraction` (default 0.05) trade minority-host **sensitivity**
   for tractability and noise rejection. Low-abundance second hosts (< ~5% of marker reads) can be
   missed; very permissive values over-split noisy reads. These need calibration against
   **mixed-host positive controls** before quantitative host-fraction claims.
   *Safeguards now implemented:* the pipeline records `host_fractions_benchmarked: false` in the
   run manifest and the report; the report explicitly states host fractions are evidence summaries.
   Both *single-host* and *lab-prepared mixed-host* controls are now checked automatically against
   their declared `expected_host_scientific_name` (`positive_control_check.tsv`), reporting per-host
   recovery, missing, and unexpected hosts plus a `mixed_host_recovery_rate`. The lab-prepared
   `mf`-noted mixed-host controls (RUN09 in the demonstration set) are constructed mixtures and can
   be used directly as the calibration set once their known composition is declared — see
   [`benchmarking.md`](benchmarking.md). This benchmarks
   **detection/recovery**; quantitative host-*fraction* accuracy still needs controls prepared at
   **known mixing ratios**.
2. **UMAP/HDBSCAN clusters reads, not error-corrected haplotypes.** Clustering raw ONT reads can
   fragment a single host across clusters; downstream identity-consensus and taxonomy collapse
   these, but cluster *counts* are not host counts. Trust `host_call_table.tsv`, not raw cluster
   numbers.
3. **Taxonomy breadth is bounded by the curated panel.** The bundled panel is 20 vertebrate
   mitogenomes (common Ghanaian/peridomestic hosts). Hosts absent from the panel will be
   unassigned or mis-assigned to the nearest relative unless you supply an `nt` fallback
   (`--fallback_blast_db`). The panel still needs publication-scale expansion and versioning.
   The `nt`-fallback and taxdump-backed-LCA paths are implemented and unit-tested but were **not**
   exercised in the demonstration run (which used `--taxonomy_strategy curated_only` with the
   conservative top-hit LCA and no `--taxdump_dir`).
4. **Taxids come from the curated sidecar (not embedded in the BLAST DB) — but taxid-LCA still
   works.** The panel's descriptive deflines (`Species|Accession|...`) exceed BLAST's 50-char
   `-parse_seqids` local-id limit, so the DB is built without `-parse_seqids`/`-taxid_map`.
   Assignment joins the BLAST `sseqid` to the curated taxonomy sidecar, which **backfills the
   taxids**, so `--taxonomy_assignment_method taxid_lca --taxdump_dir <dir>` performs genuine
   taxdump-backed lowest-common-ancestor assignment against the curated panel (verified by
   `tests/test_taxid_assignment.py`; the `nt` fallback supplies its own native `staxids`).
   Re-keying the panel to accession seqids is therefore **not required** for taxid-LCA.
   **Default-path caveat:** with no `--taxdump_dir` (the default), single-taxon hits still get their
   exact sidecar taxid, but *ambiguous multi-taxon* hits are collapsed by reference-defline **genus
   string**, not by taxid lineage escalation — the pipeline prints a `WARNING` to this effect at
   runtime. Supply `--taxdump_dir` (a pinned NCBI taxonomy snapshot) for true ancestor-based
   escalation. *Remaining:* a frozen, version-labelled DB build and a bundled `taxdump` snapshot.
5. **Single basecalling assumption.** Input is already-basecalled MinKNOW FASTQ (Dorado SUP,
   R10.4.1). POD5/Dorado re-basecalling and raw-signal demultiplexing are not implemented.
6. **Depth.** Performance/feasibility examples in the docs use depth-capped subsets on a laptop;
   full-depth runs need more RAM/time or HPC. numba JIT recompiles per task (per-work-dir cache),
   which adds fixed overhead to each denoising task.

## BOLD database gap (added 2026-06-30)

The remote fallback queries NCBI `nt`, not the Barcode of Life Data System (BOLD). For COI markers
(`co1_short`, `co1_long`), BOLD contains ~9.87M curated COI barcodes (Kipp et al. 2023; Reeves
et al. 2018) and may resolve identities that NCBI `nt` cannot — particularly for avian and
reptilian hosts, which are common in some settings (Santos et al. 2019; Ogola et al. 2017).

**Current state:** COI queries fall through to NCBI `nt` only. BOLD is not queried.

**Staged:** Mode D (BOLD-aware COI) — user-supplied BOLD-derived COI FASTA database support is
planned as a reproducible local-first option. Live BOLD API integration is lower priority because
it introduces non-reproducible runtime dependence on an external service.

**Workaround:** Users may download BOLD COI sequences (or other vertebrate mitochondrial data) and
supply them via `--blast_db` in `nt_only` mode, or combine them with the curated panel into a
custom reference FASTA.

## NUMT risk (added 2026-06-30)

Nuclear mitochondrial pseudogenes (NUMTs) are non-functional mitochondrial DNA copies integrated
into the nuclear genome. They can be co-amplified alongside genuine mitochondrial DNA, producing
misleading host identifications.

**Current state:** No NUMT detection or filtering is implemented. NUMT risk is not flagged in
outputs.

**Risk by marker:**

| Marker        | Target size | NUMT risk          | Rationale                                                            |
| ------------- | ----------- | ------------------ | -------------------------------------------------------------------- |
| `co1_short` | 234 bp      | **Moderate** | Borderline above the <200 bp risk threshold (Reeves et al. 2018)     |
| `co1_long`  | 359 bp      | **Low**      | Longer amplicon reduces NUMT co-amplification probability            |
| `cyt_b`     | 359 bp      | **Low**      | Mitochondrial multicopy advantage confirmed (Hadj-Henni et al. 2015) |

**Safeguard:** NUMT risk will be reported as a per-marker caution flag in `run_manifest.json` and
the taxonomy output. This is a reporting flag, not a hard filter — NUMT removal requires database-
confirmed NUMT sequences, which are not currently available.

**Interpretation:** A NUMT risk flag does **not** prove NUMT contamination. It indicates that the
marker's amplicon length falls in a range where NUMT co-amplification is theoretically possible.
Discordance between CytB and COI host calls for the same sample may warrant NUMT investigation.

## Detection window limitation

Following Kent & Norris (2005), host DNA is reliably detectable within 24–30 hours post-feeding.
Heavily digested blood meals (equivalent to Sella scale 5–6; Santos et al. 2019) may produce
false-negative identifications even when the host is present in the curated panel.

**Current state:** Digestion degree is not tracked in the samplesheet schema. Identification rates
should be interpreted cautiously when per-stratum median read depths differ substantially.

**Staged:** Optional `digestion_class` metadata column (Recommendation 1) — not yet implemented.

## Reference-dependent sensitivity

Following Channumsin (2021) and Kipp (2023), identifications are limited by curated panel and
NCBI `nt` coverage. Taxa absent from both databases will return `no_confident_blast_hit`.

**Current state:** The curated panel contains 20 vertebrate mitogenomes covering common
Ghanaian/peridomestic hosts. The NCBI `nt` fallback extends coverage but is not BOLD-aware.

**Recommendation:** Expand and version the curated panel with checksums, accessions, and retrieval
dates. Consider adding marker-specific extracted references (e.g., COI-only, CytB-only subsets).

## Per-marker identity threshold context (added 2026-06-30)

The pipeline uses a **global** 97% minimum identity threshold for all markers. This is a
conservative choice justified by the literature:

| Marker     | BOLD/standard species-level threshold                   | Pipeline global threshold | Assessment                                                                  |
| ---------- | ------------------------------------------------------- | ------------------------- | --------------------------------------------------------------------------- |
| COI (BOLD) | ≥98% (Reeves et al. 2018; Santos et al. 2019)          | 97%                       | Conservative; high-confidence calls (pident ≥98%) can be filtered post-hoc |
| CytB       | ≥97–98% (Townzen et al. 2008; Hadj-Henni et al. 2015) | 97%                       | Appropriate                                                                 |
| 16S rRNA   | ≥90% (Logue et al. 2016)                               | 97%                       | More conservative than standard                                             |

**Recommendation:** Document per-marker interpretation in `docs/methods.md` and in the
`run_manifest.json` output. Do not change the global threshold — the 97% floor protects against
overcalling, and the `pident` field in `bloodmeal_master_endpoint.tsv` enables post-hoc filtering.

## Length-fallback marker assignment risk (added 2026-06-30)

When primer evidence is absent, reads are assigned to a marker by amplicon length window
(`allow_length_fallback = true`, the default). These reads carry a **higher marker misassignment
risk** than reads with explicit primer evidence.

**Safeguard:** The `assignment_method` field in the master endpoint table records
`assignment_method = "length_fallback"` for these reads, enabling downstream filtering.

**Recommendation:** Document this risk in `docs/methods.md`. Consider making length-fallback
opt-in (`allow_length_fallback = false`) for production runs where marker accuracy is critical.

## RAMBO threshold comparison (added 2026-06-30)

The RAMBO evidence model uses ≥3 reads AND ≥1% fraction per host per marker. This is **more
sensitive** than the >10% read fraction threshold used by Logue et al. (2016) for minor host
detection.

**Implication:** The pipeline may report mixed meals more frequently than gel-based or HRM methods.
This sensitivity is appropriate for a discovery screen but means mixed-meal rates should not be
directly compared to studies using stricter thresholds.

**Safeguard:** The RAMBO output `summary.tsv` and the run manifest record
`host_fractions_benchmarked: false`. Stricter thresholds can be applied post-hoc by filtering the
`host_fraction` column in `rambo_host_call_table.tsv`.

## Host-use framing (added 2026-06-30)

> Blood-meal host identifications describe **realised vertebrate blood-meal host use** among
> blood-fed mosquitoes. In the absence of host availability (abundance) data, **no inference
> about host preference, host selection, or diversion from host availability is warranted**
> (Gueye et al. 2023; Altahir et al. 2022).

All figures, tables, and reports must use "host use" or "blood-meal composition" language.
Labels such as "host preference profile" are **prohibited** without host-availability data.

Read fractions are **supporting evidence only** and do **not** represent proportional ingested
blood volumes (Logue et al. 2016).

## Staged / external (not bundled)

- Pattern-aware demultiplexers (Barbell, Deepbinner) — runnable via the command-template wrapper if
  you supply the tool, container, and command; not validated here.
- Medaka polishing — implemented and model-checked at runtime, but mixed-host polishing has not
  been validated against positive controls.
- Production custom images (`haema-python`, `haema-r`) must be built or pulled from a registry.
- **BOLD-aware COI mode (Mode D)** — user-supplied BOLD-derived COI FASTA support is planned.
  Live BOLD API integration is lower priority due to reproducibility concerns.
- **Multi-marker concordance analysis** — planned (Recommendation 10). Not yet implemented.
- **Host ecology statistical comparisons** — planned (Recommendation 9). Not yet implemented.
- **NUMT risk reporting** — planned (Recommendation 8). Not yet implemented.
- **Digestion class metadata** — planned (Recommendation 1). Not yet implemented.

## What would close these gaps (next scientific steps)

- Sequence and analyse **mixed-host positive controls** to calibrate denoising thresholds and
  confirm Medaka does not merge hosts.
- Expand and **version** the curated panel (checksums, accessions, retrieval dates); add
  marker-specific extracted references.
- Add a pinned NCBI **taxdump** snapshot and migrate panel seqids to accessions to enable native
  taxid LCA.
- Benchmark host-fraction estimates against known mixtures before reporting quantitative
  zoophagic indices.
- Implement **BOLD-aware COI mode** (Mode D) for expanded avian/reptilian host coverage.
- Implement **multi-marker concordance analysis** to identify potential NUMT or reference gaps.
- Implement **NUMT risk reporting** per marker.
- Implement **host ecology statistical comparisons** (Fisher's exact + Holm correction).

## Feature status (implemented vs staged)

A concise map of what the workflow does today versus what is intentionally staged for later. "Off by
default" features are real and tested but gated behind a flag/profile because they need extra inputs
or resources.

| Area                                | Implemented now                                                                                                                                                                                                                                                                                               | Staged / external                                                                                                                                                                                                                                                                                                                                             |
| ----------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Input**                     | Already-basecalled MinKNOW`fastq_pass/barcodeXX` (Dorado SUP, R10.4.1) — the standard and only ingest path                                                                                                                                                                                                 | POD5/Dorado re-basecalling and raw-signal recovery (not an accepted input type; HÆMA ingests already-basecalled FASTQ)                                                                                                                                                                                                                                       |
| **Demultiplexing**            | Trusted pre-demuxed MinKNOW folders (default);`header_tag` pooled-FASTQ demux; external-tool command-template wrapper                                                                                                                                                                                       | Barbell/Deepbinner are**not bundled** — runnable only via `--advanced_demux_command_template` with a user-supplied tool/container                                                                                                                                                                                                                    |
| **Mixed-template denoising**  | UMAP/HDBSCAN k-mer clustering with a deterministic greedy fallback; cluster FASTQs, membership tables, summaries; RAMBO-style host-evidence/host-call layer                                                                                                                                                   | Probabilistic mixed-template modelling; calibration of`min_cluster_size`/`min_cluster_fraction` against lab-prepared mixed controls (framework shipped, wet-lab controls needed)                                                                                                                                                                          |
| **Consensus / polishing**     | Greedy cluster consensus and exact dereplication; Medaka polishing**on by default** (needs the pinned ONT Medaka image; model verified at runtime by `MEDAKA_MODEL_PREFLIGHT`; the `test` profile disables it)                                                                                      | Validation that Medaka polishing does not merge co-eluting hosts (needs mixed-host controls)                                                                                                                                                                                                                                                                  |
| **Taxonomy / LCA**            | Curated-panel BLASTn; optional local`nt` fallback; taxids backfilled from the curated sidecar; conservative top-hit and exact-taxid LCA without a taxdump, and true taxdump-backed LCA with `--taxdump_dir`                                                                                               | Panel expansion/versioning, marker-specific extracted references, and a bundled NCBI taxdump snapshot. Note:`makeblastdb` is run **without** `-parse_seqids`/`-taxid_map` (the panel's descriptive deflines exceed BLAST's 50-char local-id limit); the taxid map is emitted as a provenance artifact only and taxids are joined from the sidecar |
| **Contamination & R ecology** | Negative-control background thresholds; decontam (formal with`haema-r`, documented fallback otherwise); phyloseq + decontaminated RDS/TSV endpoints                                                                                                                                                         | Quantitative host-fraction validation (needs controls at known mixing ratios)                                                                                                                                                                                                                                                                                 |
| **Ecological indices**        | Human Blood Index, animal/zoophily index, feeding-type partition, mixed-feeding rate, host-specific blood indices, and host diversity (richness/Shannon/Simpson/Pielou) with Wilson 95% CIs, stratified by ecological zone, vector sibling species, and collection period/season (`ecological_indices.tsv`) | **Planned:** Fisher's exact + Holm pairwise comparisons between strata (Recommendation 9). **Planned:** Kruskal-Wallis for host richness (Recommendation K). **Planned:** Multi-marker concordance analysis (Recommendation 10).                                                                                                            |
| **Reporting**                 | Custom HÆMA HTML report (primary); optional MultiQC custom-content summary tables; ten publication figures (`07_figures/`, PDF/SVG/PNG) built from the endpoint tables; machine-readable `run_manifest.json`                                                                                             | **Planned:** Concordance heatmap figure (Recommendation 12). **Planned:** Database-source contribution plot. **Planned:** NUMT risk summary plot. **Planned:** HBI/ABI forest plot with corrected denominators.                                                                                                                       |
