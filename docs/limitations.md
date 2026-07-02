Local AI Agent Prompt: Review, Refactor, and Release-Prepare a Nextflow Pipeline for ONT Blood-Meal Amplicon Metabarcoding

You are a senior bioinformatician, Nextflow DSL2 engineer, and molecular epidemiology reviewer. You have been asked to review, refine, and professionalise a custom Nextflow pipeline that processes Oxford Nanopore amplicon-sequencing data from blood-fed Anopheles mosquitoes.

The pipeline implements Objective 1 of the study: identifying vertebrate blood-meal host species from blood-fed female Anopheles gambiae s.l. using amplicon metabarcoding. The final pipeline should be scientifically accurate, reproducible, modular, well documented, and suitable for public release.

The skills package is installed globally. Use relevant skills whenever they are better suited than raw shell commands for inspecting files, reviewing code, running tests, modifying files, generating documentation, or checking workflow quality.

Your task is not only to review the repository, but to actively improve it. You may read files, run commands, edit files, create tests, update documentation, refactor modules, and propose or implement missing components. Work carefully and iteratively.

⸻

1. Working principles

Follow these rules throughout:

1. Be meticulous and systematic.
2. Prefer evidence from the repository over assumptions.
3. Do not pretend the existing implementation is correct.
4. Do not over-engineer unnecessarily.
5. Preserve scientifically valid existing functionality.
6. Refactor, merge, rename, or remove components where needed.
7. Keep the workflow compatible with Nextflow DSL2.
8. Keep user-facing parameters in configuration files, not hard-coded inside scripts.
9. Pin software versions wherever possible.
10. Record scientific assumptions and justify parameter choices.
11. After major edits, rerun relevant checks and tests.
12. Use -resume when rerunning Nextflow workflows to avoid redundant computation.
13. Produce a clear final report summarising findings, changes, tests, and remaining limitations.

⸻

2. Repository exploration

First, inspect the full repository.

Use skills, shell, or Python to:

1. List the complete directory tree.
2. Identify the main workflow file, usually main.nf.
3. Identify supporting files and directories, including:
   * nextflow.config
   * modules/
   * subworkflows/
   * conf/ or configs/
   * bin/
   * assets/
   * docs/
   * README.md
   * CHANGELOG.md
   * LICENSE
   * environment files
   * container definitions
   * test data
   * nf-test or other test files
4. Determine whether the pipeline uses Nextflow DSL2.
5. Determine whether the pipeline structure resembles nf-core style or a custom DSL2 layout.
6. Identify any missing standard files needed for a professional public pipeline.

Create an initial repository map with:

* File or directory
* Purpose
* Status
* Immediate concern, if any

⸻

3. Code and workflow review

Read main.nf, all module files, all subworkflow files, and relevant helper scripts.

Evaluate:

1. Whether DSL2 is enabled and used correctly.
2. Whether workflow logic is separated from configuration.
3. Whether processes are modular and reusable.
4. Whether channel construction is clear and correct.
5. Whether each process has clear input:, output:, and script: blocks.
6. Whether outputs use named emit: channels where appropriate.
7. Whether publishDir is used consistently and safely.
8. Whether parameters are passed cleanly rather than hard-coded.
9. Whether repeated logic should be moved into modules, subworkflows, or helper scripts.
10. Whether embedded bash is too complex and should be converted into scripts in bin/.
11. Whether the pipeline supports -resume.
12. Whether errors are likely to be understandable to users.
13. Whether the workflow can be run from a clean checkout using documented commands.

Where needed, refactor the pipeline to improve clarity, modularity, and maintainability.

⸻

4. Configuration review

Open and review nextflow.config and any profile-specific configuration files.

Check for:

1. A valid manifest block with pipeline name, description, author, version, homepage, and relevant metadata.
2. Sensible default parameters.
3. User-facing parameters defined in configuration rather than hard-coded inside workflow scripts.
4. Standard profiles such as:
   * standard
   * test
   * docker
   * singularity or apptainer
   * conda
   * local
   * hpc or cluster-specific profiles, where appropriate
5. Resource declarations in the process scope:
   * cpus
   * memory
   * time
   * errorStrategy
   * maxRetries
6. Per-process resource overrides using withName: or labels.
7. Fixed container versions or pinned conda package versions.
8. Clear executor settings that can be adapted for local machines and HPC.
9. Avoidance of machine-specific absolute paths.
10. Good defaults for output directories and temporary directories.

Improve the configuration where necessary.

⸻

5. Scientific objective review

The pipeline must implement Objective 1:

Characterise realised vertebrate blood-meal host use of wild-caught, blood-fed female Anopheles gambiae s.l. across Ghanaian bioclimatic zones using Oxford Nanopore amplicon metabarcoding.

Evaluate whether the pipeline correctly supports this objective.

The pipeline should support, or clearly document limitations around, the following steps:

1. Input validation from samplesheet and FASTQ files.
2. Demultiplexed ONT reads or raw/barcoded reads, depending on available data.
3. Barcode and adapter trimming where needed.
4. Primer trimming.
5. Read quality filtering.
6. Amplicon length filtering for expected marker sizes.
7. Marker-aware processing for targets such as:
   * COI_long
   * COI_short
   * Cyt b
8. Host mitochondrial DNA assignment using a curated vertebrate reference database.
9. Conservative taxonomic assignment using BLAST, minimap2, vsearch, kraken-like tools, QIIME2, or another justified method.
10. Lowest-common-ancestor or conservative assignment logic where species-level calls are uncertain.
11. Negative-control-aware filtering.
12. Cross-barcode contamination or barcode hopping assessment.
13. Multi-marker consolidation per specimen.
14. Mixed-host detection.
15. Primary host call assignment.
16. Per-specimen host-call table.
17. Summary tables for HBI, ABI, host spectrum, and mixed feeding.
18. Outputs suitable for thesis, manuscript, and downstream integration.

Critically assess whether the chosen tools and thresholds are appropriate for ONT amplicon metabarcoding. If they are weak or missing, implement improvements or document clear recommendations.

Do not treat read fractions as blood-volume estimates. Do not frame the pipeline as measuring host preference unless host-availability data are included. The correct interpretation is realised host use.

⸻

6. Input and metadata validation

Check whether the pipeline validates required inputs.

At minimum, the pipeline should validate:

1. Samplesheet presence and format.
2. Unique sample or specimen IDs.
3. Correct mapping between specimens and FASTQ files.
4. Required columns, such as:
   * specimen_id
   * barcode or sample identifier
   * sequencing run
   * marker or primer information, if applicable
   * bioclimatic zone, if used in summaries
   * sibling species, if used in stratified analysis
5. Existence of FASTQ files.
6. Existence of primer files.
7. Existence and indexing of the host reference database.
8. Output directory writability.
9. Negative and positive control labelling.
10. Consistent zone/site naming.

Implement validation scripts if missing.

⸻

7. Output and reporting review

Check whether outputs are clear, reproducible, and useful.

The pipeline should ideally produce:

1. Per-sample or per-specimen read QC summaries.
2. Per-marker read counts.
3. Primer/amplicon filtering summaries.
4. Taxonomic assignment tables.
5. LCA or conservative host-call tables.
6. Mixed-host evidence table.
7. Negative-control and contamination-filtering report.
8. Final specimen-level host-call table.
9. HBI and ABI summary tables.
10. Host spectrum tables by zone and/or sibling species.
11. Mixed-feeding summary tables.
12. MultiQC report or equivalent QC report.
13. Publication-quality plots.
14. A pipeline provenance report with tool versions and parameters.

Ensure outputs are written in standard formats such as TSV, CSV, JSON, HTML, PNG, SVG, and PDF where appropriate.

⸻

8. Visualisation review and implementation

Review existing visualisation scripts and outputs.

The pipeline should generate interpretable figures such as:

1. Read yield and QC plots.
2. Per-sample assignment success plot.
3. Host spectrum stacked bar chart.
4. HBI/ABI forest plot with confidence intervals.
5. Mixed-feeding frequency bar plot.
6. Marker recovery summary.
7. Contamination or negative-control summary.
8. Optional vector-host network or diversity summaries, if scientifically justified.

Figures should be:

* Automatically generated by the pipeline.
* Clearly labelled.
* Colourblind-safe.
* Saved in publication-quality formats.
* Accompanied by source data tables.
* Not overcomplicated.

If visualisations are scattered across many modules, consider centralising final reporting and plotting in one reporting module or subworkflow. It is acceptable for upstream QC tools to produce their own plots, but thesis/manuscript summary figures should be generated in a single coherent reporting step.

⸻

9. Documentation review and improvement

Review all documentation.

At minimum, the repository should contain:

1. A clear README.md.
2. A pipeline overview.
3. A workflow diagram or text-based workflow schematic.
4. Installation requirements.
5. Quick-start command.
6. Test command.
7. Full input description.
8. Full output description.
9. Parameter documentation.
10. Example samplesheet.
11. Example expected output structure.
12. Troubleshooting section.
13. Citation section.
14. License information.
15. Scientific interpretation boundaries.
16. Known limitations.

The pipeline should also support a useful --help command that lists required and optional parameters, defaults, and example usage.

Implement or improve missing documentation.

⸻

10. Testing and quality control

Check whether the repository includes automated tests.

Look for:

1. nf-test tests.
2. Minimal test data.
3. GitHub Actions or other CI configuration.
4. nf-core lint compatibility, if relevant.
5. Checks for expected output files.
6. Checksums or content-based assertions for key outputs.
7. Tests for invalid input handling.

If tests are absent or weak:

1. Create minimal test data where feasible.
2. Add nf-test tests or a practical equivalent.
3. Add a test profile.
4. Ensure the test workflow runs quickly.
5. Add assertions that confirm the pipeline completes and produces expected outputs.

Run relevant tests after implementation.

⸻

11. nf-core and Nextflow best-practice alignment

Assess the pipeline against nf-core-inspired standards, even if it is not submitted to nf-core.

Check and improve:

1. DSL2 modular structure.
2. Clear parameter schema or equivalent parameter documentation.
3. Standard profile design.
4. Containerisation.
5. Version pinning.
6. CI testing.
7. Linting.
8. README quality.
9. Changelog.
10. License.
11. Consistent output directory structure.
12. Reproducible test profile.
13. Workflow diagram.
14. Tool citation reporting.
15. Semantic versioning.

If the current repository is far from nf-core structure, do not blindly rewrite everything unless justified. First decide whether to:

* Refactor within the current custom DSL2 structure, or
* Scaffold a new nf-core-style structure and migrate logic into it.

Choose the option that gives a professional, maintainable result with the least unnecessary disruption. Explain the trade-off.

⸻

12. Error handling and robustness

Improve robustness by checking:

1. Missing input files.
2. Empty FASTQ files.
3. Malformed samplesheet rows.
4. Duplicate sample IDs.
5. Missing reference database files.
6. Missing primer definitions.
7. Failed external tools.
8. Very low read counts.
9. No host assignment.
10. Negative-control contamination.
11. Barcode cross-talk.
12. Small sample sizes in summary statistics.

Use clear error messages. Where appropriate, fail early.

Implement errorStrategy, maxRetries, and resource escalation for processes likely to fail because of memory or time limits.

⸻

13. Scientific logic and parameter review

Critically evaluate the scientific logic.

Specifically check:

1. Are the read-length filters appropriate for the expected amplicons?
2. Are quality thresholds appropriate for Nanopore reads?
3. Are primer sequences correctly trimmed?
4. Is marker assignment handled correctly?
5. Is the host reference database curated and versioned?
6. Is the taxonomic assignment threshold justified?
7. Is LCA used where species-level certainty is weak?
8. Are mixed-host calls supported by sufficient evidence?
9. Are negative controls used to suppress contamination?
10. Are read fractions reported cautiously?
11. Are HBI and ABI calculated with the correct denominators?
12. Are host-use conclusions framed as realised feeding rather than preference?
13. Are outputs compatible with downstream integration through specimen_id?

If the current approach is scientifically weak, redesign the relevant part of the workflow and document the reason.

⸻

14. Release preparation

Prepare the repository for public release.

Check or implement:

1. LICENSE, preferably MIT unless another license is required.
2. CHANGELOG.md.
3. Semantic versioning.
4. Tool citation documentation.
5. Contributor acknowledgements.
6. GitHub Actions or CI workflow, if appropriate.
7. README badges, if CI exists.
8. Zenodo DOI instructions or .zenodo.json, if useful.
9. Example data and example command.
10. Clear limitations and scope statement.

Do not publish or tag a release unless explicitly instructed. Instead, prepare the repository so it is ready for release.

⸻

15. Execution workflow

Use this work plan:

Phase 1: Inspect

* Map the repository.
* Read core workflow, modules, config, docs, and tests.
* Produce an initial issue list.

Phase 2: Diagnose

Classify issues into:

* Scientific logic
* Nextflow engineering
* Reproducibility
* Configuration
* Documentation
* Testing
* Visualisation
* Release readiness

Phase 3: Implement

Make targeted improvements. Possible changes include:

* Refactoring modules
* Improving channel logic
* Adding validation scripts
* Updating config profiles
* Adding containers or pinned environments
* Improving output structure
* Adding MultiQC or reporting
* Adding plotting scripts
* Adding tests
* Updating README and docs
* Adding license/changelog files

Phase 4: Test

Run the most relevant checks available, such as:

* nextflow run . -profile test
* nextflow run . [test parameters] -resume
* nf-test test
* nf-core lint, if installed and appropriate
* syntax checks for helper scripts
* documentation sanity checks

If a command fails, inspect the error, fix the cause where feasible, and rerun.

Phase 5: Report

At the end, provide a structured final report.

⸻

16. Final report format

Return the final report in this exact structure:

Final Pipeline Review and Refactor Report

1. Executive summary

Briefly state the original condition of the pipeline and the final state after your work.

2. Repository assessment

Summarise the repository structure and key files reviewed.

3. Scientific assessment

Explain whether the pipeline now correctly implements ONT blood-meal host-use metabarcoding for Anopheles mosquitoes.

4. Engineering assessment

Summarise Nextflow DSL2 structure, modularity, configuration, containers, resources, and reproducibility.

5. Major issues found

Table with:
Issue | Category | Severity | Action taken

6. Files modified or created

Table with:
File | Change made | Reason

7. Tests and checks run

Table with:
Command/check | Result | Notes

8. Outputs produced by the pipeline

List key output files and directories.

9. Remaining limitations

List any known limitations that remain, especially scientific assumptions, missing external data, incomplete tests, reference database limitations, or unresolved nf-core compliance issues.

10. Recommended next steps

List concrete next actions before public release.

Be honest. If something could not be completed, say so clearly and explain what blocked it.

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

| Marker | Target size | NUMT risk | Rationale |
|--------|-------------|-----------|-----------|
| `co1_short` | 234 bp | **Moderate** | Borderline above the <200 bp risk threshold (Reeves et al. 2018) |
| `co1_long` | 359 bp | **Low** | Longer amplicon reduces NUMT co-amplification probability |
| `cyt_b` | 359 bp | **Low** | Mitochondrial multicopy advantage confirmed (Hadj-Henni et al. 2015) |

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

| Marker | BOLD/standard species-level threshold | Pipeline global threshold | Assessment |
|--------|--------------------------------------|--------------------------|------------|
| COI (BOLD) | ≥98% (Reeves et al. 2018; Santos et al. 2019) | 97% | Conservative; high-confidence calls (pident ≥98%) can be filtered post-hoc |
| CytB | ≥97–98% (Townzen et al. 2008; Hadj-Henni et al. 2015) | 97% | Appropriate |
| 16S rRNA | ≥90% (Logue et al. 2016) | 97% | More conservative than standard |

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
| **Ecological indices**        | Human Blood Index, animal/zoophily index, feeding-type partition, mixed-feeding rate, host-specific blood indices, and host diversity (richness/Shannon/Simpson/Pielou) with Wilson 95% CIs, stratified by ecological zone, vector sibling species, and collection period/season (`ecological_indices.tsv`) | **Planned:** Fisher's exact + Holm pairwise comparisons between strata (Recommendation 9). **Planned:** Kruskal-Wallis for host richness (Recommendation K). **Planned:** Multi-marker concordance analysis (Recommendation 10).                                                                                                         |
| **Reporting**                 | Custom HÆMA HTML report (primary); optional MultiQC custom-content summary tables; ten publication figures (`07_figures/`, PDF/SVG/PNG) built from the endpoint tables; machine-readable `run_manifest.json`                                                                                             | **Planned:** Concordance heatmap figure (Recommendation 12). **Planned:** Database-source contribution plot. **Planned:** NUMT risk summary plot. **Planned:** HBI/ABI forest plot with corrected denominators.                                                                                                                                             |
