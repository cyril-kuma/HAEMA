# Limitations & scientific caveats

HÆMA is research software. This page is the honest, consolidated statement of what is and is not
yet validated, so results are not over-interpreted. Implemented-vs-planned features are tracked in
[`DELIBERATE_LIMITATIONS.md`](DELIBERATE_LIMITATIONS.md); this page focuses on **scientific
defensibility**.

## Implemented and exercised on data
- Metadata/control validation, marker-aware primer trimming, Q20/length filtering, and
  length-based marker splitting (cyt b / short COI / long COI). Verified on real R10.4.1 data
  (~97% of reads passed Q20; primers detected).
- Pre-consensus **mixed-template denoising** with UMAP + HDBSCAN (k-mer spectra), with a
  deterministic greedy fallback for low-read groups.
- Curated **local-first BLAST** with a conservative LCA and an optional external `nt` fallback.
- Control-aware contamination flagging and a RAMBO-style mixed-host evidence/host-call layer.
- phyloseq/decontam endpoints (formal with `haema-r`; documented fallback otherwise).
- End-to-end real-data demonstration (Ghana RUN01): mixed feeding detected in multiple samples
  with biologically plausible hosts (e.g. Human + Goat corroborated across two markers), a clean
  negative control, and a Human-dominant positive control.

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
   recovery, missing, and unexpected hosts plus a `mixed_host_recovery_rate`. The bundled RUN01
   `mf` samples are lab-prepared mixtures and can be used directly as the calibration set once their
   known composition is declared — see [`benchmarking.md`](benchmarking.md). This benchmarks
   **detection/recovery**; quantitative host-*fraction* accuracy still needs controls prepared at
   **known mixing ratios**.

2. **UMAP/HDBSCAN clusters reads, not error-corrected haplotypes.** Clustering raw ONT reads can
   fragment a single host across clusters; downstream identity-consensus and taxonomy collapse
   these, but cluster *counts* are not host counts. Trust `host_call_table.tsv`, not raw cluster
   numbers.

3. **Taxonomy breadth is bounded by the curated panel.** The bundled panel is ~20 vertebrate
   mitogenomes (common Ghanaian/peridomestic hosts). Hosts absent from the panel will be
   unassigned or mis-assigned to the nearest relative unless you supply an `nt` fallback
   (`--fallback_blast_db`). The panel still needs publication-scale expansion and versioning.

4. **Taxids come from the curated sidecar (not embedded in the BLAST DB) — but taxid-LCA still
   works.** The panel's descriptive deflines (`Species|Accession|...`) exceed BLAST's 50-char
   `-parse_seqids` local-id limit, so the DB is built without `-parse_seqids`/`-taxid_map`.
   Assignment joins the BLAST `sseqid` to the curated taxonomy sidecar, which **backfills the
   taxids**, so `--taxonomy_assignment_method taxid_lca --taxdump_dir <dir>` performs genuine
   taxdump-backed lowest-common-ancestor assignment against the curated panel (verified by
   `tests/test_taxid_assignment.py`; the `nt` fallback supplies its own native `staxids`).
   Re-keying the panel to accession seqids is therefore **not required** for taxid-LCA. *Remaining:*
   a frozen, version-labelled DB build and a pinned `taxdump` snapshot for full reproducibility.

5. **Single basecalling assumption.** Input is already-basecalled MinKNOW FASTQ (Dorado SUP,
   R10.4.1). POD5/Dorado re-basecalling and raw-signal demultiplexing are not implemented.

6. **Depth.** Performance/feasibility examples in the docs use depth-capped subsets on a laptop;
   full-depth runs need more RAM/time or HPC. numba JIT recompiles per task (per-work-dir cache),
   which adds fixed overhead to each denoising task.

## Staged / external (not bundled)
- Pattern-aware demultiplexers (Barbell, Deepbinner) — runnable via the command-template wrapper if
  you supply the tool, container, and command; not validated here.
- Medaka polishing — implemented and model-checked at runtime, but mixed-host polishing has not
  been validated against positive controls.
- Production custom images (`haema-python`, `haema-r`) must be built or pulled from a registry.

## What would close these gaps (next scientific steps)
- Sequence and analyse **mixed-host positive controls** to calibrate denoising thresholds and
  confirm Medaka does not merge hosts.
- Expand and **version** the curated panel (checksums, accessions, retrieval dates); add
  marker-specific extracted references.
- Add a pinned NCBI **taxdump** snapshot and migrate panel seqids to accessions to enable native
  taxid LCA.
- Benchmark host-fraction estimates against known mixtures before reporting quantitative
  zoophagic indices.
