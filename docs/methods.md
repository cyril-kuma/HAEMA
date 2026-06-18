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
evidence and amplicon length window.

To resolve mixed blood meals, marker reads were denoised before consensus generation: 5-mer
spectra were embedded with UMAP (McInnes et al., 2018) and clustered with HDBSCAN
(McInnes et al., 2017), retaining clusters with ≥`[mixed_denoise_min_cluster_size]` reads and
≥`[mixed_denoise_min_cluster_fraction]` of the sample/marker reads; groups below
`[mixed_denoise_min_reads_for_umap]` reads used a deterministic greedy identity-clustering
fallback. Per retained cluster, a consensus/representative sequence was generated
(`[consensus_method]`; optional Medaka polishing, Oxford Nanopore Technologies, with model
`[medaka_model]`).

Consensus sequences were assigned taxonomy by BLASTn (Camacho et al., 2009; v`[blast_version]`)
against a curated Ghana/West-African vertebrate mitogenome panel, accepting hits with
≥`[min_blast_identity]`% identity and ≥`[min_blast_coverage]`% query coverage and resolving ties
with a conservative lowest-common-ancestor rule (`[taxonomy_assignment_method]`); taxids were taken
from the curated taxonomy sidecar. Features unresolved by the curated panel were optionally queried
against a local NCBI `nt` database.

Laboratory background was assessed from extraction blanks and PCR negatives: features detected in
negative controls were flagged, and (where R/Bioconductor was available) decontam
(Davis et al., 2018) prevalence calls and negative-control background thresholds were computed.
Controls were retained and flagged in all tables, never silently removed. Host assemblages per
sample/marker were summarised with an abundance/evidence model that preserves multiple hosts above
configurable read and fraction thresholds, classifying each as single- or mixed-host.

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

Final outputs included a per-feature master table, host-call tables, an ASV count matrix, the
ecological-indices table, phyloseq (McMurdie & Holmes, 2013) objects for downstream ecological
analysis, a MultiQC report (Ewels et al., 2016), ten publication figures, and a machine-readable
run manifest of parameters and container digests.

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
