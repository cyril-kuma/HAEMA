#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INPUT_VALIDATION } from './modules/local/input_validation/main.nf'
include { PRODUCTION_PREFLIGHT } from './modules/local/production_preflight/main.nf'
include { MEDAKA_MODEL_PREFLIGHT } from './modules/local/medaka_preflight/main.nf'
include { R_BIOCONDUCTOR_PREFLIGHT } from './modules/local/r_preflight/main.nf'
include { ADVANCED_DEMUX } from './modules/local/advanced_demux/main.nf'
include { PREPROCESS_READS } from './subworkflows/local/preprocess/main.nf'
include { DENOISE_MIXED_TEMPLATES } from './modules/local/denoise_mixed_templates/main.nf'
include { DEREPLICATE_ASVS } from './modules/local/dereplicate_asvs/main.nf'
include { CLUSTER_CONSENSUS } from './modules/local/cluster_consensus/main.nf'
include { MEDAKA_POLISH } from './modules/local/medaka/main.nf'
include { TAXONOMY_LOCAL_BLAST } from './subworkflows/local/taxonomy/main.nf'
include { AGGREGATE_RESULTS } from './modules/local/aggregate_results/main.nf'
include { BUILD_REPORT } from './modules/local/report/main.nf'
include { BUILD_R_OUTPUTS } from './modules/local/r_outputs/main.nf'
include { MULTIQC_REPORT } from './modules/local/multiqc/main.nf'
include { RAMBO_MIXED_MODEL } from './modules/local/rambo_model/main.nf'
include { BUILD_RUN_MANIFEST } from './modules/local/run_manifest/main.nf'

// Schema-based validation, grouped --help, and a run-parameter summary (nf-schema plugin).
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

workflow {
    // Catch unknown/misspelled/mistyped parameters against nextflow_schema.json before doing work.
    // (--help is handled automatically by nf-schema via the validation.help config.)
    validateParameters()
    log.info paramsSummaryLog(workflow)

    log.info "HÆMA blood-meal metabarcoding pipeline v${workflow.manifest.version ?: 'dev'}"
    log.info "Input mode        : ${params.input_type}"
    log.info "Samplesheet       : ${params.input}"
    log.info "Raw data directory: ${params.raw_data_dir}"
    log.info "Results directory : ${params.outdir}"
    log.info "Production mode    : ${params.production_mode}"
    log.info "Feature gates      : taxonomy=${!params.skip_taxonomy} denoise=${params.enable_mixed_denoising}(${params.mixed_denoise_backend}) " +
             "medaka=${params.enable_medaka} rambo=${params.enable_rambo_model} r_outputs=${params.enable_r_outputs} multiqc=${params.enable_multiqc}"

    // Make reference-panel divergence provenance-visible (its sha256 is recorded in run_manifest.json).
    def bundled_ref = "${projectDir}/assets/references/vertebrate_dna_ref_panel.fasta".toString()
    if (!params.skip_taxonomy && params.reference_fasta?.toString() != bundled_ref) {
        log.warn "Using a non-default --reference_fasta (${params.reference_fasta}); its sha256 is recorded in run_manifest.json for provenance."
    }

    if (!params.input) {
        error "Missing --input samplesheet CSV (or use '-profile test' for the bundled demo; run '--help' for usage)"
    }
    if (!params.primers) {
        error "Missing --primers CSV"
    }
    if (!params.raw_data_dir && !(params.enable_advanced_demux || params.input_type == 'pooled_fastq')) {
        error "Missing --raw_data_dir"
    }
    if (!['fastq', 'pooled_fastq'].contains(params.input_type as String)) {
        error "Unsupported --input_type '${params.input_type}'. Supported: fastq (already-basecalled MinKNOW fastq_pass/barcodeXX) or pooled_fastq."
    }
    if ((params.enable_advanced_demux || params.input_type == 'pooled_fastq') && !params.pooled_fastq) {
        error "Advanced demultiplexing requires --pooled_fastq with one or more pooled FASTQ files."
    }
    if (params.enable_medaka && !(params.consensus_method in ['cluster', 'dereplicate'])) {
        error "Medaka polishing requires consensus_method cluster or dereplicate."
    }
    if (params.enable_r_outputs && !params.enable_rambo_model) {
        error "R/decontam endpoint generation requires --enable_rambo_model true so host calls can be cleaned and exported."
    }

    if (params.enable_advanced_demux || params.input_type == 'pooled_fastq') {
        ch_pooled_fastqs = channel.fromPath(params.pooled_fastq, checkIfExists: true).collect()
        ADVANCED_DEMUX(ch_pooled_fastqs)
        ch_raw_data_dir = ADVANCED_DEMUX.out.demux_root.map { demux_root -> demux_root.toString() }
    } else {
        ch_raw_data_dir = channel.value(file(params.raw_data_dir).toString())
    }

    PRODUCTION_PREFLIGHT(
        file(params.input),
        file(params.primers),
        ch_raw_data_dir
    )
    MEDAKA_MODEL_PREFLIGHT()
    R_BIOCONDUCTOR_PREFLIGHT()

    INPUT_VALIDATION(
        file(params.input),
        file(params.primers),
        ch_raw_data_dir
    )

    ch_samples = INPUT_VALIDATION.out.sample_fastqs
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [:]
            row.each { key, value -> meta[key] = value }
            meta.id = row.sample_uid
            meta.run_name = row.run_id
            meta.barcode = row.barcode_id
            meta.marker = 'mixed'

            def fastqs = row.fastq_paths
                .split(/\|/)
                .findAll { it?.trim() }
                .collect { file(it) }

            tuple(meta, fastqs)
        }

    PREPROCESS_READS(ch_samples, INPUT_VALIDATION.out.marker_config)

    if (params.enable_mixed_denoising) {
        DENOISE_MIXED_TEMPLATES(PREPROCESS_READS.out.marker_reads)
        ch_consensus_input = DENOISE_MIXED_TEMPLATES.out.cluster_reads
            .flatMap { meta, reads ->
                def read_list = reads instanceof List ? reads : [reads]
                read_list.collect { read_file ->
                    def matcher = (read_file.getName() =~ /(cluster[0-9]+)/)
                    def cluster_id = matcher.find() ? matcher.group(1) : 'cluster001'
                    tuple(meta + [cluster_id: cluster_id], read_file)
                }
            }
        // NB: only per-sample/per-cluster *summaries* feed the aggregated qc_summary.tsv.
        // The per-read cluster_membership tables (one row per read => millions on real data)
        // are deliberately excluded here: they bloated qc_summary.tsv to ~1 GB and OOM-killed
        // AGGREGATE_RESULTS / BUILD_R_OUTPUTS (exit 137). They remain published per sample under
        // 03_consensus_variants/mixed_denoising/qc/*.cluster_membership.tsv for diagnostics.
        ch_preprocess_qc_for_aggregate = PREPROCESS_READS.out.qc_tables
            .mix(DENOISE_MIXED_TEMPLATES.out.summaries)
    } else {
        ch_consensus_input = PREPROCESS_READS.out.marker_reads
            .map { meta, reads -> tuple(meta + [cluster_id: 'all'], reads) }
        ch_preprocess_qc_for_aggregate = PREPROCESS_READS.out.qc_tables
    }

    def consensus_method = params.consensus_method as String
    if (!['dereplicate', 'cluster'].contains(consensus_method)) {
        error "Unsupported --consensus_method '${consensus_method}'. Use dereplicate or cluster."
    }

    if (consensus_method == 'cluster') {
        CLUSTER_CONSENSUS(ch_consensus_input)
        ch_feature_asvs = CLUSTER_CONSENSUS.out.asvs
        ch_feature_counts = CLUSTER_CONSENSUS.out.counts
        ch_feature_summaries = CLUSTER_CONSENSUS.out.summaries
    } else {
        DEREPLICATE_ASVS(ch_consensus_input)
        ch_feature_asvs = DEREPLICATE_ASVS.out.asvs
        ch_feature_counts = DEREPLICATE_ASVS.out.counts
        ch_feature_summaries = DEREPLICATE_ASVS.out.summaries
    }

    if (params.enable_medaka) {
        ch_features_keyed = ch_feature_asvs
            .map { meta, fasta, counts -> tuple("${meta.id}|${meta.marker}|${meta.cluster_id ?: 'all'}", meta, fasta, counts) }
        ch_reads_keyed = ch_consensus_input
            .map { meta, reads -> tuple("${meta.id}|${meta.marker}|${meta.cluster_id ?: 'all'}", reads) }
        ch_medaka_input = ch_features_keyed
            .join(ch_reads_keyed)
            .map { _key, meta, fasta, counts, reads -> tuple(meta, fasta, counts, reads) }
        MEDAKA_POLISH(ch_medaka_input)
        ch_taxonomy_asvs = MEDAKA_POLISH.out.asvs
        ch_feature_summaries_for_aggregate = ch_feature_summaries
    } else {
        ch_taxonomy_asvs = ch_feature_asvs
        ch_feature_summaries_for_aggregate = ch_feature_summaries
    }

    TAXONOMY_LOCAL_BLAST(ch_taxonomy_asvs)

    AGGREGATE_RESULTS(
        INPUT_VALIDATION.out.validated_samplesheet,
        ch_feature_counts.collect(),
        TAXONOMY_LOCAL_BLAST.out.assignments.collect(),
        ch_preprocess_qc_for_aggregate.collect(),
        ch_feature_summaries_for_aggregate.collect()
    )

    RAMBO_MIXED_MODEL(
        AGGREGATE_RESULTS.out.master_endpoint,
        AGGREGATE_RESULTS.out.marker_summary
    )

    // Positive-control check feeds the report when the RAMBO host-call model ran; otherwise pass a
    // placeholder so BUILD_REPORT still runs (it degrades to "no positive-control check available").
    ch_pc_check = params.enable_rambo_model
        ? RAMBO_MIXED_MODEL.out.control_check
        : channel.fromPath("${projectDir}/assets/NO_FILE")

    BUILD_REPORT(
        AGGREGATE_RESULTS.out.master_endpoint,
        AGGREGATE_RESULTS.out.sample_summary,
        AGGREGATE_RESULTS.out.marker_summary,
        AGGREGATE_RESULTS.out.contamination_flags,
        ch_pc_check
    )

    BUILD_R_OUTPUTS(
        AGGREGATE_RESULTS.out.master_endpoint,
        AGGREGATE_RESULTS.out.asv_count_table,
        AGGREGATE_RESULTS.out.sample_summary,
        AGGREGATE_RESULTS.out.marker_summary,
        AGGREGATE_RESULTS.out.qc_summary,
        AGGREGATE_RESULTS.out.contamination_flags,
        RAMBO_MIXED_MODEL.out.host_calls
    )

    MULTIQC_REPORT(
        AGGREGATE_RESULTS.out.master_endpoint,
        AGGREGATE_RESULTS.out.sample_summary,
        AGGREGATE_RESULTS.out.marker_summary,
        AGGREGATE_RESULTS.out.qc_summary,
        AGGREGATE_RESULTS.out.contamination_flags
    )

    BUILD_RUN_MANIFEST(
        AGGREGATE_RESULTS.out.manifest,
        INPUT_VALIDATION.out.report,
        PRODUCTION_PREFLIGHT.out.report
    )
}
