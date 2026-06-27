include { MERGE_FASTQ_CHUNKS } from '../../../modules/local/fastq_ingest/main.nf'
include { TRIM_FILTER_SPLIT_MARKERS } from '../../../modules/local/trim_filter_split/main.nf'

workflow PREPROCESS_READS {
    take:
    ch_samples
    marker_config

    main:
    MERGE_FASTQ_CHUNKS(ch_samples)
    TRIM_FILTER_SPLIT_MARKERS(MERGE_FASTQ_CHUNKS.out.reads, marker_config)

    ch_marker_reads = TRIM_FILTER_SPLIT_MARKERS.out.cyt_b
        .mix(TRIM_FILTER_SPLIT_MARKERS.out.co1_short)
        .mix(TRIM_FILTER_SPLIT_MARKERS.out.co1_long)
        .map { meta, marker, reads ->
            def marker_meta = meta + [marker: marker]
            tuple(marker_meta, reads)
        }

    // Only per-sample/per-marker summary tables feed the aggregated qc_summary.tsv. The per-read
    // read_decisions tables (one row per read => millions on real data) are excluded here to keep
    // qc_summary.tsv small enough for AGGREGATE_RESULTS / BUILD_R_OUTPUTS to load without OOM (exit
    // 137). They remain published per sample under 02_trimmed_filtered/qc/*.read_decisions.tsv.
    ch_qc_tables = MERGE_FASTQ_CHUNKS.out.stats
        .mix(TRIM_FILTER_SPLIT_MARKERS.out.summary)

    emit:
    marker_reads = ch_marker_reads
    qc_tables = ch_qc_tables
    // Per-read decision tables for the read-length figure (Fig 2A). Safe to collect here: the
    // figure step aggregates them into a streamed histogram, never a single in-memory dataframe.
    read_decisions = TRIM_FILTER_SPLIT_MARKERS.out.decisions
}
