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

    ch_qc_tables = MERGE_FASTQ_CHUNKS.out.stats
        .mix(TRIM_FILTER_SPLIT_MARKERS.out.summary)
        .mix(TRIM_FILTER_SPLIT_MARKERS.out.decisions)

    emit:
    marker_reads = ch_marker_reads
    qc_tables = ch_qc_tables
}
