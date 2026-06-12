process TRIM_FILTER_SPLIT_MARKERS {
    label 'process_medium'
    tag "${meta.run_id}:${meta.barcode_id}:${meta.sample_id}"

    publishDir "${params.outdir}/02_trimmed_filtered/${meta.run_id}/${meta.barcode_id}", mode: 'copy', pattern: '*.fastq.gz'
    publishDir "${params.outdir}/02_trimmed_filtered/qc", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(meta), path(reads)
    path marker_config

    output:
    tuple val(meta), val('cyt_b'), path("${meta.id}.cyt_b.fastq.gz"), emit: cyt_b
    tuple val(meta), val('co1_short'), path("${meta.id}.co1_short.fastq.gz"), emit: co1_short
    tuple val(meta), val('co1_long'), path("${meta.id}.co1_long.fastq.gz"), emit: co1_long
    path "${meta.id}.trim_filter_split_summary.tsv", emit: summary
    path "${meta.id}.read_decisions.tsv", emit: decisions
    path 'versions.yml', emit: versions

    script:
    """
    python3 ${projectDir}/bin/trim_filter_split.py \\
        --sample-uid '${meta.id}' \\
        --input '${reads}' \\
        --marker-config '${marker_config}' \\
        --output-prefix '${meta.id}' \\
        --min-mean-q '${params.min_mean_q}' \\
        --min-read-length '${params.min_read_length}' \\
        --primer-search-window '${params.primer_search_window}' \\
        --primer-max-error-rate '${params.primer_max_error_rate}' \\
        --allow-length-fallback '${params.allow_length_fallback}'

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
    END_VERSIONS
    """

    stub:
    """
    cp ${reads} ${meta.id}.cyt_b.fastq.gz
    : > ${meta.id}.co1_short.fastq.gz
    : > ${meta.id}.co1_long.fastq.gz
    cat > ${meta.id}.trim_filter_split_summary.tsv <<-END
    sample_uid	marker	raw_reads	pass_quality	assigned_by_primer	assigned_by_length	written_reads	filtered_length	no_marker
    ${meta.id}	cyt_b	1	1	1	0	1	0	0
    END
    cat > ${meta.id}.read_decisions.tsv <<-END
    sample_uid	read_id	status	marker	assignment_method	raw_length	trimmed_length	mean_q
    ${meta.id}	stub	pass	cyt_b	primer_pair	1	1	40
    END
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
