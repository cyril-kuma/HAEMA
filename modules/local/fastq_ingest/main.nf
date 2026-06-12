process MERGE_FASTQ_CHUNKS {
    label 'process_medium'
    tag "${meta.run_id}:${meta.barcode_id}:${meta.sample_id}"

    publishDir "${params.outdir}/01_ingested/${meta.run_id}/${meta.barcode_id}", mode: 'copy', pattern: '*.fastq.gz'
    publishDir "${params.outdir}/01_ingested/qc", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), path("${meta.id}.merged.fastq.gz"), emit: reads
    path "${meta.id}.merge_stats.tsv", emit: stats
    path 'versions.yml', emit: versions

    script:
    """
    python3 ${projectDir}/bin/merge_fastqs.py \\
        --sample-uid '${meta.id}' \\
        --run-id '${meta.run_id}' \\
        --sample-id '${meta.sample_id}' \\
        --barcode-id '${meta.barcode_id}' \\
        --output '${meta.id}.merged.fastq.gz' \\
        --stats '${meta.id}.merge_stats.tsv' \\
        ${fastqs}

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
    END_VERSIONS
    """

    stub:
    """
    cp ${fastqs[0]} ${meta.id}.merged.fastq.gz
    cat > ${meta.id}.merge_stats.tsv <<-END
    sample_uid	run_id	sample_id	barcode_id	n_files	reads	bases	min_len	max_len	mean_len	mean_q
    ${meta.id}	${meta.run_id}	${meta.sample_id}	${meta.barcode_id}	1	1	1	1	1	1	40
    END
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
