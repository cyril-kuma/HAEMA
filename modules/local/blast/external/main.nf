process BLASTN_EXTERNAL_ASVS {
    label 'process_medium'
    tag "${meta.run_id}:${meta.barcode_id}:${meta.sample_id}:${meta.marker}:${meta.cluster_id ?: 'all'}"

    publishDir "${params.outdir}/04_taxonomy/raw_blast_external/${meta.run_id}/${meta.barcode_id}/${meta.marker}", mode: 'copy', pattern: '*.blast.tsv'

    input:
    tuple val(meta), path(asv_fasta), path(asv_counts)
    val db_prefix
    val db_label

    output:
    tuple val(meta), path(asv_counts), path("${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.external.blast.tsv"), emit: blast
    path 'versions.yml', emit: versions

    when:
    !params.skip_taxonomy

    script:
    """
    if [[ ! -s '${asv_fasta}' ]]; then
        : > ${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.external.blast.tsv
    else
        blastn \\
            -query '${asv_fasta}' \\
            -db '${db_prefix}' \\
            -max_target_seqs ${params.blast_max_target_seqs} \\
            -evalue ${params.blast_evalue} \\
            ${params.blast_extra_args} \\
            -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore stitle staxids' \\
            -out ${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.external.blast.tsv
    fi

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      blast: "\$(blastn -version | head -n 1 | sed 's/^blastn: //')"
      database_label: "${db_label}"
      database_prefix: "${db_prefix}"
    END_VERSIONS
    """

    stub:
    """
    : > ${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.external.blast.tsv
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
      database_label: "${db_label}"
      database_prefix: "${db_prefix}"
    END
    """
}
