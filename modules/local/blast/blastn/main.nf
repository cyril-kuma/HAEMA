process BLASTN_ASVS {
    label 'process_medium'
    tag "${meta.run_id}:${meta.barcode_id}:${meta.sample_id}:${meta.marker}:${meta.cluster_id ?: 'all'}"

    publishDir "${params.outdir}/04_taxonomy/raw_blast/${meta.run_id}/${meta.barcode_id}/${meta.marker}", mode: 'copy', pattern: '*.blast.tsv'

    input:
    tuple val(meta), path(asv_fasta), path(asv_counts)
    tuple path(db_dir), val(db_name)

    output:
    tuple val(meta), path(asv_counts), path("${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.blast.tsv"), emit: blast
    path 'versions.yml', emit: versions

    when:
    !params.skip_taxonomy

    script:
    """
    if [[ ! -s '${asv_fasta}' ]]; then
        : > ${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.blast.tsv
    else
        blastn \\
            -query '${asv_fasta}' \\
            -db '${db_dir}/${db_name}' \\
            -max_target_seqs ${params.blast_max_target_seqs} \\
            -evalue ${params.blast_evalue} \\
            ${params.blast_extra_args} \\
            -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore stitle staxids' \\
            -out ${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.blast.tsv
    fi

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      blast: "\$(blastn -version | head -n 1 | sed 's/^blastn: //')"
    END_VERSIONS
    """

    stub:
    """
    : > ${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.blast.tsv
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}

process BLASTN_BOLD_ASVS {
    // BLAST against a runtime-built, project-local BOLD-derived COI database (Mode D).
    // Distinct '.bold.blast.tsv' output so it never clashes with the curated primary
    // when both are staged into ASSIGN_TAXONOMY_WITH_FALLBACK. Effective only for COI
    // markers; CytB queries against a COI DB simply return no hits (harmless).
    label 'process_medium'
    tag "${meta.run_id}:${meta.barcode_id}:${meta.sample_id}:${meta.marker}:${meta.cluster_id ?: 'all'}"

    publishDir "${params.outdir}/04_taxonomy/raw_blast_bold/${meta.run_id}/${meta.barcode_id}/${meta.marker}", mode: 'copy', pattern: '*.blast.tsv'

    input:
    tuple val(meta), path(asv_fasta), path(asv_counts)
    tuple path(db_dir), val(db_name)

    output:
    tuple val(meta), path(asv_counts), path("${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.bold.blast.tsv"), emit: blast
    path 'versions.yml', emit: versions

    when:
    !params.skip_taxonomy

    script:
    """
    if [[ ! -s '${asv_fasta}' ]]; then
        : > ${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.bold.blast.tsv
    else
        blastn \\
            -query '${asv_fasta}' \\
            -db '${db_dir}/${db_name}' \\
            -max_target_seqs ${params.blast_max_target_seqs} \\
            -evalue ${params.blast_evalue} \\
            ${params.blast_extra_args} \\
            -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore stitle staxids' \\
            -out ${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.bold.blast.tsv
    fi

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      blast: "\$(blastn -version | head -n 1 | sed 's/^blastn: //')"
      database_label: "bold_coi"
    END_VERSIONS
    """

    stub:
    """
    : > ${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.bold.blast.tsv
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
