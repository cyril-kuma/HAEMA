process NO_TAXONOMY {
    label 'process_low'
    tag "${meta.run_id}:${meta.barcode_id}:${meta.sample_id}:${meta.marker}:${meta.cluster_id ?: 'all'}"

    publishDir "${params.outdir}/04_taxonomy/assignments/${meta.run_id}/${meta.barcode_id}/${meta.marker}", mode: 'copy', pattern: '*.taxonomy.tsv'

    input:
    tuple val(meta), path(asv_fasta), path(asv_counts)

    output:
    path "${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.taxonomy.tsv", emit: assignments
    path 'versions.yml', emit: versions

    when:
    params.skip_taxonomy

    script:
    """
    python3 ${projectDir}/bin/no_taxonomy.py \\
        --counts '${asv_counts}' \\
        --output '${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.taxonomy.tsv'

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
    END_VERSIONS
    """

    stub:
    """
    python3 ${projectDir}/bin/no_taxonomy.py \\
        --counts '${asv_counts}' \\
        --output '${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.taxonomy.tsv'
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
