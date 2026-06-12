process CLUSTER_CONSENSUS {
    label 'process_medium'
    tag "${meta.run_id}:${meta.barcode_id}:${meta.sample_id}:${meta.marker}:${meta.cluster_id ?: 'all'}"

    publishDir "${params.outdir}/03_consensus_variants/${meta.run_id}/${meta.barcode_id}/${meta.marker}", mode: 'copy', pattern: '*.fasta'
    publishDir "${params.outdir}/03_consensus_variants/qc", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(meta), path(marker_reads)

    output:
    tuple val(meta), path("${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.consensus.fasta"), path("${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.consensus_counts.tsv"), emit: asvs
    path "${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.consensus_counts.tsv", emit: counts
    path "${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.cluster_summary.tsv", emit: summaries
    path 'versions.yml', emit: versions

    script:
    """
    python3 ${projectDir}/bin/cluster_consensus.py \\
        --sample-uid '${meta.id}' \\
        --run-id '${meta.run_id}' \\
        --sample-id '${meta.sample_id}' \\
        --barcode-id '${meta.barcode_id}' \\
        --marker '${meta.marker}' \\
        --cluster-id '${meta.cluster_id ?: "all"}' \\
        --input '${marker_reads}' \\
        --min-identity '${params.min_cluster_identity}' \\
        --min-reads '${params.min_asv_reads}' \\
        --min-fraction '${params.min_asv_fraction}' \\
        --output-fasta '${meta.id}.${meta.marker}.${meta.cluster_id ?: "all"}.consensus.fasta' \\
        --output-counts '${meta.id}.${meta.marker}.${meta.cluster_id ?: "all"}.consensus_counts.tsv' \\
        --summary '${meta.id}.${meta.marker}.${meta.cluster_id ?: "all"}.cluster_summary.tsv'

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
    END_VERSIONS
    """

    stub:
    """
    cat > ${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.consensus.fasta <<-END
    >${meta.id}|${meta.marker}|${meta.cluster_id ?: 'all'}|CONS0001
    ACGTACGTACGT
    END
    cat > ${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.consensus_counts.tsv <<-END
    sample_uid	run_id	sample_id	barcode_id	marker	cluster_id	asv_id	sequence	count	fraction	retained
    ${meta.id}	${meta.run_id}	${meta.sample_id}	${meta.barcode_id}	${meta.marker}	${meta.cluster_id ?: 'all'}	${meta.id}|${meta.marker}|${meta.cluster_id ?: 'all'}|CONS0001	ACGTACGTACGT	1	1.0	true
    END
    cat > ${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.cluster_summary.tsv <<-END
    sample_uid	run_id	sample_id	barcode_id	marker	cluster_id	n_reads	n_unique	n_retained	mixed_template_warning	method	min_cluster_identity
    ${meta.id}	${meta.run_id}	${meta.sample_id}	${meta.barcode_id}	${meta.marker}	${meta.cluster_id ?: 'all'}	1	1	1	false	greedy_cluster_consensus_no_medaka	${params.min_cluster_identity}
    END
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
