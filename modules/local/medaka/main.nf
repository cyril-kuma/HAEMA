process MEDAKA_POLISH {
    label 'process_high'
    tag "${meta.run_id}:${meta.barcode_id}:${meta.sample_id}:${meta.marker}:${meta.cluster_id ?: 'all'}"

    publishDir "${params.outdir}/03_consensus_variants/medaka/${meta.run_id}/${meta.barcode_id}/${meta.marker}", mode: 'copy', pattern: '*.fasta'
    publishDir "${params.outdir}/03_consensus_variants/medaka/qc", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(meta), path(feature_fasta), path(feature_counts), path(marker_reads)

    output:
    tuple val(meta), path("${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.medaka.fasta"), path(feature_counts), emit: asvs
    path "${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.medaka_summary.tsv", emit: summary
    path 'versions.yml', emit: versions

    when:
    params.enable_medaka

    script:
    """
    if ! command -v medaka_consensus >/dev/null 2>&1; then
        cat >&2 <<'END_ERROR'
Medaka polishing is enabled but medaka_consensus is not available in the selected container/environment.
Set --enable_medaka false, choose a container with Medaka, or set --medaka_container appropriately.
END_ERROR
        exit 127
    fi

    if [[ ! -s '${feature_fasta}' ]]; then
        : > ${meta.id}.${meta.marker}.${meta.cluster_id ?: "all"}.medaka.fasta
        status="empty_input"
    else
        medaka_consensus \\
            -i '${marker_reads}' \\
            -d '${feature_fasta}' \\
            -o medaka_work \\
            -m '${params.medaka_model}' \\
            -t ${task.cpus} \\
            ${params.medaka_extra_args}

        if [[ -s medaka_work/consensus.fasta ]]; then
            cp medaka_work/consensus.fasta ${meta.id}.${meta.marker}.${meta.cluster_id ?: "all"}.medaka.fasta
            status="polished"
        else
            echo "Medaka completed but did not create medaka_work/consensus.fasta" >&2
            exit 2
        fi
    fi

    cat > ${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.medaka_summary.tsv <<-END
    sample_uid	run_id	sample_id	barcode_id	marker	cluster_id	medaka_model	status	input_fasta	output_fasta
    ${meta.id}	${meta.run_id}	${meta.sample_id}	${meta.barcode_id}	${meta.marker}	${meta.cluster_id ?: 'all'}	${params.medaka_model}	\${status}	${feature_fasta}	${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.medaka.fasta
    END

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      medaka: "\$(medaka --version 2>/dev/null || medaka_consensus --version 2>/dev/null || echo unknown)"
      medaka_model: "${params.medaka_model}"
    END_VERSIONS
    """

    stub:
    """
    cp '${feature_fasta}' ${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.medaka.fasta
    cat > ${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.medaka_summary.tsv <<-END
    sample_uid	run_id	sample_id	barcode_id	marker	cluster_id	medaka_model	status	input_fasta	output_fasta
    ${meta.id}	${meta.run_id}	${meta.sample_id}	${meta.barcode_id}	${meta.marker}	${meta.cluster_id ?: 'all'}	${params.medaka_model}	stub	${feature_fasta}	${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.medaka.fasta
    END
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
