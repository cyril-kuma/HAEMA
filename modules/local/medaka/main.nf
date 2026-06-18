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
    def out_fa = "${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.medaka.fasta"
    def out_sum = "${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.medaka_summary.tsv"
    def cluster = "${meta.cluster_id ?: 'all'}"
    """
    if ! command -v medaka_consensus >/dev/null 2>&1; then
        cat >&2 <<'END_ERROR'
Medaka polishing is enabled but medaka_consensus is not available in the selected container/environment.
Set --enable_medaka false, choose a container with Medaka, or set --medaka_container appropriately.
END_ERROR
        exit 127
    fi

    # Reads available to polish this cluster (marker_reads may be gzipped).
    n_reads=\$(( \$(zcat -f '${marker_reads}' 2>/dev/null | wc -l 2>/dev/null || echo 0) / 4 ))
    status="polished"

    if [[ ! -s '${feature_fasta}' ]]; then
        : > ${out_fa}
        status="empty_input"
    elif [[ "\${n_reads}" -lt ${params.medaka_min_reads} ]]; then
        # Polishing a handful of reads (e.g. negative-control clusters) is meaningless and
        # crash-prone; emit the unpolished consensus unchanged.
        cp '${feature_fasta}' ${out_fa}
        status="skipped_low_reads"
    elif medaka_consensus -i '${marker_reads}' -d '${feature_fasta}' -o medaka_work \\
             -m '${params.medaka_model}' -t ${task.cpus} ${params.medaka_extra_args} > medaka.log 2>&1 \\
         && [[ -s medaka_work/consensus.fasta ]]; then
        cp medaka_work/consensus.fasta ${out_fa}
        status="polished"
    else
        # Medaka failed or produced no consensus (degenerate input, resource kill, model issue):
        # fall back to the unpolished consensus so the run continues — equivalent to the
        # --enable_medaka false result for this cluster, and recorded in the summary status.
        cp '${feature_fasta}' ${out_fa}
        status="fallback_unpolished"
        echo "WARN: medaka polishing failed for ${meta.id} ${meta.marker} ${cluster} (n_reads=\${n_reads}); using unpolished consensus." >&2
        tail -n 8 medaka.log >&2 2>/dev/null || true
    fi

    {
      printf 'sample_uid\\trun_id\\tsample_id\\tbarcode_id\\tmarker\\tcluster_id\\tmedaka_model\\tn_reads\\tstatus\\tinput_fasta\\toutput_fasta\\n'
      printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        '${meta.id}' '${meta.run_id}' '${meta.sample_id}' '${meta.barcode_id}' '${meta.marker}' '${cluster}' '${params.medaka_model}' "\${n_reads}" "\${status}" '${feature_fasta}' '${out_fa}'
    } > ${out_sum}

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
