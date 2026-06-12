process DENOISE_MIXED_TEMPLATES {
    label 'process_high'
    tag "${meta.run_id}:${meta.barcode_id}:${meta.sample_id}:${meta.marker}"

    publishDir "${params.outdir}/03_consensus_variants/mixed_denoising/${meta.run_id}/${meta.barcode_id}/${meta.marker}", mode: 'copy', pattern: '*.clusters/*.fastq.gz'
    publishDir "${params.outdir}/03_consensus_variants/mixed_denoising/qc", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(meta), path(marker_reads)

    output:
    tuple val(meta), path("${meta.id}.${meta.marker}.clusters/*.fastq.gz"), emit: cluster_reads
    path "${meta.id}.${meta.marker}.mixed_denoise_summary.tsv", emit: summaries
    path "${meta.id}.${meta.marker}.cluster_membership.tsv", emit: membership
    path 'versions.yml', emit: versions

    when:
    params.enable_mixed_denoising

    script:
    """
    # UMAP relies on numba JIT, which must cache compiled functions to a writable location.
    # Under Nextflow's non-root container execution (-u uid:gid) the default in-tree/HOME cache
    # is unwritable, which makes umap import fail ("cannot cache function 'rdist'...") and forces
    # the greedy fallback. Point numba/HOME at the writable task work dir so real UMAP/HDBSCAN runs.
    export NUMBA_CACHE_DIR="\${PWD}/.numba_cache"
    export HOME="\${PWD}"
    mkdir -p "\${NUMBA_CACHE_DIR}" '${meta.id}.${meta.marker}.clusters'
    python3 ${projectDir}/bin/denoise_mixed_templates.py \\
        --sample-uid '${meta.id}' \\
        --run-id '${meta.run_id}' \\
        --sample-id '${meta.sample_id}' \\
        --barcode-id '${meta.barcode_id}' \\
        --marker '${meta.marker}' \\
        --input '${marker_reads}' \\
        --backend '${params.mixed_denoise_backend}' \\
        --kmer-size '${params.mixed_denoise_kmer_size}' \\
        --min-cluster-size '${params.mixed_denoise_min_cluster_size}' \\
        --min-cluster-fraction '${params.mixed_denoise_min_cluster_fraction}' \\
        --min-reads-for-umap '${params.mixed_denoise_min_reads_for_umap}' \\
        --greedy-min-identity '${params.min_cluster_identity}' \\
        --allow-greedy-fallback '${params.mixed_denoise_allow_greedy_fallback}' \\
        --output-dir '${meta.id}.${meta.marker}.clusters' \\
        --summary '${meta.id}.${meta.marker}.mixed_denoise_summary.tsv' \\
        --membership '${meta.id}.${meta.marker}.cluster_membership.tsv'

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
      backend_requested: "${params.mixed_denoise_backend}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p '${meta.id}.${meta.marker}.clusters'
    cp '${marker_reads}' '${meta.id}.${meta.marker}.clusters/${meta.id}.${meta.marker}.cluster001.fastq.gz'
    cat > ${meta.id}.${meta.marker}.mixed_denoise_summary.tsv <<-END
    sample_uid	run_id	sample_id	barcode_id	marker	input_reads	retained_reads	noise_reads	n_clusters	cluster_read_counts	backend_requested	backend_used	fallback_used	fallback_reason	mixed_template_warning	min_cluster_size	min_cluster_fraction	kmer_size
    ${meta.id}	${meta.run_id}	${meta.sample_id}	${meta.barcode_id}	${meta.marker}	1	1	0	1	cluster001:1	stub	stub	false		false	${params.mixed_denoise_min_cluster_size}	${params.mixed_denoise_min_cluster_fraction}	${params.mixed_denoise_kmer_size}
    END
    cat > ${meta.id}.${meta.marker}.cluster_membership.tsv <<-END
    sample_uid	run_id	sample_id	barcode_id	marker	read_id	raw_cluster_label	cluster_id	retained	reason	read_length	backend_requested	backend_used
    ${meta.id}	${meta.run_id}	${meta.sample_id}	${meta.barcode_id}	${meta.marker}	stub_read	0	cluster001	true	retained	12	stub	stub
    END
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
