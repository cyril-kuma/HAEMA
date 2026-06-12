process ASSIGN_TAXONOMY {
    label 'process_low'
    tag "${meta.run_id}:${meta.barcode_id}:${meta.sample_id}:${meta.marker}:${meta.cluster_id ?: 'all'}"

    publishDir "${params.outdir}/04_taxonomy/assignments/${meta.run_id}/${meta.barcode_id}/${meta.marker}", mode: 'copy', pattern: '*.taxonomy.tsv'

    input:
    tuple val(meta), path(asv_counts), path(blast_tsv), path(reference_taxonomy)

    output:
    path "${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.taxonomy.tsv", emit: assignments
    path 'versions.yml', emit: versions

    script:
    """
    python3 ${projectDir}/bin/parse_blast_assignments.py \\
        --counts '${asv_counts}' \\
        --blast '${blast_tsv}' \\
        --blast-source '${params.reference_db_label}' \\
        --assignment-method '${params.taxonomy_assignment_method}' \\
        --taxdump-dir '${params.taxdump_dir}' \\
        --reference-taxonomy '${reference_taxonomy}' \\
        --min-identity '${params.min_blast_identity}' \\
        --min-coverage '${params.min_blast_coverage}' \\
        --top-bitscore-delta '${params.top_bitscore_delta}' \\
        --output '${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.taxonomy.tsv'

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
    END_VERSIONS
    """

    stub:
    """
    python3 ${projectDir}/bin/parse_blast_assignments.py \\
        --counts '${asv_counts}' \\
        --blast '${blast_tsv}' \\
        --blast-source '${params.reference_db_label}' \\
        --assignment-method '${params.taxonomy_assignment_method}' \\
        --taxdump-dir '${params.taxdump_dir}' \\
        --reference-taxonomy '${reference_taxonomy}' \\
        --min-identity '${params.min_blast_identity}' \\
        --min-coverage '${params.min_blast_coverage}' \\
        --top-bitscore-delta '${params.top_bitscore_delta}' \\
        --output '${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.taxonomy.tsv'
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}

process ASSIGN_TAXONOMY_WITH_FALLBACK {
    label 'process_low'
    tag "${meta.run_id}:${meta.barcode_id}:${meta.sample_id}:${meta.marker}:${meta.cluster_id ?: 'all'}"

    publishDir "${params.outdir}/04_taxonomy/assignments/${meta.run_id}/${meta.barcode_id}/${meta.marker}", mode: 'copy', pattern: '*.taxonomy.tsv'

    input:
    tuple val(meta), path(asv_counts), path(primary_blast_tsv), path(fallback_blast_tsv), path(reference_taxonomy)

    output:
    path "${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.taxonomy.tsv", emit: assignments
    path 'versions.yml', emit: versions

    script:
    """
    python3 ${projectDir}/bin/parse_blast_assignments.py \\
        --counts '${asv_counts}' \\
        --blast '${primary_blast_tsv}' \\
        --blast-source '${params.reference_db_label}' \\
        --fallback-blast '${fallback_blast_tsv}' \\
        --fallback-source '${params.fallback_blast_db_label}' \\
        --assignment-method '${params.taxonomy_assignment_method}' \\
        --taxdump-dir '${params.taxdump_dir}' \\
        --reference-taxonomy '${reference_taxonomy}' \\
        --min-identity '${params.min_blast_identity}' \\
        --min-coverage '${params.min_blast_coverage}' \\
        --top-bitscore-delta '${params.top_bitscore_delta}' \\
        --output '${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.taxonomy.tsv'

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
    END_VERSIONS
    """

    stub:
    """
    python3 ${projectDir}/bin/parse_blast_assignments.py \\
        --counts '${asv_counts}' \\
        --blast '${primary_blast_tsv}' \\
        --blast-source '${params.reference_db_label}' \\
        --fallback-blast '${fallback_blast_tsv}' \\
        --fallback-source '${params.fallback_blast_db_label}' \\
        --assignment-method '${params.taxonomy_assignment_method}' \\
        --taxdump-dir '${params.taxdump_dir}' \\
        --reference-taxonomy '${reference_taxonomy}' \\
        --min-identity '${params.min_blast_identity}' \\
        --min-coverage '${params.min_blast_coverage}' \\
        --top-bitscore-delta '${params.top_bitscore_delta}' \\
        --output '${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.taxonomy.tsv'
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}

process ASSIGN_TAXONOMY_EXTERNAL_DB {
    label 'process_low'
    tag "${meta.run_id}:${meta.barcode_id}:${meta.sample_id}:${meta.marker}:${meta.cluster_id ?: 'all'}"

    publishDir "${params.outdir}/04_taxonomy/assignments/${meta.run_id}/${meta.barcode_id}/${meta.marker}", mode: 'copy', pattern: '*.taxonomy.tsv'

    input:
    tuple val(meta), path(asv_counts), path(blast_tsv)

    output:
    path "${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.taxonomy.tsv", emit: assignments
    path 'versions.yml', emit: versions

    script:
    """
    python3 ${projectDir}/bin/parse_blast_assignments.py \\
        --counts '${asv_counts}' \\
        --blast '${blast_tsv}' \\
        --blast-source '${params.blast_db_label}' \\
        --assignment-method '${params.taxonomy_assignment_method}' \\
        --taxdump-dir '${params.taxdump_dir}' \\
        --min-identity '${params.min_blast_identity}' \\
        --min-coverage '${params.min_blast_coverage}' \\
        --top-bitscore-delta '${params.top_bitscore_delta}' \\
        --output '${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.taxonomy.tsv'

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
    END_VERSIONS
    """

    stub:
    """
    python3 ${projectDir}/bin/parse_blast_assignments.py \\
        --counts '${asv_counts}' \\
        --blast '${blast_tsv}' \\
        --blast-source '${params.blast_db_label}' \\
        --assignment-method '${params.taxonomy_assignment_method}' \\
        --taxdump-dir '${params.taxdump_dir}' \\
        --min-identity '${params.min_blast_identity}' \\
        --min-coverage '${params.min_blast_coverage}' \\
        --top-bitscore-delta '${params.top_bitscore_delta}' \\
        --output '${meta.id}.${meta.marker}.${meta.cluster_id ?: 'all'}.taxonomy.tsv'
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
