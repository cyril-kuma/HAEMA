process MARKER_CONCORDANCE {
    label 'process_low'
    tag 'marker_concordance'

    publishDir "${params.outdir}/05_endpoint_files", mode: 'copy', pattern: 'marker_concordance.tsv'
    publishDir "${params.outdir}/06_reports", mode: 'copy', pattern: 'marker_concordance_summary.json'

    input:
    path host_calls
    path master_endpoint

    output:
    path 'marker_concordance.tsv', emit: concordance
    path 'marker_concordance_summary.json', emit: summary
    path 'versions.yml', emit: versions

    when:
    params.enable_marker_concordance

    script:
    """
    python3 ${projectDir}/bin/compute_marker_concordance.py \\
        --host-call-table '${host_calls}' \\
        --endpoint-manifest '${master_endpoint}' \\
        --output marker_concordance.tsv \\
        --summary-output marker_concordance_summary.json

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
    END_VERSIONS
    """

    stub:
    """
    printf 'specimen_id\\tsample_id\\tmarkers_with_signal\\thost_calls_by_marker\\tspecies_level_concordance\\tgenus_level_concordance\\tconcordance_status\\tdiscordance_reason\\tpossible_numt_flag\\tpossible_mixed_meal_flag\\n' > marker_concordance.tsv
    printf '{"total_samples": 0}\\n' > marker_concordance_summary.json
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
