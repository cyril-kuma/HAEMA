process HOST_ECOLOGY_COMPARISONS {
    label 'process_low'
    tag 'host_ecology_comparisons'

    // Exploratory statistical comparisons (Fisher's exact + Holm) — see docs/methods.md.
    // These AUGMENT, and never replace, the descriptive ecological indices.
    publishDir "${params.outdir}/05_endpoint_files/host_ecology_comparisons", mode: 'copy', pattern: '*.tsv'
    publishDir "${params.outdir}/06_reports", mode: 'copy', pattern: 'host_use_statistical_tests_summary.json'

    input:
    path host_calls
    path master_endpoint

    output:
    path '*.tsv', emit: comparisons
    path 'host_use_statistical_tests_summary.json', emit: summary
    path 'versions.yml', emit: versions

    when:
    params.enable_host_ecology_comparisons

    script:
    """
    python3 ${projectDir}/bin/compute_host_ecology_comparisons.py \\
        --host-call-table '${host_calls}' \\
        --master-endpoint '${master_endpoint}' \\
        --zone-column '${params.eco_index_zone_column}' \\
        --species-column '${params.eco_index_species_column}' \\
        --date-column '${params.eco_index_date_column}' \\
        --wet-months '${params.eco_index_wet_months}' \\
        --output-dir .

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
    END_VERSIONS
    """

    stub:
    """
    for f in pairwise_hbi_comparisons pairwise_hbi_species_comparisons pairwise_mixed_feeding_comparisons host_richness_by_zone; do
        printf 'stub\\n' > \${f}.tsv
    done
    printf '{"field_specimens_identified": 0}\\n' > host_use_statistical_tests_summary.json
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
