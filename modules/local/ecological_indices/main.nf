process ECOLOGICAL_INDICES {
    label 'process_low'
    tag 'eco_indices'

    publishDir "${params.outdir}/05_endpoint_files", mode: 'copy', pattern: 'ecological_indices.tsv'
    publishDir "${params.outdir}/06_reports", mode: 'copy', pattern: 'ecological_indices_summary.json'

    input:
    path host_calls
    path master_endpoint

    output:
    path 'ecological_indices.tsv', emit: indices
    path 'ecological_indices_summary.json', emit: summary
    path 'versions.yml', emit: versions

    script:
    """
    python3 ${projectDir}/bin/compute_ecological_indices.py \\
        --host-calls '${host_calls}' \\
        --master-endpoint '${master_endpoint}' \\
        --zone-column '${params.eco_index_zone_column}' \\
        --species-column '${params.eco_index_species_column}' \\
        --date-column '${params.eco_index_date_column}' \\
        --wet-months '${params.eco_index_wet_months}' \\
        --output-tsv ecological_indices.tsv \\
        --output-json ecological_indices_summary.json

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
    END_VERSIONS
    """

    stub:
    """
    printf 'stratum_type\\tstratum\\tmetric\\tvalue\\tci_low\\tci_high\\tn\\tdetail\\n' > ecological_indices.tsv
    printf 'overall\\tall_field_samples\\thuman_blood_index\\t0\\t0\\t0\\t0\\tstub\\n' >> ecological_indices.tsv
    printf '{"indices": {}}\\n' > ecological_indices_summary.json
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
