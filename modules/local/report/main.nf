process BUILD_REPORT {
    label 'process_low'
    tag 'report'

    publishDir "${params.outdir}/06_reports", mode: 'copy'

    input:
    path master_endpoint
    path sample_summary
    path marker_summary
    path contamination_flags
    path control_check

    output:
    path 'bloodmeal_pipeline_report.html', emit: html
    path 'versions.yml', emit: versions

    script:
    def control_check_arg = control_check.name != 'NO_FILE' ? "--control-check '${control_check}'" : ''
    """
    python3 ${projectDir}/bin/build_report.py \\
        --master-endpoint '${master_endpoint}' \\
        --sample-summary '${sample_summary}' \\
        --marker-summary '${marker_summary}' \\
        --contamination-flags '${contamination_flags}' \\
        ${control_check_arg} \\
        --output bloodmeal_pipeline_report.html

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
    END_VERSIONS
    """

    stub:
    """
    printf '<html><body><h1>HÆMA blood-meal pipeline stub report</h1></body></html>\\n' > bloodmeal_pipeline_report.html
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
