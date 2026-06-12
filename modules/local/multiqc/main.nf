process MULTIQC_REPORT {
    label 'process_low'
    tag 'multiqc'

    publishDir "${params.outdir}/06_reports", mode: 'copy', pattern: 'multiqc_report.html'
    publishDir "${params.outdir}/06_reports/multiqc_data", mode: 'copy', pattern: 'multiqc_data/**'

    input:
    path master_endpoint
    path sample_summary
    path marker_summary
    path qc_summary
    path contamination_flags

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_data', emit: data
    path 'versions.yml', emit: versions

    when:
    params.enable_multiqc

    script:
    """
    mkdir -p multiqc_input
    cp '${master_endpoint}' multiqc_input/
    cp '${sample_summary}' multiqc_input/
    cp '${marker_summary}' multiqc_input/
    cp '${qc_summary}' multiqc_input/
    cp '${contamination_flags}' multiqc_input/

    set +e
    multiqc \\
        multiqc_input \\
        --outdir . \\
        --filename multiqc_report.html \\
        --force
    status=\$?
    set -e

    if [[ ! -f multiqc_report.html ]]; then
        printf '<html><body><h1>MultiQC did not find native tool logs yet</h1><p>Status: %s</p></body></html>\\n' "\$status" > multiqc_report.html
    fi
    mkdir -p multiqc_data

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      multiqc: "\$(multiqc --version | sed 's/multiqc, version //')"
    END_VERSIONS
    """

    stub:
    """
    printf '<html><body><h1>HÆMA MultiQC stub report</h1></body></html>\\n' > multiqc_report.html
    mkdir -p multiqc_data
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
