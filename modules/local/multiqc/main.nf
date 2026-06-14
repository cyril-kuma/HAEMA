process MULTIQC_REPORT {
    label 'process_low'
    tag 'multiqc'

    publishDir "${params.outdir}/06_reports", mode: 'copy', pattern: 'multiqc_report.html'
    publishDir "${params.outdir}/06_reports", mode: 'copy', pattern: 'multiqc_data'

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

    # MultiQC does not natively parse HÆMA's custom TSVs, so the compact per-sample and per-marker
    # summaries are exposed as MultiQC *custom-content* tables: a file ending in *_mqc.tsv whose
    # first lines are a '# key: value' comment header is auto-ingested as a report section. This
    # makes the MultiQC report a genuine overview instead of an empty placeholder. The large master
    # endpoint and per-read qc_summary are deliberately NOT fed to MultiQC (too wide/long); they
    # remain the primary tables under 05_endpoint_files/.
    {
        echo '# id: "haema_sample_summary"'
        echo '# section_name: "HÆMA: per-sample summary"'
        echo '# description: "Reads, retained features, and host calls per sample and control."'
        echo '# plot_type: "table"'
        cat '${sample_summary}'
    } > multiqc_input/haema_sample_summary_mqc.tsv

    # marker_level_summary has one row per sample/marker, so prepend a unique sample.marker key
    # (MultiQC tables key rows on the first column).
    {
        echo '# id: "haema_marker_summary"'
        echo '# section_name: "HÆMA: per-marker summary"'
        echo '# description: "Per-sample/marker read and feature counts (cyt b / short COI / long COI)."'
        echo '# plot_type: "table"'
        awk -F'\\t' 'NR==1{print "sample_marker\\t"\$0} NR>1{print \$1"."\$2"\\t"\$0}' '${marker_summary}'
    } > multiqc_input/haema_marker_summary_mqc.tsv

    set +e
    multiqc multiqc_input --outdir . --filename multiqc_report.html --title 'HÆMA report' --force
    status=\$?
    set -e

    # MultiQC names its data dir after the report file; normalise it to the published 'multiqc_data'.
    for d in multiqc_report_data multiqc_report.html_data; do
        [ -d "\$d" ] && mv "\$d" multiqc_data
    done
    mkdir -p multiqc_data

    if [ ! -f multiqc_report.html ]; then
        printf '<html><body><h1>HÆMA: MultiQC report</h1><p>MultiQC exited %s without producing a report; use the primary HÆMA report (bloodmeal_pipeline_report.html) and the endpoint tables under 05_endpoint_files/.</p></body></html>\\n' "\$status" > multiqc_report.html
    fi

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
