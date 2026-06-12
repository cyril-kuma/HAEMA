process RAMBO_MIXED_MODEL {
    label 'process_low'
    tag 'mixed_host_model'

    publishDir "${params.outdir}/05_endpoint_files", mode: 'copy', pattern: 'host_call_table.tsv'
    publishDir "${params.outdir}/04_taxonomy/evidence", mode: 'copy', pattern: 'mixed_host_evidence.tsv'
    publishDir "${params.outdir}/06_reports", mode: 'copy', pattern: 'rambo_model_summary.tsv'
    publishDir "${params.outdir}/06_reports", mode: 'copy', pattern: 'positive_control_check.tsv'

    input:
    path master_endpoint
    path marker_summary

    output:
    path 'mixed_host_evidence.tsv', emit: evidence
    path 'host_call_table.tsv', emit: host_calls
    path 'rambo_model_summary.tsv', emit: summary
    path 'positive_control_check.tsv', emit: control_check
    path 'versions.yml', emit: versions

    when:
    params.enable_rambo_model

    script:
    """
    python3 ${projectDir}/bin/rambo_mixed_model.py \\
        --master-endpoint '${master_endpoint}' \\
        --marker-summary '${marker_summary}' \\
        --min-host-reads '${params.rambo_min_host_reads}' \\
        --min-host-fraction '${params.rambo_min_host_fraction}' \\
        --include-contaminants '${params.rambo_include_contaminants}' \\
        --output-evidence mixed_host_evidence.tsv \\
        --output-host-calls host_call_table.tsv \\
        --output-summary rambo_model_summary.tsv \\
        --output-control-check positive_control_check.tsv

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
      model: "rambo_style_abundance_evidence"
    END_VERSIONS
    """

    stub:
    """
    cat > mixed_host_evidence.tsv <<-END
    sample_uid	run_id	sample_id	barcode_id	control_status	marker	cluster_id	asv_id	host_assignment	taxon_rank	count	feature_fraction	sample_marker_fraction	confidence	assignment_status	blast_source	contamination_flag	retained_for_mixed_model	model_reason
    RUN_TEST__barcode01__TST001	RUN_TEST	TST001	barcode01	sample	cyt_b	cluster001	stub	Homo sapiens	species	1	1.0	1.0	high	assigned	stub	false	true	assigned_noncontaminant
    END
    cat > host_call_table.tsv <<-END
    sample_uid	run_id	sample_id	barcode_id	control_status	marker	best_cluster_id	host_assignment	host_rank	host_reads	host_fraction	n_supporting_features	best_feature_id	best_confidence	best_assignment_status	mixed_status	total_marker_reads	mixed_template_warning	model
    RUN_TEST__barcode01__TST001	RUN_TEST	TST001	barcode01	sample	cyt_b	cluster001	Homo sapiens	1	1	1.0	1	stub	high	assigned	single_host	1	false	rambo_style_abundance_evidence
    END
    cat > rambo_model_summary.tsv <<-END
    metric	value
    sample_marker_groups	1
    single_host_groups	1
    mixed_host_groups	0
    positive_controls_total	0
    host_fractions_benchmarked	false
    END
    cat > positive_control_check.tsv <<-END
    sample_uid	control_kind	control_status	expected_hosts	observed_hosts	n_expected	n_recovered	missing_hosts	unexpected_hosts	markers_with_signal	status	note
    END
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
