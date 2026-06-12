process AGGREGATE_RESULTS {
    label 'process_low'
    tag 'endpoint'

    publishDir "${params.outdir}/05_endpoint_files", mode: 'copy'

    input:
    path validated_samplesheet
    path asv_counts
    path taxonomy_assignments
    path preprocess_qc
    path derep_summaries

    output:
    path 'bloodmeal_master_endpoint.tsv', emit: master_endpoint
    path 'host_assignments.tsv', emit: host_assignments
    path 'asv_count_table.tsv', emit: asv_count_table
    path 'sample_level_summary.tsv', emit: sample_summary
    path 'marker_level_summary.tsv', emit: marker_summary
    path 'qc_summary.tsv', emit: qc_summary
    path 'contamination_flags.tsv', emit: contamination_flags
    path 'endpoint_manifest.json', emit: manifest
    path 'versions.yml', emit: versions

    script:
    """
    python3 ${projectDir}/bin/aggregate_results.py \\
        --samplesheet '${validated_samplesheet}' \\
        --counts ${asv_counts} \\
        --assignments ${taxonomy_assignments} \\
        --preprocess-qc ${preprocess_qc} \\
        --derep-summaries ${derep_summaries} \\
        --negative-control-multiplier '${params.negative_control_multiplier}' \\
        --output-dir .

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
    END_VERSIONS
    """

    stub:
    """
    cat > bloodmeal_master_endpoint.tsv <<-END
    sample_uid	run_id	sample_id	barcode_id	marker	asv_id	host_assignment	assignment_status	count	control_status	contamination_flag
    RUN_TEST__barcode01__TST001	RUN_TEST	TST001	barcode01	cyt_b	RUN_TEST__barcode01__TST001|cyt_b|ASV0001	unassigned	unassigned	1	sample	false
    END
    cp bloodmeal_master_endpoint.tsv host_assignments.tsv
    cat > asv_count_table.tsv <<-END
    feature_id	marker	sequence	RUN_TEST__barcode01__TST001
    cyt_b_000000000000	cyt_b	ACGTACGTACGT	1
    END
    cat > sample_level_summary.tsv <<-END
    sample_uid	run_id	sample_id	barcode_id	sample_type	total_asv_reads	retained_asvs	assigned_hosts	contamination_flags
    RUN_TEST__barcode01__TST001	RUN_TEST	TST001	barcode01	sample	1	1	0	0
    END
    cat > marker_level_summary.tsv <<-END
    sample_uid	marker	total_asv_reads	retained_asvs	assigned_hosts	mixed_template_warning
    RUN_TEST__barcode01__TST001	cyt_b	1	1	0	false
    END
    cat > qc_summary.tsv <<-END
    source_file	metric	value
    stub	stub	true
    END
    cat > contamination_flags.tsv <<-END
    sample_uid	marker	asv_id	contamination_flag	reason
    END
    printf '{"stub": true}\\n' > endpoint_manifest.json
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
