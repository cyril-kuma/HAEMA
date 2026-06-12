process BUILD_R_OUTPUTS {
    label 'process_low'
    tag 'r_outputs'

    publishDir "${params.outdir}/05_endpoint_files", mode: 'copy', pattern: '*.rds'
    publishDir "${params.outdir}/05_endpoint_files", mode: 'copy', pattern: '*_decontaminated.tsv'
    publishDir "${params.outdir}/06_reports", mode: 'copy', pattern: '*.tsv'

    input:
    path master_endpoint
    path asv_count_table
    path sample_summary
    path marker_summary
    path qc_summary
    path contamination_flags
    path host_calls

    output:
    path 'bloodmeal_ecology_data.rds', emit: rds
    path 'bloodmeal_ecology_data_decontaminated.rds', emit: decontaminated_rds
    path 'bloodmeal_phyloseq.rds', emit: phyloseq
    path 'decontam_results.tsv', emit: decontam
    path 'host_calls_decontaminated.tsv', emit: decontaminated_host_calls
    path 'asv_count_table_decontaminated.tsv', emit: decontaminated_asv_counts
    path 'qc_background_thresholds.tsv', emit: thresholds
    path 'contamination_model_summary.tsv', emit: contamination_summary
    path 'r_outputs_manifest.tsv', emit: manifest
    path 'versions.yml', emit: versions

    when:
    params.enable_r_outputs

    script:
    """
    Rscript ${projectDir}/bin/build_r_outputs.R \\
        --master-endpoint '${master_endpoint}' \\
        --asv-count-table '${asv_count_table}' \\
        --sample-summary '${sample_summary}' \\
        --marker-summary '${marker_summary}' \\
        --qc-summary '${qc_summary}' \\
        --contamination-flags '${contamination_flags}' \\
        --host-calls '${host_calls}' \\
        --enable-phyloseq '${params.enable_phyloseq}' \\
        --enable-decontam '${params.enable_decontam}' \\
        --strict-bioconductor '${params.strict_bioconductor}' \\
        --decontam-threshold '${params.decontam_threshold}' \\
        --output-dir .

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      R: "\$(Rscript --version 2>&1 | sed 's/R scripting front-end version //')"
    END_VERSIONS
    """

    stub:
    """
    Rscript -e 'saveRDS(list(stub=TRUE), "bloodmeal_ecology_data.rds")'
    Rscript -e 'saveRDS(list(stub=TRUE), "bloodmeal_ecology_data_decontaminated.rds")'
    Rscript -e 'saveRDS(list(stub=TRUE), "bloodmeal_phyloseq.rds")'
    cat > decontam_results.tsv <<-END
    feature_id	contaminant	method	p	neg_prevalence	pos_prevalence	max_negative_count
    END
    cat > qc_background_thresholds.tsv <<-END
    marker	host_assignment	sequence	max_negative_control_count
    END
    cat > host_calls_decontaminated.tsv <<-END
    sample_uid	run_id	sample_id	barcode_id	marker	host_assignment
    END
    cat > asv_count_table_decontaminated.tsv <<-END
    feature_id	marker	sequence
    END
    cat > contamination_model_summary.tsv <<-END
    metric	value
    stub	true
    END
    cat > r_outputs_manifest.tsv <<-END
    output	description
    bloodmeal_ecology_data.rds	stub
    bloodmeal_ecology_data_decontaminated.rds	stub
    bloodmeal_phyloseq.rds	stub
    decontam_results.tsv	stub
    host_calls_decontaminated.tsv	stub
    asv_count_table_decontaminated.tsv	stub
    qc_background_thresholds.tsv	stub
    END
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
