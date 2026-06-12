process PRODUCTION_PREFLIGHT {
    label 'process_low'
    tag 'production_preflight'

    publishDir "${params.outdir}/pipeline_info/production_preflight", mode: 'copy'

    input:
    path samplesheet
    path primers
    val raw_data_dir

    output:
    path 'production_preflight_report.json', emit: report
    path 'production_preflight_summary.tsv', emit: summary
    path 'versions.yml', emit: versions

    script:
    """
    python3 ${projectDir}/bin/production_preflight.py \\
        --samplesheet '${samplesheet}' \\
        --primers '${primers}' \\
        --raw-data-dir '${raw_data_dir}' \\
        --reference-fasta '${params.reference_fasta}' \\
        --curated-reference-metadata '${params.curated_reference_metadata}' \\
        --reference-targets '${params.reference_targets}' \\
        --blast-db '${params.blast_db}' \\
        --fallback-blast-db '${params.fallback_blast_db}' \\
        --taxdump-dir '${params.taxdump_dir}' \\
        --barcode-kit '${params.barcode_kit}' \\
        --basecalling-model '${params.basecalling_model}' \\
        --flowcell '${params.flowcell}' \\
        --demux-strategy '${params.demux_strategy}' \\
        --production-mode '${params.production_mode}' \\
        --enable-medaka '${params.enable_medaka}' \\
        --medaka-model '${params.medaka_model}' \\
        --enable-r-outputs '${params.enable_r_outputs}' \\
        --strict-bioconductor '${params.strict_bioconductor}' \\
        --r-container '${params.r_container}' \\
        --python-container '${params.python_container}' \\
        --medaka-container '${params.medaka_container}' \\
        --output-dir .

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
    END_VERSIONS
    """

    stub:
    """
    cat > production_preflight_report.json <<-END
    {"stub": true, "production_mode": ${params.production_mode}}
    END
    cat > production_preflight_summary.tsv <<-END
    check	status	severity	message	evidence
    stub	pass	info	stub	preflight
    END
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
