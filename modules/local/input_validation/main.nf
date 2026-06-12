process INPUT_VALIDATION {
    label 'process_low'
    tag 'metadata'

    publishDir "${params.outdir}/pipeline_info/input_validation", mode: 'copy'

    input:
    path samplesheet
    path primers
    val raw_data_dir

    output:
    path 'validated_samplesheet.tsv', emit: validated_samplesheet
    path 'sample_fastqs.tsv', emit: sample_fastqs
    path 'marker_config.tsv', emit: marker_config
    path 'control_manifest.tsv', emit: control_manifest
    path 'barcode_filename_mismatches.tsv', emit: barcode_filename_mismatches
    path 'input_validation_report.json', emit: report
    path 'versions.yml', emit: versions

    script:
    """
    python3 ${projectDir}/bin/validate_inputs.py \\
        --samplesheet ${samplesheet} \\
        --primers ${primers} \\
        --raw-data-dir ${raw_data_dir} \\
        --fastq-extensions '${params.fastq_extensions}' \\
        --marker-windows '${params.marker_windows}' \\
        --marker-window-padding '${params.marker_window_padding}' \\
        --require-fastqs '${params.require_fastqs}' \\
        --strict-metadata '${params.strict_metadata}' \\
        --strict-barcode-filenames '${params.strict_barcode_filenames}' \\
        --production-mode '${params.production_mode}' \\
        --barcode-kit '${params.barcode_kit}' \\
        --basecalling-model '${params.basecalling_model}' \\
        --output-dir .

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
    END_VERSIONS
    """

    stub:
    """
    cat > validated_samplesheet.tsv <<-END
    run_id	barcode_id	sample_id	specimen_id	sample_type	species	sibling_species	feeding_status	collection_date	collection_time	collection_location	bioclimatic_zone	collection_region	collection_cordinates	collection_context	collection_method	specimen_sex	sample_uid
    RUN_TEST	barcode01	TST001	TST001	sample	Anopheles_gambiae_s.l	Anopheles_coluzzii	Blood_fed	2026-01-01	00:00	Test	Forest	Test	0,0	Indoor	LTC	Female	RUN_TEST__barcode01__TST001
    END
    cat > sample_fastqs.tsv <<-END
    run_id	barcode_id	sample_id	specimen_id	sample_type	species	sibling_species	feeding_status	collection_date	collection_time	collection_location	bioclimatic_zone	collection_region	collection_cordinates	collection_context	collection_method	specimen_sex	sample_uid	n_fastq_files	fastq_paths
    RUN_TEST	barcode01	TST001	TST001	sample	Anopheles_gambiae_s.l	Anopheles_coluzzii	Blood_fed	2026-01-01	00:00	Test	Forest	Test	0,0	Indoor	LTC	Female	RUN_TEST__barcode01__TST001	1	${projectDir}/assets/test_data/runs/RUN_TEST/20260101_TEST/fastq_pass/barcode01/test.fastq
    END
    cat > marker_config.tsv <<-END
    marker	forward_primer	reverse_primer	expected_size	min_len	max_len
    cyt_b	GAGGMCAAATATCATTCTGAGG	TAGGGCVAGGACTCCTCCTAGT	450	40	140
    co1_short	GCAGGAACAGGWTGAACCG	AATCAGAAYAGGTGTTGGTATAG	324	40	140
    co1_long	AACCACAAAGACATTGGCAC	AAGAATCAGAATARGTGTTG	663	40	160
    END
    cat > control_manifest.tsv <<-END
    sample_uid	run_id	barcode_id	sample_id	sample_type
    END
    cat > barcode_filename_mismatches.tsv <<-END
    run_id	folder_barcode	file_barcode	fastq_path
    END
    printf '{"stub": true}\\n' > input_validation_report.json
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
