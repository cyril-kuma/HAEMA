process BUILD_RUN_MANIFEST {
    label 'process_low'
    tag 'run_manifest'

    publishDir "${params.outdir}/05_endpoint_files", mode: 'copy'

    input:
    path endpoint_manifest
    path input_validation_report
    path production_preflight_report

    output:
    path 'run_manifest.json', emit: manifest
    path 'versions.yml', emit: versions

    script:
    def gpu_active = (workflow.profile ?: '').tokenize(',').contains('gpu')
    // Content hashes of the resolved reference assets so two runs with different reference data are
    // distinguishable from provenance (computed host-side; the files are not staged into the task).
    def fileSha256 = { p ->
        try {
            def f = new File(p?.toString() ?: '')
            (f.exists() && f.isFile()) ? java.security.MessageDigest.getInstance('SHA-256').digest(f.bytes).encodeHex().toString() : ''
        } catch (Exception e) {
            ''
        }
    }
    def reference_checksums = [
        reference_fasta: fileSha256(params.reference_fasta),
        curated_reference_metadata: fileSha256(params.curated_reference_metadata),
        reference_targets: fileSha256(params.reference_targets),
    ]
    def parameters = [
        reference_checksums_sha256: reference_checksums,
        production_mode: params.production_mode,
        input_type: params.input_type,
        barcode_kit: params.barcode_kit,
        flowcell: params.flowcell,
        basecalling_model: params.basecalling_model,
        demux_strategy: params.demux_strategy,
        // Execution-mode and scientific-status flags for provenance / methods reporting.
        medaka_execution: params.enable_medaka ? (gpu_active ? 'enabled_gpu_profile' : 'enabled_cpu') : 'disabled',
        gpu_profile_active: gpu_active,
        mixed_host_thresholds_benchmarked: false,
        host_fractions_benchmarked: false,
        enable_advanced_demux: params.enable_advanced_demux,
        advanced_demux_tool: params.advanced_demux_tool,
        enable_mixed_denoising: params.enable_mixed_denoising,
        mixed_denoise_backend: params.mixed_denoise_backend,
        mixed_denoise_kmer_size: params.mixed_denoise_kmer_size,
        mixed_denoise_min_cluster_size: params.mixed_denoise_min_cluster_size,
        mixed_denoise_min_cluster_fraction: params.mixed_denoise_min_cluster_fraction,
        mixed_denoise_min_reads_for_umap: params.mixed_denoise_min_reads_for_umap,
        mixed_denoise_allow_greedy_fallback: params.mixed_denoise_allow_greedy_fallback,
        consensus_method: params.consensus_method,
        enable_medaka: params.enable_medaka,
        medaka_model: params.medaka_model,
        enable_rambo_model: params.enable_rambo_model,
        taxonomy_strategy: params.taxonomy_strategy,
        taxonomy_assignment_method: params.taxonomy_assignment_method,
        require_curated_taxids: params.require_curated_taxids,
        curated_reference_metadata: params.curated_reference_metadata,
        taxdump_dir: params.taxdump_dir,
        reference_fasta: params.reference_fasta,
        reference_targets: params.reference_targets,
        blast_db: params.blast_db,
        fallback_blast_db: params.fallback_blast_db,
        blast_db_mount: params.blast_db_mount,
        enable_phyloseq: params.enable_phyloseq,
        enable_decontam: params.enable_decontam,
        strict_bioconductor: params.strict_bioconductor,
        decontam_threshold: params.decontam_threshold,
        python_container: params.python_container,
        blast_container: params.blast_container,
        advanced_demux_container: params.advanced_demux_container,
        medaka_container: params.medaka_container,
        r_container: params.r_container,
        multiqc_container: params.multiqc_container,
        outdir: params.outdir,
        log_dir: params.log_dir
    ]
    def outputs = [
        master_endpoint: '05_endpoint_files/bloodmeal_master_endpoint.tsv',
        host_assignments: '05_endpoint_files/host_assignments.tsv',
        host_call_table: params.enable_rambo_model ? '05_endpoint_files/host_call_table.tsv' : '',
        mixed_host_evidence: params.enable_rambo_model ? '04_taxonomy/evidence/mixed_host_evidence.tsv' : '',
        asv_count_table: '05_endpoint_files/asv_count_table.tsv',
        phyloseq_object: params.enable_r_outputs ? '05_endpoint_files/bloodmeal_phyloseq.rds' : '',
        ecology_object: params.enable_r_outputs ? '05_endpoint_files/bloodmeal_ecology_data.rds' : '',
        ecology_object_decontaminated: params.enable_r_outputs ? '05_endpoint_files/bloodmeal_ecology_data_decontaminated.rds' : '',
        host_calls_decontaminated: params.enable_r_outputs ? '05_endpoint_files/host_calls_decontaminated.tsv' : '',
        asv_count_table_decontaminated: params.enable_r_outputs ? '05_endpoint_files/asv_count_table_decontaminated.tsv' : '',
        decontam_results: params.enable_r_outputs ? '06_reports/decontam_results.tsv' : '',
        production_preflight: 'pipeline_info/production_preflight/production_preflight_report.json',
        multiqc_report: params.enable_multiqc ? '06_reports/multiqc_report.html' : '',
        custom_report: '06_reports/bloodmeal_pipeline_report.html'
    ]
    def parameters_b64 = groovy.json.JsonOutput.toJson(parameters).bytes.encodeBase64().toString()
    def outputs_b64 = groovy.json.JsonOutput.toJson(outputs).bytes.encodeBase64().toString()
    """
    python3 - <<'PY'
import base64
from pathlib import Path
Path("parameters.json").write_bytes(base64.b64decode("${parameters_b64}"))
Path("outputs.json").write_bytes(base64.b64decode("${outputs_b64}"))
PY

    python3 ${projectDir}/bin/build_run_manifest.py \\
        --endpoint-manifest '${endpoint_manifest}' \\
        --input-validation-report '${input_validation_report}' \\
        --production-preflight-report '${production_preflight_report}' \\
        --pipeline-version '${workflow.manifest.version ?: "dev"}' \\
        --workflow-run-name '${workflow.runName ?: ""}' \\
        --workflow-session-id '${workflow.sessionId ?: ""}' \\
        --workflow-profile '${workflow.profile ?: ""}' \\
        --workflow-command-line '${workflow.commandLine ?: ""}' \\
        --parameters-json-file parameters.json \\
        --outputs-json-file outputs.json \\
        --output run_manifest.json

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
    END_VERSIONS
    """

    stub:
    """
    printf '{"stub": true}\\n' > run_manifest.json
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
