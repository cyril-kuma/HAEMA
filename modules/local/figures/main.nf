process BUILD_FIGURES {
    label 'process_low'
    tag 'figures'

    publishDir "${params.outdir}/07_figures", mode: 'copy'

    input:
    // A single collected bag of endpoint/report tables. Each is staged under its canonical
    // filename (bloodmeal_master_endpoint.tsv, host_call_table.tsv, rambo_model_summary.tsv, ...),
    // so build_figures.py finds them by name via --endpoint-dir . --reports-dir . — and any table
    // whose upstream step was disabled is simply absent and its figure is skipped.
    path endpoint_files

    output:
    path 'figure_*.pdf', emit: pdf
    path 'figure_*.png', emit: png
    path 'figure_*.svg', emit: svg, optional: true
    path 'figure_captions.md', emit: captions
    path 'figure_manifest.tsv', emit: manifest
    path 'versions.yml', emit: versions

    script:
    def formats = (params.figure_formats ?: 'pdf,svg,png')
    """
    python3 ${projectDir}/bin/build_figures.py \\
        --endpoint-dir . \\
        --reports-dir . \\
        --outdir . \\
        --formats '${formats}'

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
      matplotlib: "\$(python3 -c 'import matplotlib; print(matplotlib.__version__)')"
    END_VERSIONS
    """

    stub:
    """
    for n in 01_workflow 02_sequencing_qc 03_depth_denoising 04_host_assignment \\
             05_host_composition 06_mixed_host 07_controls_contamination 08_ecology; do
        printf 'stub' > figure_\${n}.pdf
        printf 'stub' > figure_\${n}.png
    done
    printf '# stub captions\\n' > figure_captions.md
    printf 'figure\\tfiles\\tinputs\\tcaption\\n' > figure_manifest.tsv
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
