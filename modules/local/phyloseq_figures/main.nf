process PHYLOSEQ_FIGURES {
    label 'process_low'
    tag 'phyloseq_figures'

    publishDir "${params.outdir}/07_figures", mode: 'copy'

    input:
    path phyloseq_rds
    path decontam_results

    output:
    // phyloseq figures are vector PDF + 300 dpi PNG (R-side; SVG is produced for the Python figures).
    // All outputs are optional so the step no-ops cleanly when the .rds is an R-fallback object
    // (no formal phyloseq) — it then emits only phyloseq_figures_SKIPPED.txt + versions.yml.
    path 'figure_1*_phyloseq_*.pdf', emit: pdf, optional: true
    path 'figure_1*_phyloseq_*.png', emit: png, optional: true
    path 'phyloseq_figure_captions.md', emit: captions, optional: true
    path 'phyloseq_figure_manifest.tsv', emit: manifest, optional: true
    path 'phyloseq_figures_SKIPPED.txt', optional: true
    path 'versions.yml', emit: versions

    script:
    """
    Rscript ${projectDir}/bin/build_phyloseq_figures.R \\
        --phyloseq '${phyloseq_rds}' \\
        --decontam-results '${decontam_results}' \\
        --outdir . \\
        --formats 'pdf,png'

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      rscript: "\$(Rscript --version 2>&1 | sed -n 's/.*version \\([0-9.]*\\).*/\\1/p' | head -1)"
    END_VERSIONS
    """

    stub:
    """
    printf 'stub' > figure_11_phyloseq_composition.pdf
    printf 'stub' > figure_11_phyloseq_composition.png
    printf '# stub phyloseq captions\\n' > phyloseq_figure_captions.md
    printf 'figure\\tfiles\\tinputs\\tcaption\\n' > phyloseq_figure_manifest.tsv
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
