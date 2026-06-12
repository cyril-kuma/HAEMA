process R_BIOCONDUCTOR_PREFLIGHT {
    label 'process_low'
    tag 'r_bioconductor_preflight'

    publishDir "${params.outdir}/pipeline_info/production_preflight", mode: 'copy'

    output:
    path 'r_bioconductor_preflight.tsv', emit: report
    path 'versions.yml', emit: versions

    when:
    params.enable_r_outputs && params.strict_bioconductor

    script:
    """
    Rscript - <<'RSCRIPT'
packages <- c("phyloseq", "decontam", "ape", "BiocManager")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
write.table(
  data.frame(package = packages, available = !(packages %in% missing)),
  file = "r_bioconductor_preflight.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
if (length(missing) > 0) {
  stop(sprintf("Missing required R/Bioconductor packages: %s", paste(missing, collapse = ", ")), call. = FALSE)
}
RSCRIPT

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      R: "\$(Rscript --version 2>&1 | sed 's/R scripting front-end version //')"
    END_VERSIONS
    """

    stub:
    """
    cat > r_bioconductor_preflight.tsv <<-END
    package	available
    phyloseq	TRUE
    decontam	TRUE
    ape	TRUE
    BiocManager	TRUE
    END
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
