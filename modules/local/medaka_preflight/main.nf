process MEDAKA_MODEL_PREFLIGHT {
    label 'process_low'
    tag 'medaka_model_preflight'

    publishDir "${params.outdir}/pipeline_info/production_preflight", mode: 'copy'

    output:
    path 'medaka_model_preflight.txt', emit: report
    path 'medaka_models.txt', emit: models
    path 'versions.yml', emit: versions

    when:
    params.enable_medaka

    script:
    """
    if ! command -v medaka_consensus >/dev/null 2>&1; then
        echo "Medaka is enabled but medaka_consensus is not available." >&2
        exit 127
    fi

    if medaka tools list_models > medaka_models.txt 2>/dev/null; then
        :
    elif medaka tools list-models > medaka_models.txt 2>/dev/null; then
        :
    else
        echo "Could not list Medaka models with 'medaka tools list_models' or 'list-models'." >&2
        exit 2
    fi

    if ! grep -Fqx '${params.medaka_model}' medaka_models.txt && ! grep -Fq '${params.medaka_model}' medaka_models.txt; then
        echo "Requested Medaka model '${params.medaka_model}' was not found in the selected Medaka environment." >&2
        exit 3
    fi

    cat > medaka_model_preflight.txt <<-END
    status	pass
    requested_model	${params.medaka_model}
    END

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      medaka: "\$(medaka --version 2>/dev/null || medaka_consensus --version 2>/dev/null || echo unknown)"
      medaka_model: "${params.medaka_model}"
    END_VERSIONS
    """

    stub:
    """
    cat > medaka_model_preflight.txt <<-END
    status	stub
    requested_model	${params.medaka_model}
    END
    echo '${params.medaka_model}' > medaka_models.txt
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
