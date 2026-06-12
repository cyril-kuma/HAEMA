process BUILD_CURATED_TAXID_MAP {
    label 'process_low'
    tag 'curated_taxids'

    publishDir "${params.outdir}/04_taxonomy/curated_reference", mode: 'copy'

    input:
    path reference_fasta
    path metadata

    output:
    path 'curated_taxid_map.tsv', emit: taxid_map
    path 'curated_reference_taxonomy.tsv', emit: taxonomy_table
    path 'curated_reference_taxid_report.json', emit: report
    path 'versions.yml', emit: versions

    when:
    !params.skip_taxonomy

    script:
    """
    python3 ${projectDir}/bin/build_curated_taxid_map.py \\
        --reference-fasta '${reference_fasta}' \\
        --metadata '${metadata}' \\
        --require-taxids '${params.require_curated_taxids}' \\
        --taxid-map curated_taxid_map.tsv \\
        --taxonomy-table curated_reference_taxonomy.tsv \\
        --report curated_reference_taxid_report.json

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
    END_VERSIONS
    """

    stub:
    """
    cat > curated_taxid_map.tsv <<-END
    test_host|taxid_9606	9606
    END
    cat > curated_reference_taxonomy.tsv <<-END
    seqid	taxid	scientific_name	rank	common_name	source_accession	provenance
    test_host|taxid_9606	9606	Homo sapiens	species	human	test	stub
    END
    printf '{"stub": true}\\n' > curated_reference_taxid_report.json
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
