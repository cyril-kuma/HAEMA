process MAKE_BLAST_DB {
    label 'process_medium'
    tag 'reference'

    publishDir "${params.outdir}/04_taxonomy/blast_db", mode: 'copy'

    input:
    path reference_fasta
    path taxid_map

    output:
    tuple path('blastdb'), val('bloodmeal_ref'), emit: db
    path 'versions.yml', emit: versions

    when:
    !params.skip_taxonomy

    script:
    """
    mkdir -p blastdb
    # Taxonomy/seqid design note:
    # The curated panel uses long, pipe-delimited deflines (Species|Accession|...|CommonName).
    # These are fundamentally incompatible with makeblastdb -parse_seqids/-taxid_map:
    #   (1) -parse_seqids fatally fails when a parsed local id exceeds 50 chars
    #       ("the local id is too long ... No volumes were created"), which the real panel hits; and
    #   (2) even under 50 chars, -taxid_map silently attaches no taxids for pipe-delimited ids.
    # Taxonomy is therefore resolved by parse_blast_assignments.py joining the BLAST sseqid to the
    # curated reference-taxonomy sidecar (which carries authoritative taxids). Without -parse_seqids,
    # sseqid is emitted as the full defline, which matches the sidecar key exactly, so -parse_seqids
    # and -taxid_map are intentionally omitted. The taxid_map remains a published provenance artifact.
    # NB: this does NOT disable taxid-backed LCA — the sidecar backfills taxids, so
    # --taxonomy_assignment_method taxid_lca + --taxdump_dir works against the curated panel
    # (proven by tests/test_taxid_assignment.py). The nt fallback supplies its own native staxids.
    makeblastdb -in ${reference_fasta} -dbtype nucl -out blastdb/bloodmeal_ref

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      blast: "\$(makeblastdb -version | head -n 1 | sed 's/^makeblastdb: //')"
      taxid_source: "curated_reference_sidecar"
      taxid_map_records: "\$( [[ -s '${taxid_map}' ]] && (wc -l < '${taxid_map}' | tr -d ' ') || echo 0 )"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p blastdb
    touch blastdb/bloodmeal_ref.nhr blastdb/bloodmeal_ref.nin blastdb/bloodmeal_ref.nsq
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
