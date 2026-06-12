include { MAKE_BLAST_DB } from '../../../modules/local/blast/makeblastdb/main.nf'
include { BUILD_CURATED_TAXID_MAP } from '../../../modules/local/curated_reference_metadata/main.nf'
include { BLASTN_ASVS } from '../../../modules/local/blast/blastn/main.nf'
include { BLASTN_EXTERNAL_ASVS as BLASTN_EXTERNAL_PRIMARY; BLASTN_EXTERNAL_ASVS as BLASTN_EXTERNAL_FALLBACK } from '../../../modules/local/blast/external/main.nf'
include { ASSIGN_TAXONOMY as ASSIGN_TAXONOMY_PRIMARY; ASSIGN_TAXONOMY_WITH_FALLBACK; ASSIGN_TAXONOMY_EXTERNAL_DB } from '../../../modules/local/taxonomy_assign/main.nf'
include { NO_TAXONOMY } from '../../../modules/local/no_taxonomy/main.nf'

workflow TAXONOMY_LOCAL_BLAST {
    take:
    ch_asvs

    main:
    if (params.skip_taxonomy) {
        NO_TAXONOMY(ch_asvs)
        ch_assignments = NO_TAXONOMY.out.assignments
    } else {
        def strategy = params.taxonomy_strategy as String
        if (!['curated_then_fallback', 'curated_only', 'nt_only'].contains(strategy)) {
            error "Unsupported --taxonomy_strategy '${strategy}'. Use curated_then_fallback, curated_only, or nt_only."
        }

        if (strategy == 'nt_only') {
            if (!params.blast_db) {
                error "Taxonomy strategy nt_only requires --blast_db with a BLAST database prefix."
            }
            BLASTN_EXTERNAL_PRIMARY(ch_asvs, params.blast_db as String, params.blast_db_label as String)
            ASSIGN_TAXONOMY_EXTERNAL_DB(BLASTN_EXTERNAL_PRIMARY.out.blast)
            ch_assignments = ASSIGN_TAXONOMY_EXTERNAL_DB.out.assignments
        } else {
            if (!params.reference_fasta) {
                error "Taxonomy strategy ${strategy} requires --reference_fasta for the curated local panel."
            }
            if (params.require_curated_taxids && !params.curated_reference_metadata) {
                error "Curated taxid-backed classification requires --curated_reference_metadata when --require_curated_taxids true."
            }
            if (params.require_taxdump_for_lca && (params.taxonomy_assignment_method in ['taxid_lca', 'conservative_lca']) && !params.taxdump_dir) {
                error "Taxid-backed LCA requires --taxdump_dir when --require_taxdump_for_lca true."
            }
            def curated_metadata = params.curated_reference_metadata
                ? file(params.curated_reference_metadata)
                : file("${projectDir}/assets/empty_curated_reference_metadata.tsv")
            BUILD_CURATED_TAXID_MAP(file(params.reference_fasta), curated_metadata)
            MAKE_BLAST_DB(file(params.reference_fasta), BUILD_CURATED_TAXID_MAP.out.taxid_map)
            ch_db = MAKE_BLAST_DB.out.db
            ch_curated_taxonomy = BUILD_CURATED_TAXID_MAP.out.taxonomy_table
            BLASTN_ASVS(ch_asvs, ch_db)

            if (strategy == 'curated_then_fallback' && params.fallback_blast_db) {
                BLASTN_EXTERNAL_FALLBACK(ch_asvs, params.fallback_blast_db as String, params.fallback_blast_db_label as String)
                ch_primary = BLASTN_ASVS.out.blast
                    .map { meta, counts, blast -> tuple("${meta.id}|${meta.marker}|${meta.cluster_id ?: 'all'}", meta, counts, blast) }
                ch_fallback = BLASTN_EXTERNAL_FALLBACK.out.blast
                    .map { meta, _counts, blast -> tuple("${meta.id}|${meta.marker}|${meta.cluster_id ?: 'all'}", blast) }
                ch_joined = ch_primary
                    .join(ch_fallback)
                    .map { _key, meta, counts, primary_blast, fallback_blast -> tuple(meta, counts, primary_blast, fallback_blast) }
                ch_joined_with_taxonomy = ch_joined
                    .combine(ch_curated_taxonomy)
                    .map { meta, counts, primary_blast, fallback_blast, taxonomy_table -> tuple(meta, counts, primary_blast, fallback_blast, taxonomy_table) }
                ASSIGN_TAXONOMY_WITH_FALLBACK(ch_joined_with_taxonomy)
                ch_assignments = ASSIGN_TAXONOMY_WITH_FALLBACK.out.assignments
            } else {
                if (strategy == 'curated_then_fallback' && !params.fallback_blast_db) {
                    log.warn "Taxonomy strategy curated_then_fallback selected but --fallback_blast_db is empty; running curated_only."
                }
                ch_primary_with_taxonomy = BLASTN_ASVS.out.blast
                    .combine(ch_curated_taxonomy)
                    .map { meta, counts, blast, taxonomy_table -> tuple(meta, counts, blast, taxonomy_table) }
                ASSIGN_TAXONOMY_PRIMARY(ch_primary_with_taxonomy)
                ch_assignments = ASSIGN_TAXONOMY_PRIMARY.out.assignments
            }
        }
    }

    emit:
    assignments = ch_assignments
}
