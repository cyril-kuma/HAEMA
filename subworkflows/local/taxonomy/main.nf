include { MAKE_BLAST_DB as MAKE_CURATED_DB; MAKE_BLAST_DB as MAKE_BOLD_DB } from '../../../modules/local/blast/makeblastdb/main.nf'
include { BUILD_CURATED_TAXID_MAP } from '../../../modules/local/curated_reference_metadata/main.nf'
include { BLASTN_ASVS; BLASTN_BOLD_ASVS } from '../../../modules/local/blast/blastn/main.nf'
include { BLASTN_EXTERNAL_ASVS as BLASTN_EXTERNAL_PRIMARY; BLASTN_EXTERNAL_ASVS as BLASTN_EXTERNAL_FALLBACK } from '../../../modules/local/blast/external/main.nf'
include { ASSIGN_TAXONOMY as ASSIGN_TAXONOMY_PRIMARY; ASSIGN_TAXONOMY_WITH_FALLBACK; ASSIGN_TAXONOMY_EXTERNAL_DB } from '../../../modules/local/taxonomy_assign/main.nf'
include { NO_TAXONOMY } from '../../../modules/local/no_taxonomy/main.nf'

// ---------------------------------------------------------------------------
// Reference database architecture (re-engineered).
//
// Four reference modes, selected by params.reference_mode. The curated vertebrate
// panel is preserved as an optional first-line reference, NOT as the sole limiting
// reference strategy:
//
//   Mode A  curated_panel   - curated vertebrate panel only (fast, offline).
//   Mode B  broad_blast     - user-supplied broad local BLAST DB (--blast_db),
//                             optionally with the curated panel as a pre-check.
//   Mode C  remote_fallback - curated panel primary, NCBI nt remote fallback for
//                             unresolved features (reproducibility caveat; off by
//                             default via --enable_ncbi_remote_fallback).
//   Mode D  bold_aware      - curated panel primary, BOLD-derived COI FASTA DB
//                             fallback (reproducible local-first; effective for COI).
//
// reference_mode is the canonical selector. The legacy --taxonomy_strategy
// (curated_only | curated_then_fallback | nt_only) is still honoured and mapped
// onto a reference_mode when reference_mode is left at its default, so existing
// configs and tests keep working unchanged.
// ---------------------------------------------------------------------------

def keyOf(meta) {
    "${meta.id}|${meta.marker}|${meta.cluster_id ?: 'all'}"
}

workflow TAXONOMY_LOCAL_BLAST {
    take:
    ch_asvs

    main:
    def ch_assignments

    if (params.skip_taxonomy) {
        NO_TAXONOMY(ch_asvs)
        ch_assignments = NO_TAXONOMY.out.assignments
    } else {
        // --- Resolve the effective reference mode --------------------------
        def mode = params.reference_mode as String
        def curated_check = params.enable_curated_panel_check as boolean
        def legacy = params.taxonomy_strategy as String

        // Back-compat: only consult the legacy strategy when reference_mode is default.
        if (mode == 'curated_panel') {
            if (legacy == 'nt_only') {
                mode = 'broad_blast'
                curated_check = false            // legacy nt_only never ran the curated panel
            } else if (legacy == 'curated_then_fallback' && params.fallback_blast_db) {
                mode = 'broad_blast'             // curated primary + local external fallback DB
                curated_check = true
            }
            // curated_only, or curated_then_fallback without a fallback DB -> curated_panel
        }

        if (!['curated_panel', 'broad_blast', 'remote_fallback', 'bold_aware'].contains(mode)) {
            error "Unsupported --reference_mode '${mode}'. Use curated_panel, broad_blast, remote_fallback, or bold_aware."
        }
        if (params.require_taxdump_for_lca && (params.taxonomy_assignment_method in ['taxid_lca', 'conservative_lca']) && !params.taxdump_dir) {
            error "Taxid-backed LCA requires --taxdump_dir when --require_taxdump_for_lca true."
        }

        // Mode B without a curated pre-check is the only branch that does NOT build the
        // curated panel; every other path runs the curated panel as the primary reference.
        def broad_only = (mode == 'broad_blast' && !curated_check)
        def needs_curated = !broad_only

        if (needs_curated && !params.reference_fasta) {
            error "reference_mode '${mode}' uses the curated panel but --reference_fasta is not set."
        }

        // Whether unresolved curated hits are escalated to a fallback DB, and how.
        def use_fallback = false
        def fallback_label = ''

        if (needs_curated) {
            if (params.require_curated_taxids && !params.curated_reference_metadata) {
                error "Curated taxid-backed classification requires --curated_reference_metadata when --require_curated_taxids true."
            }
            def curated_metadata = params.curated_reference_metadata
                ? file(params.curated_reference_metadata)
                : file("${projectDir}/assets/empty_curated_reference_metadata.tsv")
            BUILD_CURATED_TAXID_MAP(file(params.reference_fasta), curated_metadata)
            MAKE_CURATED_DB(file(params.reference_fasta), BUILD_CURATED_TAXID_MAP.out.taxid_map, 'curated')
            BLASTN_ASVS(ch_asvs, MAKE_CURATED_DB.out.db)
        }

        // Produce the fallback BLAST channel (curated-primary modes only).
        def ch_fallback_blast
        if (mode == 'broad_blast' && curated_check) {
            def broad_db = (params.blast_db ?: params.fallback_blast_db) as String
            if (!broad_db) {
                error "reference_mode 'broad_blast' requires --blast_db (a local BLAST database prefix)."
            }
            BLASTN_EXTERNAL_FALLBACK(ch_asvs, broad_db, params.blast_db_label as String, 'false')
            ch_fallback_blast = BLASTN_EXTERNAL_FALLBACK.out.blast
            fallback_label = params.blast_db_label as String
            use_fallback = true
        } else if (mode == 'remote_fallback') {
            if (!params.enable_ncbi_remote_fallback) {
                log.warn "reference_mode 'remote_fallback' selected but --enable_ncbi_remote_fallback is false; running curated panel only."
            } else {
                def remote_db = (params.remote_blast_db ?: 'core_nt') as String
                BLASTN_EXTERNAL_FALLBACK(ch_asvs, remote_db, 'ncbi_nt_remote', 'true')
                ch_fallback_blast = BLASTN_EXTERNAL_FALLBACK.out.blast
                fallback_label = 'ncbi_nt_remote'
                use_fallback = true
            }
        } else if (mode == 'bold_aware') {
            if (!params.bold_fasta) {
                error "reference_mode 'bold_aware' requires --bold_fasta (a BOLD-derived COI FASTA)."
            }
            if (params.bold_mode == 'api_query') {
                log.warn "bold_mode 'api_query' (live BOLD API) is not implemented; using the reproducible local BOLD FASTA (--bold_fasta)."
            }
            def bold_metadata = params.bold_taxonomy
                ? file(params.bold_taxonomy)
                : file("${projectDir}/assets/empty_curated_reference_metadata.tsv")
            MAKE_BOLD_DB(file(params.bold_fasta), bold_metadata, 'bold_coi')
            BLASTN_BOLD_ASVS(ch_asvs, MAKE_BOLD_DB.out.db)
            ch_fallback_blast = BLASTN_BOLD_ASVS.out.blast
            fallback_label = 'bold_coi'
            use_fallback = true
        }

        // --- Assignment ----------------------------------------------------
        if (broad_only) {
            // Mode B, broad DB only (legacy nt_only): external primary, no curated panel.
            def broad_db = (params.blast_db ?: params.fallback_blast_db) as String
            if (!broad_db) {
                error "reference_mode 'broad_blast' requires --blast_db (a local BLAST database prefix)."
            }
            BLASTN_EXTERNAL_PRIMARY(ch_asvs, broad_db, params.blast_db_label as String, 'false')
            ASSIGN_TAXONOMY_EXTERNAL_DB(BLASTN_EXTERNAL_PRIMARY.out.blast)
            ch_assignments = ASSIGN_TAXONOMY_EXTERNAL_DB.out.assignments
        } else if (use_fallback) {
            // Curated primary + one fallback DB (broad / nt-remote / BOLD).
            def ch_primary = BLASTN_ASVS.out.blast
                .map { meta, counts, blast -> tuple(keyOf(meta), meta, counts, blast) }
            def ch_fb = ch_fallback_blast
                .map { meta, _counts, blast -> tuple(keyOf(meta), blast) }
            def ch_joined = ch_primary
                .join(ch_fb)
                .map { _key, meta, counts, primary_blast, fallback_blast -> tuple(meta, counts, primary_blast, fallback_blast) }
                .combine(BUILD_CURATED_TAXID_MAP.out.taxonomy_table)
                .map { meta, counts, primary_blast, fallback_blast, taxonomy_table -> tuple(meta, counts, primary_blast, fallback_blast, taxonomy_table) }
            ASSIGN_TAXONOMY_WITH_FALLBACK(ch_joined, fallback_label)
            ch_assignments = ASSIGN_TAXONOMY_WITH_FALLBACK.out.assignments
        } else {
            // Mode A (or a fallback-disabled degrade): curated panel only.
            def ch_primary_with_taxonomy = BLASTN_ASVS.out.blast
                .combine(BUILD_CURATED_TAXID_MAP.out.taxonomy_table)
                .map { meta, counts, blast, taxonomy_table -> tuple(meta, counts, blast, taxonomy_table) }
            ASSIGN_TAXONOMY_PRIMARY(ch_primary_with_taxonomy)
            ch_assignments = ASSIGN_TAXONOMY_PRIMARY.out.assignments
        }
    }

    emit:
    assignments = ch_assignments
}
