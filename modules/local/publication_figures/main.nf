process PUBLICATION_FIGURES {
    label 'process_medium'
    tag 'publication_figures'

    // Curated, self-contained publication figures (Objective 1): vector PDF/EPS to results/figures/.
    publishDir params.publication_figures_dir, mode: 'copy'

    input:
    path host_assignments       // host_assignments.tsv (per-ASV; Fig 2 marker recovery/confidence)
    path host_calls             // host_call_table.tsv (rambo-filtered; ecology)
    path master_endpoint        // bloodmeal_master_endpoint.tsv (metadata for bioclim indices)
    path eco_indices            // ecological_indices.tsv (overall stratum)
    path samplesheet            // sample list + metadata (vector composition over ALL mosquitoes)
    path read_decisions         // collected *.read_decisions.tsv (Fig 2A read lengths)

    output:
    path 'main/*',          emit: main
    path 'supplementary/*', emit: supplementary
    path 'figure_data/*',   emit: tables
    path 'versions.yml',    emit: versions

    when:
    params.enable_publication_figures

    script:
    """
    mkdir -p figure_data main supplementary

    # bioclimatic-zone-stratified indices (same validated Wilson-CI computation as the overall run)
    python3 '${projectDir}/bin/compute_ecological_indices.py' \\
        --host-calls '${host_calls}' \\
        --master-endpoint '${master_endpoint}' \\
        --zone-column '${params.figure_bioclim_column}' \\
        --species-column '${params.eco_index_species_column}' \\
        --date-column '${params.eco_index_date_column}' \\
        --output-tsv figure_data/ecological_indices_bioclim.tsv

    # plot-ready tables (all aggregation happens here, not in the plotting scripts)
    python3 '${projectDir}/bin/figure_data_prep.py' \\
        --host-assignments '${host_assignments}' \\
        --host-calls '${host_calls}' \\
        --read-decisions-glob '*.read_decisions.tsv' \\
        --samplesheet '${samplesheet}' \\
        --outdir figure_data

    # vector host-use ecology indices (Fig 5): zooprophylaxis, Levins' niche breadth, Pianka
    # overlap, Bray-Curtis/Jaccard beta-diversity, bipartite network (connectance, H2' proxy).
    python3 '${projectDir}/bin/compute_host_ecology_indices.py' \\
        --host-calls '${host_calls}' \\
        --samplesheet '${samplesheet}' \\
        --outdir figure_data

    # 5 main figures (incl. GADM bioclimatic-zone map)
    PYTHONPATH='${projectDir}/bin' python3 '${projectDir}/bin/build_main_figures.py' \\
        --figure-data figure_data \\
        --eco-overall '${eco_indices}' \\
        --eco-bioclim figure_data/ecological_indices_bioclim.tsv \\
        --geo-dir '${projectDir}/assets/geo' \\
        --outdir main

    # supplementary figures
    PYTHONPATH='${projectDir}/bin' python3 '${projectDir}/bin/build_supplementary_figures.py' \\
        --figure-data figure_data \\
        --eco-bioclim figure_data/ecological_indices_bioclim.tsv \\
        --outdir supplementary

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      python: "\$(python3 --version | sed 's/Python //')"
      geopandas: "\$(python3 -c 'import geopandas; print(geopandas.__version__)')"
      matplotlib: "\$(python3 -c 'import matplotlib; print(matplotlib.__version__)')"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p main supplementary figure_data
    for f in figure_1_study_system figure_2_metabarcoding_efficacy \\
             figure_3_host_use_by_zone figure_4_mixed_feeding figure_5_host_use_ecology; do
        printf 'stub' > main/\${f}.pdf
    done
    printf 'stub' > supplementary/figure_S1_rarefaction.pdf
    printf 'stub' > supplementary/figure_S2_host_zone_matrix.pdf
    printf 'host\\tn\\n' > figure_data/host_overall.tsv
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
