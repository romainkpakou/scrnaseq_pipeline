// =============================================================================
// MODULE NEXTFLOW : MARQUEURS DIFFERENTIELS
// =============================================================================
//
// Encapsule bin/04_find_markers.R.
//
// Particularite : ce process ne produit aucun objet Seurat. Il analyse les
// donnees sans les modifier, et emet uniquement des tables et des figures. Les
// deux tables CSV sont emises sur des canaux distincts afin de rester
// referencables individuellement en aval.
// =============================================================================

process FIND_MARKERS {

    tag "Find_markers"

    publishDir "${params.output_path}/figures", mode: 'copy', pattern: '*.png'
    publishDir "${params.output_path}/tables",  mode: 'copy', pattern: '*.csv'
    publishDir "${params.output_path}/logs",    mode: 'copy', pattern: '*.txt'

    input:
    path input_rds        // objet Seurat avec clusters
    path config
    path script

    output:
    path "all_markers.csv",            emit: markers
    path "top_markers_by_cluster.csv", emit: top_markers
    path "*.png",                      emit: figures
    path "markers_log.txt",            emit: log

    script:
    """
    Rscript ${script} \\
        --input_rds ${input_rds} \\
        --config ${config} \\
        --output_markers all_markers.csv \\
        --output_top top_markers_by_cluster.csv \\
        --output_figures . \\
        --log markers_log.txt
    """
}
