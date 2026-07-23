// =============================================================================
// MODULE NEXTFLOW : CLUSTERING ET PROJECTION UMAP
// =============================================================================
//
// Encapsule bin/03_clustering_umap.R. Produit l'objet central du pipeline :
// les cellules y portent desormais une identite de cluster et des coordonnees
// UMAP, informations exploitees par tous les modules en aval.
// =============================================================================

process CLUSTERING_UMAP {

    tag "Clustering_UMAP"

    publishDir "${params.output_path}/figures", mode: 'copy', pattern: '*.png'
    publishDir "${params.output_path}/objects", mode: 'copy', pattern: '*.rds'
    publishDir "${params.output_path}/logs",    mode: 'copy', pattern: '*.txt'

    input:
    path input_rds        // objet Seurat avec reduction ACP
    path config
    path script

    output:
    path "pbmc_clustered.rds", emit: rds
    path "*.png",              emit: figures
    path "clustering_log.txt", emit: log

    script:
    """
    Rscript ${script} \\
        --input_rds ${input_rds} \\
        --config ${config} \\
        --output_rds pbmc_clustered.rds \\
        --output_figures . \\
        --log clustering_log.txt
    """
}
