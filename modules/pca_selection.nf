// =============================================================================
// MODULE NEXTFLOW : ANALYSE EN COMPOSANTES PRINCIPALES
// =============================================================================
//
// Encapsule bin/02_pca_selection.R. Consomme l'objet Seurat normalise produit
// par QC_NORMALIZATION et l'enrichit de la reduction ACP.
//
// Le nombre de figures produites varie selon la configuration : le trace
// JackStraw n'est genere que si l'option correspondante est activee. Le motif
// de sortie *.png accepte donc indifferemment quatre ou cinq fichiers.
// =============================================================================

process PCA_SELECTION {

    tag "PCA_selection"

    publishDir "${params.output_path}/figures", mode: 'copy', pattern: '*.png'
    publishDir "${params.output_path}/objects", mode: 'copy', pattern: '*.rds'
    publishDir "${params.output_path}/logs",    mode: 'copy', pattern: '*.txt'

    input:
    path input_rds        // objet Seurat normalise
    path config
    path script

    output:
    path "pbmc_pca.rds", emit: rds
    path "*.png",        emit: figures
    path "pca_log.txt",  emit: log

    script:
    """
    Rscript ${script} \\
        --input_rds ${input_rds} \\
        --config ${config} \\
        --output_rds pbmc_pca.rds \\
        --output_figures . \\
        --log pca_log.txt
    """
}
