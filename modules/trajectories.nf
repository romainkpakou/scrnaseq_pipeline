// =============================================================================
// MODULE NEXTFLOW : TRAJECTOIRES CELLULAIRES
// =============================================================================
//
// Encapsule bin/06_trajectories.R.
//
// Ce module est doublement conditionnel : il requiert que l'analyse de
// trajectoire soit activee, mais aussi que l'annotation ait ete executee, car
// la selection du compartiment d'interet s'appuie sur les identites
// biologiques et non sur les numeros de cluster.
//
// L'objet emis est un SingleCellExperiment et non un objet Seurat : la
// conversion est imposee par slingshot, qui appartient a l'ecosysteme
// Bioconductor.
// =============================================================================

process TRAJECTORIES {

    tag "Trajectories"

    publishDir "${params.output_path}/figures", mode: 'copy', pattern: '*.png'
    publishDir "${params.output_path}/objects", mode: 'copy', pattern: '*.rds'
    publishDir "${params.output_path}/tables",  mode: 'copy', pattern: '*.csv'
    publishDir "${params.output_path}/logs",    mode: 'copy', pattern: '*.txt'

    input:
    path input_rds        // objet Seurat annote
    path config
    path script

    output:
    path "sce_slingshot.rds",      emit: rds
    path "dynamic_genes.csv",      emit: dynamic
    path "trajectory_summary.csv", emit: summary
    path "*.png",                  emit: figures
    path "trajectories_log.txt",   emit: log

    script:
    """
    Rscript ${script} \\
        --input_rds ${input_rds} \\
        --config ${config} \\
        --output_rds sce_slingshot.rds \\
        --output_dynamic dynamic_genes.csv \\
        --output_summary trajectory_summary.csv \\
        --output_figures . \\
        --log trajectories_log.txt
    """
}
