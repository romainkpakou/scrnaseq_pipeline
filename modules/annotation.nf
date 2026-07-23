// =============================================================================
// MODULE NEXTFLOW : ANNOTATION CELLULAIRE AUTOMATISEE
// =============================================================================
//
// Encapsule bin/05_annotation.R.
//
// Particularite : ce process recoit une entree supplementaire, la base de
// marqueurs canoniques. Celle-ci doit imperativement etre declaree en entree
// plutot que referencee par un chemin absolu dans la commande : le process
// s'execute dans un repertoire isole ou seuls les fichiers explicitement
// declares sont mis a disposition.
//
// Ce module est conditionnel : le workflow ne l'invoque que si l'annotation
// est activee dans la configuration.
// =============================================================================

process ANNOTATION {

    tag "Annotation"

    publishDir "${params.output_path}/figures", mode: 'copy', pattern: '*.png'
    publishDir "${params.output_path}/objects", mode: 'copy', pattern: '*.rds'
    publishDir "${params.output_path}/tables",  mode: 'copy', pattern: '*.csv'
    publishDir "${params.output_path}/logs",    mode: 'copy', pattern: '*.txt'

    input:
    path input_rds        // objet Seurat avec clusters
    path config
    path markers_db       // base de marqueurs canoniques au format YAML
    path script

    output:
    path "pbmc_annotated.rds",     emit: rds
    path "annotation_summary.csv", emit: summary
    path "*.png",                  emit: figures
    path "annotation_log.txt",     emit: log

    script:
    """
    Rscript ${script} \\
        --input_rds ${input_rds} \\
        --config ${config} \\
        --markers_db ${markers_db} \\
        --output_rds pbmc_annotated.rds \\
        --output_summary annotation_summary.csv \\
        --output_figures . \\
        --log annotation_log.txt
    """
}
