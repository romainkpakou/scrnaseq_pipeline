// =============================================================================
// MODULE NEXTFLOW : CONTROLE QUALITE ET NORMALISATION
// =============================================================================
//
// Premiere etape du pipeline. Encapsule bin/01_qc_normalization.R.
//
// Ce process est le seul a consommer les donnees brutes : tous les suivants
// travaillent sur l'objet Seurat qu'il produit.
// =============================================================================

process QC_NORMALIZATION {

    // Etiquette affichee dans le journal d'execution et les rapports de trace.
    tag "QC_normalization"

    // publishDir extrait les sorties du repertoire de travail isole vers le
    // repertoire de resultats visible par l'utilisateur. Sans cette directive,
    // les fichiers resteraient confines dans work/ sous un chemin a hachage.
    // Le tri par motif range chaque type de sortie dans son sous-repertoire.
    publishDir "${params.output_path}/figures", mode: 'copy', pattern: '*.png'
    publishDir "${params.output_path}/objects", mode: 'copy', pattern: '*.rds'
    publishDir "${params.output_path}/logs",    mode: 'copy', pattern: '*.txt'

    input:
    path input_data       // repertoire de la matrice 10X
    path config           // configuration YAML de l'analyse
    path script           // script R, declare en entree pour etre trace

    output:
    // La directive emit nomme chaque canal de sortie, ce qui permet aux
    // process en aval de referencer precisement la sortie qui les interesse
    // plutot que de dependre de l'ordre de declaration.
    path "pbmc_qc_normalized.rds", emit: rds
    path "*.png",                  emit: figures
    path "qc_log.txt",             emit: log

    script:
    """
    Rscript ${script} \\
        --input ${input_data} \\
        --config ${config} \\
        --output_rds pbmc_qc_normalized.rds \\
        --output_figures . \\
        --log qc_log.txt
    """
}
