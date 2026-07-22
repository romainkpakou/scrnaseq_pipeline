process QC_NORMALIZATION {
    tag "QC_normalization"

    // publishDir copie les sorties hors du dossier work/ vers results/
    publishDir "${params.output_path}/figures", mode: 'copy', pattern: '*.png'
    publishDir "${params.output_path}/objects", mode: 'copy', pattern: '*.rds'
    publishDir "${params.output_path}/logs",    mode: 'copy', pattern: '*.txt'

    input:
    path input_data          // le dossier matrice 10X
    path config              // le YAML

    output:
    path "pbmc_qc_normalized.rds", emit: rds      // emit = nomme le canal de sortie
    path "*.png",                  emit: figures
    path "qc_log.txt",             emit: log

    script:
    """
    Rscript ${projectDir}/bin/01_qc_normalization.R \\
        --input ${input_data} \\
        --config ${config} \\
        --output_rds pbmc_qc_normalized.rds \\
        --output_figures . \\
        --log qc_log.txt
    """
}
