process PCA_SELECTION {
    tag "PCA_selection"

    publishDir "${params.output_path}/figures", mode: 'copy', pattern: '*.png'
    publishDir "${params.output_path}/objects", mode: 'copy', pattern: '*.rds'
    publishDir "${params.output_path}/logs",    mode: 'copy', pattern: '*.txt'

    input:
    path input_rds           // le .rds emis par QC_NORMALIZATION
    path config

    output:
    path "pbmc_pca.rds", emit: rds
    path "*.png",        emit: figures
    path "pca_log.txt",  emit: log

    script:
    """
    Rscript ${projectDir}/bin/02_pca_selection.R \\
        --input_rds ${input_rds} \\
        --config ${config} \\
        --output_rds pbmc_pca.rds \\
        --output_figures . \\
        --log pca_log.txt
    """
}
