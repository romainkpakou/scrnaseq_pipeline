process CLUSTERING_UMAP {
    tag "Clustering_UMAP"

    publishDir "${params.output_path}/figures", mode: 'copy', pattern: '*.png'
    publishDir "${params.output_path}/objects", mode: 'copy', pattern: '*.rds'
    publishDir "${params.output_path}/logs",    mode: 'copy', pattern: '*.txt'

    input:
    path input_rds
    path config

    output:
    path "pbmc_clustered.rds",  emit: rds
    path "*.png",               emit: figures
    path "clustering_log.txt",  emit: log

    script:
    """
    Rscript ${projectDir}/bin/03_clustering_umap.R \\
        --input_rds ${input_rds} \\
        --config ${config} \\
        --output_rds pbmc_clustered.rds \\
        --output_figures . \\
        --log clustering_log.txt
    """
}
