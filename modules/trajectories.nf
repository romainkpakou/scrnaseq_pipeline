process TRAJECTORIES {
    tag "Trajectories"

    publishDir "${params.output_path}/figures", mode: 'copy', pattern: '*.png'
    publishDir "${params.output_path}/objects", mode: 'copy', pattern: '*.rds'
    publishDir "${params.output_path}/tables",  mode: 'copy', pattern: '*.csv'
    publishDir "${params.output_path}/logs",    mode: 'copy', pattern: '*.txt'

    input:
    path input_rds
    path config

    output:
    path "sce_slingshot.rds",       emit: rds
    path "dynamic_genes.csv",       emit: dynamic
    path "trajectory_summary.csv",  emit: summary
    path "*.png",                   emit: figures
    path "trajectories_log.txt",    emit: log

    script:
    """
    Rscript ${projectDir}/bin/06_trajectories.R \\
        --input_rds ${input_rds} \\
        --config ${config} \\
        --output_rds sce_slingshot.rds \\
        --output_dynamic dynamic_genes.csv \\
        --output_summary trajectory_summary.csv \\
        --output_figures . \\
        --log trajectories_log.txt
    """
}
