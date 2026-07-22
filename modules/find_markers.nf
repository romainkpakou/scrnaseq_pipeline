process FIND_MARKERS {
    tag "Find_markers"

    publishDir "${params.output_path}/figures", mode: 'copy', pattern: '*.png'
    publishDir "${params.output_path}/tables",  mode: 'copy', pattern: '*.csv'
    publishDir "${params.output_path}/logs",    mode: 'copy', pattern: '*.txt'

    input:
    path input_rds
    path config

    output:
    path "all_markers.csv",            emit: markers
    path "top_markers_by_cluster.csv", emit: top_markers
    path "*.png",                      emit: figures
    path "markers_log.txt",            emit: log

    script:
    """
    Rscript ${projectDir}/bin/04_find_markers.R \\
        --input_rds ${input_rds} \\
        --config ${config} \\
        --output_markers all_markers.csv \\
        --output_top top_markers_by_cluster.csv \\
        --output_figures . \\
        --log markers_log.txt
    """
}
