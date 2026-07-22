process ANNOTATION {
    tag "Annotation"

    publishDir "${params.output_path}/figures", mode: 'copy', pattern: '*.png'
    publishDir "${params.output_path}/objects", mode: 'copy', pattern: '*.rds'
    publishDir "${params.output_path}/tables",  mode: 'copy', pattern: '*.csv'
    publishDir "${params.output_path}/logs",    mode: 'copy', pattern: '*.txt'

    input:
    path input_rds
    path config
    path markers_db          // la base de marqueurs, stagee dans le work dir

    output:
    path "pbmc_annotated.rds",       emit: rds
    path "annotation_summary.csv",   emit: summary
    path "*.png",                    emit: figures
    path "annotation_log.txt",       emit: log

    script:
    """
    Rscript ${projectDir}/bin/05_annotation.R \\
        --input_rds ${input_rds} \\
        --config ${config} \\
        --markers_db ${markers_db} \\
        --output_rds pbmc_annotated.rds \\
        --output_summary annotation_summary.csv \\
        --output_figures . \\
        --log annotation_log.txt
    """
}
