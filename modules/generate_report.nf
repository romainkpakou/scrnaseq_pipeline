// =============================================================================
// MODULE NEXTFLOW : RAPPORT HTML FINAL
// =============================================================================
//
// Encapsule bin/07_generate_report.R.
//
// Particularite structurante : ce process est le point de convergence du
// pipeline. Il recoit l'integralite des figures, tables et logs produits en
// amont, agreges par le workflow au moyen d'operateurs de canaux.
//
// Les entrees figures, tables et logs sont des collections et non des fichiers
// uniques. Le workflow leur applique collect() avant transmission, ce qui
// attend que tous les modules amont soient acheves et depose l'ensemble des
// fichiers dans un meme repertoire de travail. Le template les y retrouve par
// leur nom, ce qui explique qu'aucun chemin ne lui soit transmis en argument.
//
// Le rapport HTML est publie a la racine du repertoire de resultats plutot que
// dans un sous-repertoire : c'est le livrable principal destine a
// l'utilisateur final.
// =============================================================================

process GENERATE_REPORT {

    tag "Generate_report"

    publishDir "${params.output_path}",      mode: 'copy', pattern: '*.html'
    publishDir "${params.output_path}",      mode: 'copy', pattern: '*.pdf'
    publishDir "${params.output_path}/logs", mode: 'copy', pattern: '*.txt'

    input:
    path figures          // toutes les figures produites en amont
    path tables           // toutes les tables CSV
    path logs             // tous les journaux de modules
    path config
    path template         // template RMarkdown parametrable
    // L'en-tete LaTeX doit etre stage dans le repertoire de travail : le
    // template le reference par son nom simple, et LaTeX le cherche a cote
    // du document en cours de compilation.
    path tex_header       // en-tete LaTeX pour le rendu PDF
    path script

    output:
    path "report.html",    emit: report
    // optional true : le PDF peut ne pas etre demande, ou son rendu peut
    // echouer sans que cela invalide l'execution du module.
    path "report.pdf",     emit: pdf, optional: true
    path "report_log.txt", emit: log

    script:
    """
    Rscript ${script} \\
        --template ${template} \\
        --config ${config} \\
        --output_html report.html \\
        --log report_log.txt
    """
}
