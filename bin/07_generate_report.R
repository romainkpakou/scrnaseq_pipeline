#!/usr/bin/env Rscript

# ==============================================================================
# MODULE 7 : GENERATION DU RAPPORT HTML FINAL
# ==============================================================================
#
# ROLE
#   Assemble l'ensemble des sorties du pipeline, figures, tables et logs, en un
#   document HTML unique et navigable.
#
#   Ce module est le seul livrable directement consultable par un utilisateur
#   non technique. Sans lui, les resultats resteraient disperses en une
#   vingtaine de figures et une demi-douzaine de tables sans fil conducteur.
#
# PRINCIPE DE CONCEPTION
#   Le template RMarkdown ne contient aucune valeur codee en dur. Tous les
#   chiffres cites dans le texte, effectifs, nombre de clusters, seuils
#   appliques, sont lus dynamiquement depuis les logs, les tables CSV et le
#   fichier de configuration.
#
#   Cette contrainte est ce qui distingue un veritable template d'un rapport
#   fige : le document reste exact quel que soit le jeu de donnees traite, et
#   les sections correspondant a des modules non executes sont automatiquement
#   omises plutot que d'afficher des references vides.
#
# ENTREES
#   --template        Template RMarkdown parametrable
#   --config          Fichier YAML de configuration
#
#   Note : les figures, tables et logs ne sont pas passes en arguments. Ils sont
#   deposes dans le repertoire de travail par Nextflow, et le template les y
#   trouve par leur nom de fichier.
#
# SORTIES
#   --output_html     Rapport HTML autonome, images incluses
#   --log             Metriques de generation
#
# EXEMPLE D'APPEL
#   Rscript 07_generate_report.R \
#       --template assets/report_template.Rmd \
#       --config params/example_pbmc.yaml \
#       --output_html report.html \
#       --log report_log.txt
#
# DEPENDANCES
#   optparse, yaml, rmarkdown, knitr
#   Requiert egalement pandoc (>= 1.12.3) accessible sur le PATH. Une session
#   RStudio le fournit via la variable RSTUDIO_PANDOC, mais une execution en
#   ligne de commande necessite une installation systeme ou l'environnement
#   Conda du pipeline.
#
# AUTEUR
#   Romain KPAKOU
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. CHARGEMENT DES LIBRAIRIES
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(rmarkdown)
})


# ------------------------------------------------------------------------------
# 1. INTERFACE EN LIGNE DE COMMANDE
# ------------------------------------------------------------------------------

option_list <- list(
  make_option("--template", type = "character",
              help = "Template RMarkdown parametrable"),
  make_option("--config", type = "character",
              help = "Fichier YAML de configuration"),
  make_option("--output_html", type = "character", default = "report.html",
              help = "Nom du rapport HTML en sortie"),
  make_option("--log", type = "character", default = "report_log.txt",
              help = "Fichier de log de generation")
)
opt <- parse_args(OptionParser(option_list = option_list))


# ------------------------------------------------------------------------------
# 2. CHARGEMENT DE LA CONFIGURATION
# ------------------------------------------------------------------------------
# Seule la section report est exploitee ici, pour le titre et l'auteur. Le
# template relit lui-meme le fichier complet afin d'en extraire les parametres
# d'analyse qu'il cite dans son texte.

cfg <- yaml::read_yaml(opt$config)
rep <- cfg$report


# ------------------------------------------------------------------------------
# 3. COPIE LOCALE DU TEMPLATE
# ------------------------------------------------------------------------------
# Nextflow met les fichiers d'entree a disposition sous forme de liens
# symboliques pointant vers un stockage en lecture seule. Or le rendu
# RMarkdown ecrit des fichiers intermediaires a cote du document source, ce qui
# echouerait sur un emplacement non inscriptible.
#
# La copie prealable dans le repertoire de travail resout ce point. Le nom
# distinct evite par ailleurs toute collision avec le lien symbolique d'origine.

work_rmd <- "report_working.Rmd"
file.copy(opt$template, work_rmd, overwrite = TRUE)


# ------------------------------------------------------------------------------
# 4. RENDU DU DOCUMENT
# ------------------------------------------------------------------------------
# Les parametres transmis au template alimentent son en-tete YAML. Les valeurs
# de repli garantissent que le rendu aboutit meme si la section report est
# absente ou incomplete dans la configuration.
#
# Le parametre envir isole l'evaluation du document dans un environnement neuf.
# Sans cette precaution, le template pourrait acceder aux variables definies
# dans ce script et masquer une erreur de sa propre logique de lecture.
#
# quiet supprime la sortie verbeuse de pandoc, qui saturerait le journal
# Nextflow sans apporter d'information exploitable.

# Parametres transmis au template, communs aux deux formats de sortie.
render_params <- list(
  config = opt$config,
  title  = if (!is.null(rep$title))  rep$title  else "Analyse single-cell RNA-seq",
  author = if (!is.null(rep$author)) rep$author else "scrnaseq_pipeline"
)

# --- Rendu HTML : livrable principal, toujours produit ---
rmarkdown::render(
  input         = work_rmd,
  output_format = "html_document",
  output_file   = opt$output_html,
  output_dir    = getwd(),
  params        = render_params,
  envir         = new.env(),
  quiet         = TRUE
)

# --- Rendu PDF : format secondaire, active depuis la configuration ---
# Le PDF n'est pas produit systematiquement : sa generation ajoute plusieurs
# dizaines de secondes et requiert une distribution LaTeX, dependance lourde
# que tous les environnements ne fournissent pas.
#
# L'echec du rendu PDF n'interrompt pas le module. Le rapport HTML est deja
# produit a ce stade, et une erreur LaTeX ne doit pas invalider une execution
# de pipeline par ailleurs reussie. L'incident est journalise.
pdf_done <- FALSE

if (isTRUE(rep$pdf)) {
  pdf_file <- sub("\\.html$", ".pdf", opt$output_html)

  tryCatch({
    rmarkdown::render(
      input         = work_rmd,
      output_format = "pdf_document",
      output_file   = pdf_file,
      output_dir    = getwd(),
      params        = render_params,
      envir         = new.env(),
      quiet         = TRUE
    )
    pdf_done <- TRUE
  }, error = function(e) {
    message("AVERTISSEMENT : le rendu PDF a echoue : ", conditionMessage(e))
  })
}


# ------------------------------------------------------------------------------
# 5. JOURNALISATION
# ------------------------------------------------------------------------------
# Le decompte des fichiers presents dans le repertoire de travail permet de
# verifier que Nextflow a bien collecte l'ensemble des sorties amont. Un nombre
# de figures inferieur a l'attendu signale un probleme de chainage des canaux
# plutot qu'une erreur du rendu lui-meme.

n_fig <- length(list.files(".", pattern = "\\.png$"))
n_csv <- length(list.files(".", pattern = "\\.csv$"))

log_lines <- c("=== Module Rapport HTML ===",
               paste("Date :", Sys.time()),
               paste("Figures integrees :", n_fig),
               paste("Tables integrees :", n_csv),
               paste("Sortie HTML :", opt$output_html),
               paste("Sortie PDF :",
                     if (pdf_done) "generee"
                     else if (isTRUE(rep$pdf)) "demandee mais echouee"
                     else "non demandee"))
writeLines(log_lines, opt$log)

cat("Module Rapport termine.\n")
cat(paste(log_lines, collapse = "\n"), "\n")
