#!/usr/bin/env Rscript

# ==============================================================================
# MODULE 5 : ANNOTATION CELLULAIRE AUTOMATISEE
# ==============================================================================
#
# ROLE
#   Attribue automatiquement une identite biologique a chaque cluster en
#   confrontant son profil d'expression a une base de marqueurs canoniques
#   externalisee.
#
#   Ce module remplace une etape traditionnellement manuelle, ou l'analyste
#   inspecte les marqueurs de chaque cluster et reconnait les signatures. Cette
#   automatisation est ce qui rend le pipeline reellement autonome, au prix
#   d'un choix methodologique explicite documente ci-dessous.
#
# METHODE : SCORE DE MODULE
#   Pour chaque type cellulaire de la base, un score est calcule par cellule via
#   AddModuleScore, puis moyenne par cluster. Chaque cluster recoit l'identite
#   dont le score est maximal, sous reserve de depasser un seuil de confiance.
#
#   Pourquoi un score de module plutot qu'une moyenne d'expression brute :
#   AddModuleScore compare l'expression des marqueurs d'un type a celle d'un
#   ensemble de genes temoins de niveau d'expression comparable. Cette
#   correction neutralise le biais qui ferait ressortir un type dont les
#   marqueurs sont simplement tres exprimes partout. Une moyenne brute serait
#   par ailleurs fragile face aux marqueurs partages entre types apparentes.
#
# LIMITE ASSUMEE
#   La methode reste sensible aux types cellulaires dont les signatures se
#   recouvrent largement, typiquement les etats naif et memoire d'une meme
#   lignee. Les scores rapportes dans le log permettent de reperer ces cas : un
#   ecart faible entre le meilleur score et le suivant signale une attribution
#   a verifier manuellement.
#
# ENTREES
#   --input_rds       Objet Seurat avec clusters issu du module 3
#   --config          Fichier YAML de configuration
#   --markers_db      Base de marqueurs canoniques (.yaml, un bloc par type)
#
# SORTIES
#   --output_rds      Objet Seurat avec identites annotees (.rds)
#   --output_summary  Effectifs et frequences par type cellulaire (.csv)
#   --output_figures  Repertoire de destination des figures :
#                       16_umap_annotated.png    UMAP avec etiquettes
#                       17_umap_publication.png  version avec legende laterale
#                       18_dotplot_annotated.png marqueurs par type annote
#   --log             Scores et attributions detailles par cluster
#
# EXEMPLE D'APPEL
#   Rscript 05_annotation.R \
#       --input_rds pbmc_clustered.rds \
#       --config params/example_pbmc.yaml \
#       --markers_db assets/canonical_markers/pbmc_markers_human.yaml \
#       --output_rds pbmc_annotated.rds \
#       --output_summary annotation_summary.csv \
#       --output_figures . \
#       --log annotation_log.txt
#
# DEPENDANCES
#   optparse, Seurat (>= 5.0), yaml, dplyr, ggplot2
#
# AUTEUR
#   Romain KPAKOU
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. CHARGEMENT DES LIBRAIRIES
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(yaml)
  library(dplyr)
  library(ggplot2)
})


# ------------------------------------------------------------------------------
# 1. INTERFACE EN LIGNE DE COMMANDE
# ------------------------------------------------------------------------------

option_list <- list(
  make_option("--input_rds", type = "character",
              help = "Objet Seurat avec clusters issu du module 3"),
  make_option("--config", type = "character",
              help = "Fichier YAML de configuration"),
  make_option("--markers_db", type = "character",
              help = "Base de marqueurs canoniques au format YAML"),
  make_option("--output_rds", type = "character", default = "pbmc_annotated.rds",
              help = "Nom de l'objet Seurat annote en sortie"),
  make_option("--output_summary", type = "character",
              default = "annotation_summary.csv",
              help = "Table de synthese des types cellulaires"),
  make_option("--output_figures", type = "character", default = ".",
              help = "Repertoire de destination des figures"),
  make_option("--log", type = "character", default = "annotation_log.txt",
              help = "Fichier de log des scores et attributions")
)
opt <- parse_args(OptionParser(option_list = option_list))


# ------------------------------------------------------------------------------
# 2. CHARGEMENT DE LA CONFIGURATION, DE LA BASE ET DES DONNEES
# ------------------------------------------------------------------------------
# La base de marqueurs est un fichier distinct du fichier de configuration.
# Cette separation est deliberee : la configuration decrit une execution
# particuliere, tandis que la base de marqueurs constitue une ressource
# scientifique reutilisable, versionnable et extensible independamment. Changer
# d'espece ou de tissu revient a pointer vers une autre base, sans toucher au
# code ni aux parametres d'analyse.

cfg <- yaml::read_yaml(opt$config)
ann <- cfg$annotation
db  <- yaml::read_yaml(opt$markers_db)

log_lines <- c("=== Module Annotation cellulaire ===",
               paste("Date :", Sys.time()),
               paste("Types dans la base :", length(db)))

obj <- readRDS(opt$input_rds)


# ------------------------------------------------------------------------------
# 3. CALCUL DES SCORES DE MODULE PAR TYPE CELLULAIRE
# ------------------------------------------------------------------------------
# Une passe par type cellulaire de la base. Pour chacun, AddModuleScore ajoute
# aux metadonnees une colonne contenant, pour chaque cellule, un score
# d'expression de la signature corrige du fond.
#
# Fonctionnement du calcul : les genes de l'objet sont repartis en classes de
# niveau d'expression comparable ; pour chaque marqueur de la signature, un
# ensemble de genes temoins est tire dans sa classe. Le score final est la
# moyenne d'expression des marqueurs moins la moyenne d'expression des temoins.
# Un score positif signifie que la signature est exprimee au-dela de ce que le
# niveau d'expression general laisserait attendre.
#
# Parametre ctrl : nombre de genes temoins tires par marqueur. Une valeur
# elevee stabilise le score en reduisant la variance du terme de reference.

cell_types <- names(db)
score_cols <- c()

for (i in seq_along(cell_types)) {

  ct <- cell_types[i]

  # Filtrage indispensable : un marqueur absent de la matrice ferait echouer
  # AddModuleScore. L'absence survient notamment lorsque le gene a ete elimine
  # au filtrage qualite, ou lorsque la base cible une autre espece.
  markers_present <- db[[ct]]$markers[db[[ct]]$markers %in% rownames(obj)]

  if (length(markers_present) == 0) {
    # Un type dont aucun marqueur n'est retrouve est ignore plutot que de
    # provoquer un arret : le pipeline reste fonctionnel meme avec une base
    # partiellement inadaptee, et l'avertissement est trace dans le log.
    log_lines <- c(log_lines, paste("ATTENTION : aucun marqueur trouve pour", ct))
    next
  }

  # Les noms de colonnes sont indices numeriquement plutot que nommes d'apres
  # le type cellulaire : les libelles contiennent des espaces et des caracteres
  # speciaux qui poseraient probleme dans les noms de colonnes du data frame de
  # metadonnees. La correspondance est conservee dans le vecteur score_cols.
  col_name <- paste0("score_", i)

  obj <- AddModuleScore(obj,
                        features = list(markers_present),
                        name     = col_name,
                        ctrl     = 20)

  # AddModuleScore suffixe systematiquement le nom fourni par un indice
  # numerique, car il accepte plusieurs signatures simultanement. Une seule
  # signature etant passee ici, le suffixe est toujours 1.
  score_cols[ct] <- paste0(col_name, "1")
}


# ------------------------------------------------------------------------------
# 4. AGREGATION DES SCORES A L'ECHELLE DU CLUSTER
# ------------------------------------------------------------------------------
# L'annotation porte sur les clusters et non sur les cellules individuelles :
# une cellule isolee peut presenter un profil atypique, mais la moyenne sur
# plusieurs centaines de cellules constitue un estimateur robuste de l'identite
# du groupe.
#
# La matrice resultante croise les clusters en lignes et les types cellulaires
# candidats en colonnes.

meta     <- obj@meta.data
clusters <- levels(Idents(obj))

score_matrix <- sapply(score_cols, function(col) {
  tapply(meta[[col]], Idents(obj), mean)
})
colnames(score_matrix) <- names(score_cols)


# ------------------------------------------------------------------------------
# 5. ATTRIBUTION DES IDENTITES
# ------------------------------------------------------------------------------
# Regle de decision : chaque cluster recoit le type cellulaire dont le score
# moyen est le plus eleve, a condition que ce score depasse le seuil de
# confiance defini en configuration.
#
# Justification du garde-fou : les scores de module etant centres autour de
# zero par construction, un score negatif signifie que la signature n'est pas
# exprimee au-dela du fond attendu. Attribuer une identite dans ce cas
# reviendrait a forcer une conclusion que les donnees ne soutiennent pas.
# L'etiquette Unknown est preferable a une annotation fausse : elle signale a
# l'analyste qu'une inspection manuelle s'impose, soit que la base de marqueurs
# ne couvre pas ce type, soit que le cluster correspond a un etat cellulaire
# inattendu.

assignments <- data.frame(cluster = rownames(score_matrix),
                          stringsAsFactors = FALSE)
assignments$best_type  <- NA
assignments$best_score <- NA

for (cl in rownames(score_matrix)) {

  scores     <- score_matrix[cl, ]
  best       <- which.max(scores)
  best_score <- scores[best]

  if (best_score >= ann$score_threshold) {
    assignments$best_type[assignments$cluster == cl] <- names(scores)[best]
  } else {
    assignments$best_type[assignments$cluster == cl] <- "Unknown"
  }

  assignments$best_score[assignments$cluster == cl] <- round(best_score, 3)
}

# Le score retenu est consigne pour chaque cluster. Sa lecture est informative :
# un score eleve traduit une signature nette et specifique, un score faible
# une identite plus incertaine. Les types aux signatures distinctes, comme les
# plaquettes, obtiennent des scores nettement superieurs a ceux des types
# formant un continuum, comme les etats naif et memoire d'une meme lignee.
log_lines <- c(log_lines, "Attribution par cluster :")
for (cl in rownames(score_matrix)) {
  bt      <- assignments$best_type[assignments$cluster == cl]
  bs      <- assignments$best_score[assignments$cluster == cl]
  n_cells <- sum(Idents(obj) == cl)
  log_lines <- c(log_lines,
                 paste0("  Cluster ", cl, " -> ", bt,
                        " (score ", bs, ", ", n_cells, " cellules)"))
}


# ------------------------------------------------------------------------------
# 6. RENOMMAGE DES CLUSTERS
# ------------------------------------------------------------------------------
# Les identites numeriques sont remplacees par les libelles biologiques. Le
# resultat est de plus stocke dans une colonne dediee des metadonnees : les
# identites actives d'un objet Seurat peuvent etre modifiees par des operations
# ulterieures, alors qu'une colonne de metadonnees persiste. Les modules en aval,
# notamment l'analyse de trajectoire, s'appuient sur cette colonne stable.

new_ids <- assignments$best_type
names(new_ids) <- assignments$cluster

obj <- RenameIdents(obj, new_ids)
obj$cell_type <- Idents(obj)


# ------------------------------------------------------------------------------
# 7. TABLE DE SYNTHESE
# ------------------------------------------------------------------------------
# Effectifs et frequences par type, tries par ordre decroissant. Cette table est
# le principal livrable interpretable du module : elle se compare directement
# aux proportions attendues pour le tissu etudie, ce qui constitue un controle
# de vraisemblance de l'annotation dans son ensemble.

summary_df <- as.data.frame(table(obj$cell_type))
colnames(summary_df) <- c("cell_type", "n_cells")
summary_df$pct <- round(100 * summary_df$n_cells / sum(summary_df$n_cells), 2)
summary_df <- summary_df[order(-summary_df$n_cells), ]

write.csv(summary_df, opt$output_summary, row.names = FALSE)


# ------------------------------------------------------------------------------
# 8. FIGURES DE L'ANNOTATION
# ------------------------------------------------------------------------------
# Deux variantes de l'UMAP annote sont produites pour deux usages distincts :
# la version avec etiquettes sur les amas facilite la lecture analytique, la
# version avec legende laterale convient a une figure de publication ou de
# presentation.
#
# Le parametre repel decale les etiquettes qui se chevaucheraient, situation
# frequente lorsque plusieurs populations apparentees occupent une meme region
# de la projection.

p_annot <- DimPlot(obj, reduction = "umap", label = TRUE, label.size = 4,
                   repel = TRUE) + NoLegend() +
  ggtitle("UMAP - types cellulaires annotes")
ggsave(file.path(opt$output_figures, "16_umap_annotated.png"),
       p_annot, width = 10, height = 8, dpi = 300)

p_pub <- DimPlot(obj, reduction = "umap", label = FALSE) +
  ggtitle("UMAP - types cellulaires annotes")
ggsave(file.path(opt$output_figures, "17_umap_publication.png"),
       p_pub, width = 11, height = 8, dpi = 300)

# Controle final de l'annotation : le premier marqueur de chaque type de la base
# est projete sur les clusters desormais nommes. Une annotation correcte se
# traduit par une diagonale nette, chaque marqueur s'exprimant preferentiellement
# sur le type qu'il est cense definir. Toute deviation de cette diagonale
# signale une attribution a reexaminer.
key_markers <- sapply(db, function(x) x$markers[1])
key_markers <- key_markers[key_markers %in% rownames(obj)]

p_dot <- DotPlot(obj, features = unique(key_markers)) + RotatedAxis()
ggsave(file.path(opt$output_figures, "18_dotplot_annotated.png"),
       p_dot, width = 12, height = 7, dpi = 300)


# ------------------------------------------------------------------------------
# 9. SAUVEGARDE
# ------------------------------------------------------------------------------
# L'objet annote alimente le module de trajectoire, qui selectionne les
# populations d'interet par leur nom biologique plutot que par leur numero de
# cluster. Cette dependance explique que l'analyse de trajectoire soit
# conditionnee a l'execution prealable de l'annotation dans le workflow.

saveRDS(obj, file = opt$output_rds)
writeLines(log_lines, opt$log)

cat("Module Annotation termine.\n")
cat(paste(log_lines, collapse = "\n"), "\n")
