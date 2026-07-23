#!/usr/bin/env Rscript

# ==============================================================================
# MODULE 2 : ANALYSE EN COMPOSANTES PRINCIPALES ET SELECTION DES DIMENSIONS
# ==============================================================================
#
# ROLE
#   Reduit la dimensionnalite de la matrice d'expression normalisee et determine
#   le nombre de composantes principales a conserver pour les analyses en aval.
#
#   Cette etape est le pivot du pipeline : le clustering et l'UMAP ne travaillent
#   jamais sur les genes bruts, mais sur cet espace reduit. La qualite de l'ACP
#   conditionne donc directement celle de toutes les analyses suivantes.
#
# ENTREES
#   --input_rds       Objet Seurat normalise issu du module 1
#   --config          Fichier YAML de configuration
#
# SORTIES
#   --output_rds      Objet Seurat enrichi de la reduction ACP (.rds)
#   --output_figures  Repertoire de destination des figures :
#                       05_pca_dimplot.png    projection sur PC1 et PC2
#                       06_pca_loadings.png   genes contributeurs des deux axes
#                       07_pca_heatmap.png    structure des neuf premieres PC
#                       08_elbow_plot.png     variance expliquee par composante
#                       09_jackstraw_plot.png significativite statistique
#                                             (genere uniquement si active)
#   --log             Metriques cles pour le rapport final
#
# EXEMPLE D'APPEL
#   Rscript 02_pca_selection.R \
#       --input_rds pbmc_qc_normalized.rds \
#       --config params/example_pbmc.yaml \
#       --output_rds pbmc_pca.rds \
#       --output_figures . \
#       --log pca_log.txt
#
# DEPENDANCES
#   optparse, Seurat (>= 5.0), yaml, ggplot2
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
  library(ggplot2)
})


# ------------------------------------------------------------------------------
# 1. INTERFACE EN LIGNE DE COMMANDE
# ------------------------------------------------------------------------------

option_list <- list(
  make_option("--input_rds", type = "character",
              help = "Objet Seurat normalise issu du module QC"),
  make_option("--config", type = "character",
              help = "Fichier YAML de configuration"),
  make_option("--output_rds", type = "character", default = "pbmc_pca.rds",
              help = "Nom de l'objet Seurat en sortie"),
  make_option("--output_figures", type = "character", default = ".",
              help = "Repertoire de destination des figures"),
  make_option("--log", type = "character", default = "pca_log.txt",
              help = "Fichier de log des metriques")
)
opt <- parse_args(OptionParser(option_list = option_list))


# ------------------------------------------------------------------------------
# 2. CHARGEMENT DE LA CONFIGURATION ET DES DONNEES
# ------------------------------------------------------------------------------

cfg <- yaml::read_yaml(opt$config)
pca <- cfg$pca

log_lines <- c("=== Module ACP + selection PC ===",
               paste("Date :", Sys.time()))

obj <- readRDS(opt$input_rds)


# ------------------------------------------------------------------------------
# 3. CALCUL DE L'ANALYSE EN COMPOSANTES PRINCIPALES
# ------------------------------------------------------------------------------
# Probleme resolu : la matrice d'expression compte plusieurs milliers de genes
# hautement variables pour plusieurs milliers de cellules. Cette dimensionnalite
# rend les algorithmes de clustering inoperants, phenomene connu sous le nom de
# fleau de la dimension : dans un espace de tres grande dimension, les distances
# entre points deviennent toutes comparables et perdent leur pouvoir
# discriminant.
#
# L'ACP construit de nouveaux axes, les composantes principales, chacun etant
# une combinaison lineaire de genes co-regules. Les premieres composantes
# capturent l'essentiel de la variance et donc du signal biologique, tandis que
# les dernieres ne portent plus que du bruit technique.
#
# Note : RunPCA travaille par defaut sur les genes hautement variables
# selectionnes au module precedent, ce qui est le comportement souhaite.

obj <- RunPCA(obj, npcs = pca$n_dims, verbose = FALSE)
log_lines <- c(log_lines, paste("ACP calculee :", pca$n_dims, "composantes"))


# ------------------------------------------------------------------------------
# 4. FIGURES D'INTERPRETATION BIOLOGIQUE DES COMPOSANTES
# ------------------------------------------------------------------------------
# Ces trois figures servent un objectif souvent neglige : verifier que les axes
# de variance dominants correspondent bien a des distinctions biologiques
# reelles, et non a des artefacts techniques comme un effet de lot ou un
# gradient de qualite cellulaire.
#
# Sur un jeu de donnees sain, les premieres composantes doivent s'interpreter :
# separation myeloide contre lymphoide, cellules B contre cellules cytotoxiques,
# et ainsi de suite. Cette interpretabilite est un indicateur fort de la qualite
# du pretraitement.

# Projection des cellules sur les deux premiers axes. Une structure en amas
# visible des cette projection annonce un clustering ulterieur bien separe.
p_dim <- DimPlot(obj, reduction = "pca")
ggsave(file.path(opt$output_figures, "05_pca_dimplot.png"),
       p_dim, width = 7, height = 6, dpi = 300)

# Genes contribuant le plus aux deux premiers axes. Le sens de la contribution,
# positif ou negatif, indique quelle population chaque extremite de l'axe
# represente : c'est la lecture qui permet de nommer biologiquement chaque
# composante.
p_load <- VizDimLoadings(obj, dims = 1:2, reduction = "pca")
ggsave(file.path(opt$output_figures, "06_pca_loadings.png"),
       p_load, width = 9, height = 7, dpi = 300)

# Cartes thermiques des neuf premieres composantes. Chaque colonne represente
# une cellule ordonnee selon son score sur l'axe, chaque ligne un gene
# contributeur. Un motif nettement contraste signe une composante porteuse de
# structure biologique ; un motif brouille signale une composante dominee par
# le bruit.
#
# DimHeatmap dessine directement sur le peripherique graphique actif au lieu de
# retourner un objet ggplot : il faut donc ouvrir et fermer explicitement un
# peripherique png, ggsave ne fonctionnerait pas ici. Le sous-echantillonnage
# a 500 cellules garde la figure lisible sans en alterer la lecture.
png(file.path(opt$output_figures, "07_pca_heatmap.png"),
    width = 2400, height = 2400, res = 300)
DimHeatmap(obj, dims = 1:9, cells = 500, balanced = TRUE)
dev.off()


# ------------------------------------------------------------------------------
# 5. SELECTION DU NOMBRE DE COMPOSANTES : GRAPHIQUE DU COUDE
# ------------------------------------------------------------------------------
# Le choix du nombre de composantes est un arbitrage entre deux risques
# opposes : en retenir trop peu revient a perdre du signal biologique et a
# fusionner des populations distinctes ; en retenir trop revient a injecter du
# bruit technique et a fragmenter artificiellement des populations homogenes.
#
# L'ElbowPlot represente la variance expliquee par chaque composante par ordre
# decroissant. Le point d'inflexion, ou la courbe s'aplatit, marque le seuil
# au-dela duquel chaque composante additionnelle n'apporte plus d'information
# substantielle.
#
# La ligne verticale materialise le choix retenu dans la configuration, ce qui
# permet de verifier d'un coup d'oeil que ce choix tombe bien au niveau du coude.

p_elbow <- ElbowPlot(obj, ndims = pca$n_dims) +
  geom_vline(xintercept = pca$n_pcs_selected, linetype = "dashed", color = "red")
ggsave(file.path(opt$output_figures, "08_elbow_plot.png"),
       p_elbow, width = 8, height = 5, dpi = 300)


# ------------------------------------------------------------------------------
# 6. SELECTION DU NOMBRE DE COMPOSANTES : TEST JACKSTRAW (OPTIONNEL)
# ------------------------------------------------------------------------------
# Le graphique du coude repose sur une appreciation visuelle. La procedure
# JackStraw y ajoute un fondement statistique : elle permute aleatoirement une
# fraction des donnees, recalcule l'ACP, et compare la distribution des p-values
# observees a la distribution uniforme attendue sous hypothese nulle. Une
# composante dont les p-values s'ecartent nettement de la diagonale porte un
# signal statistiquement significatif.
#
# Ce test est desactive par defaut dans la configuration car il est couteux,
# de deux a cinq minutes selon le jeu de donnees, pour une conclusion qui
# confirme generalement celle du coude. Il reste activable lorsque la rigueur
# statistique doit etre documentee, par exemple pour une publication.

if (isTRUE(pca$jackstraw)) {
  obj <- JackStraw(obj, num.replicate = 100, dims = 20)
  obj <- ScoreJackStraw(obj, dims = 1:20)

  p_js <- JackStrawPlot(obj, dims = 1:20)
  ggsave(file.path(opt$output_figures, "09_jackstraw_plot.png"),
         p_js, width = 9, height = 6, dpi = 300)

  log_lines <- c(log_lines, "JackStraw : calcule (20 PC)")
} else {
  log_lines <- c(log_lines, "JackStraw : desactive (config)")
}

log_lines <- c(log_lines,
               paste("PC retenues pour la suite :", pca$n_pcs_selected))


# ------------------------------------------------------------------------------
# 7. SAUVEGARDE
# ------------------------------------------------------------------------------
# L'objet transmis au module suivant contient desormais la reduction ACP dans
# son slot de reductions. Le nombre de composantes a utiliser n'est pas stocke
# dans l'objet : il est relu depuis la configuration par chaque module qui en a
# besoin, ce qui garantit qu'une seule source de verite gouverne ce parametre.

saveRDS(obj, file = opt$output_rds)
writeLines(log_lines, opt$log)

cat("Module ACP termine.\n")
cat(paste(log_lines, collapse = "\n"), "\n")
