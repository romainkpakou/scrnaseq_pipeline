#!/usr/bin/env Rscript

# ==============================================================================
# MODULE 4 : IDENTIFICATION DES MARQUEURS DIFFERENTIELS
# ==============================================================================
#
# ROLE
#   Determine, pour chaque cluster issu du partitionnement, les genes dont
#   l'expression le distingue significativement de l'ensemble des autres
#   clusters. Ces genes constituent la signature moleculaire du cluster.
#
#   Ce module est le pivot entre le clustering, qui produit des groupes
#   anonymes numerotes, et l'annotation, qui leur attribue une identite
#   biologique. Sans signature moleculaire, un cluster reste un artefact
#   mathematique sans interpretation possible.
#
# ENTREES
#   --input_rds       Objet Seurat avec clusters issu du module 3
#   --config          Fichier YAML de configuration
#
# SORTIES
#   --output_markers  Table complete des marqueurs, tous clusters (.csv)
#   --output_top      Table des N meilleurs marqueurs par cluster (.csv)
#   --output_figures  Repertoire de destination des figures :
#                       13_dotplot_canonical.png   expression par cluster
#                       14_vlnplot_top_markers.png distributions comparees
#                       15_heatmap_top_markers.png signatures completes
#   --log             Metriques cles pour le rapport final
#
# EXEMPLE D'APPEL
#   Rscript 04_find_markers.R \
#       --input_rds pbmc_clustered.rds \
#       --config params/example_pbmc.yaml \
#       --output_markers all_markers.csv \
#       --output_top top_markers_by_cluster.csv \
#       --output_figures . \
#       --log markers_log.txt
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
  library(dplyr)     # manipulation des tables de marqueurs
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
  make_option("--output_markers", type = "character", default = "all_markers.csv",
              help = "Table complete des marqueurs differentiels"),
  make_option("--output_top", type = "character",
              default = "top_markers_by_cluster.csv",
              help = "Table des meilleurs marqueurs par cluster"),
  make_option("--output_figures", type = "character", default = ".",
              help = "Repertoire de destination des figures"),
  make_option("--log", type = "character", default = "markers_log.txt",
              help = "Fichier de log des metriques")
)
opt <- parse_args(OptionParser(option_list = option_list))


# ------------------------------------------------------------------------------
# 2. CHARGEMENT DE LA CONFIGURATION ET DES DONNEES
# ------------------------------------------------------------------------------

cfg <- yaml::read_yaml(opt$config)
mk  <- cfg$markers

log_lines <- c("=== Module Marqueurs differentiels ===",
               paste("Date :", Sys.time()))

obj <- readRDS(opt$input_rds)


# ------------------------------------------------------------------------------
# 3. TEST D'EXPRESSION DIFFERENTIELLE
# ------------------------------------------------------------------------------
# Principe : pour chaque cluster, l'expression de chaque gene est comparee entre
# les cellules du cluster et l'ensemble des cellules des autres clusters. Un
# gene significativement surexprime devient un marqueur candidat.
#
# Parametres et leur justification :
#
#   test_use          Le test de Wilcoxon (rank sum) est le choix par defaut en
#                     single-cell. Non parametrique, il ne suppose aucune
#                     distribution particuliere des donnees, ce qui convient a
#                     l'expression genique dont la distribution est fortement
#                     asymetrique et surchargee en zeros.
#
#   min_pct           Un gene doit etre detecte dans au moins cette fraction de
#                     cellules de l'un des deux groupes compares. Ce filtre
#                     ecarte les genes exprimes dans une poignee de cellules,
#                     qui produiraient des p-values apparemment significatives
#                     sans pertinence biologique.
#
#   logfc_threshold   Seuil minimal de difference d'expression en log2. Il
#                     ecarte les differences statistiquement detectables mais
#                     d'amplitude negligeable, cas frequent lorsque le nombre
#                     de cellules est eleve.
#
#   only_positive     Restreint la recherche aux genes surexprimes. Pour une
#                     annotation cellulaire, un marqueur se definit par sa
#                     presence et non par son absence : conserver les genes
#                     sous-exprimes doublerait le volume de la table sans
#                     apporter d'information exploitable a cette etape.
#
# Avertissement de lecture : les p-values ajustees rapportees sont souvent
# extremement faibles. Cela tient au nombre de cellules, qui confere au test
# une puissance considerable, et non a une certitude biologique proportionnelle.
# La valeur de avg_log2FC et le pourcentage de cellules exprimant le gene
# restent les indicateurs les plus informatifs pour juger de la specificite
# reelle d'un marqueur.

t0 <- Sys.time()

all_markers <- FindAllMarkers(obj,
                              test.use        = mk$test_use,
                              min.pct         = mk$min_pct,
                              logfc.threshold = mk$logfc_threshold,
                              only.pos        = mk$only_positive,
                              verbose         = FALSE)

dt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)

log_lines <- c(log_lines,
               paste("Test :", mk$test_use),
               paste("Temps de calcul :", dt, "min"),
               paste("Marqueurs identifies :", nrow(all_markers)))

# La table complete est conservee : elle constitue la donnee brute de
# l'annotation et permet a l'analyste de rechercher n'importe quel gene
# d'interet sans relancer le calcul.
write.csv(all_markers, opt$output_markers, row.names = FALSE)


# ------------------------------------------------------------------------------
# 4. EXTRACTION DES MEILLEURS MARQUEURS PAR CLUSTER
# ------------------------------------------------------------------------------
# Le classement s'appuie sur avg_log2FC plutot que sur la p-value ajustee.
# Raison : avec plusieurs milliers de cellules, les p-values saturent aux
# limites de la precision numerique et cessent de discriminer les genes entre
# eux. L'amplitude du fold-change reste au contraire directement interpretable
# et reflete la specificite reelle du marqueur.
#
# Note : ce classement peut differer legerement d'une selection manuelle fondee
# sur la connaissance de la litterature. Un marqueur canonique reconnu peut
# presenter un fold-change moderé et se retrouver au-dela des premiers rangs,
# sans que cela remette en cause sa validite biologique.

top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = mk$top_n_per_cluster) %>%
  ungroup()

write.csv(top_markers, opt$output_top, row.names = FALSE)


# ------------------------------------------------------------------------------
# 5. FIGURE : EXPRESSION DES MARQUEURS PAR CLUSTER
# ------------------------------------------------------------------------------
# Le graphique en points croise deux informations pour chaque couple
# gene-cluster : la taille du point indique la proportion de cellules du cluster
# exprimant le gene, la couleur indique le niveau moyen d'expression parmi ces
# cellules.
#
# Cette double lecture est essentielle. Un gene fortement exprime mais dans une
# minorite de cellules, soit un point petit et intense, n'est pas un bon
# marqueur : il signale plutot une sous-population heterogene au sein du
# cluster. Un bon marqueur se reconnait a un point a la fois grand et intense
# sur son cluster, et petit ou pale ailleurs.
#
# La restriction aux trois meilleurs marqueurs par cluster maintient la figure
# lisible : au-dela, l'axe des genes devient trop dense pour etre exploitable.

top3 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3) %>%
  ungroup()

# unique est necessaire : un meme gene peut figurer parmi les meilleurs
# marqueurs de plusieurs clusters, cas frequent pour les populations
# apparentees comme les sous-types lymphocytaires.
feats <- unique(top3$gene)

p_dot <- DotPlot(obj, features = feats) + RotatedAxis() +
  theme(axis.text.x = element_text(size = 8))
ggsave(file.path(opt$output_figures, "13_dotplot_canonical.png"),
       p_dot, width = 16, height = 7, dpi = 300)


# ------------------------------------------------------------------------------
# 6. FIGURE : DISTRIBUTIONS D'EXPRESSION COMPAREES
# ------------------------------------------------------------------------------
# Le graphique en violons montre la distribution complete de l'expression, la ou
# le graphique en points n'en resume que la moyenne. Cette representation revele
# les distributions bimodales, signature d'une population melangee, qu'une
# moyenne masquerait.
#
# Le mode empile avec retournement des axes permet de comparer plusieurs genes
# sur une meme figure compacte, chaque ligne correspondant a un gene et chaque
# colonne a un cluster.

top1 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 1) %>%
  ungroup()

p_vln <- VlnPlot(obj, features = unique(top1$gene),
                 stack = TRUE, flip = TRUE) + NoLegend()
ggsave(file.path(opt$output_figures, "14_vlnplot_top_markers.png"),
       p_vln, width = 10, height = 8, dpi = 300)


# ------------------------------------------------------------------------------
# 7. FIGURE : CARTE THERMIQUE DES SIGNATURES
# ------------------------------------------------------------------------------
# Vue d'ensemble a l'echelle de la cellule individuelle : chaque colonne est une
# cellule, chaque ligne un gene marqueur. Une signature nette se traduit par des
# blocs diagonaux contrastes, chaque bloc correspondant a un cluster exprimant
# specifiquement son groupe de marqueurs.
#
# Des blocs flous ou des chevauchements importants entre clusters adjacents
# signalent des populations transcriptomiquement proches, ce qui constitue une
# information biologique en soi et non necessairement un defaut du clustering.
#
# La taille de police reduite sur l'axe des genes est imposee par leur nombre :
# dix marqueurs multiplies par le nombre de clusters saturent rapidement l'axe.

top10 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()

p_heat <- DoHeatmap(obj, features = top10$gene) + NoLegend() +
  theme(axis.text.y = element_text(size = 5))
ggsave(file.path(opt$output_figures, "15_heatmap_top_markers.png"),
       p_heat, width = 14, height = 12, dpi = 300)

log_lines <- c(log_lines,
               paste("Top", mk$top_n_per_cluster, "marqueurs/cluster sauvegardes"))


# ------------------------------------------------------------------------------
# 8. SAUVEGARDE
# ------------------------------------------------------------------------------
# Ce module ne produit pas d'objet Seurat : il n'apporte aucune modification aux
# donnees elles-memes et se contente de les analyser. Les tables de marqueurs
# sont consommees par le rapport final ; le module d'annotation, quant a lui,
# recalcule ses propres scores a partir de l'objet issu du clustering.

writeLines(log_lines, opt$log)

cat("Module Marqueurs termine.\n")
cat(paste(log_lines, collapse = "\n"), "\n")
