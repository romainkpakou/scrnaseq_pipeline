#!/usr/bin/env Rscript

# ==============================================================================
# MODULE 3 : CLUSTERING NON SUPERVISE ET PROJECTION UMAP
# ==============================================================================
#
# ROLE
#   Identifie les populations cellulaires par partitionnement non supervise,
#   puis produit une representation bidimensionnelle exploitable visuellement.
#
#   Deux operations de nature differente sont menees ici, et il importe de ne
#   pas les confondre :
#     - le clustering opere dans l'espace des composantes principales et
#       determine seul l'appartenance des cellules aux populations
#     - l'UMAP ne sert qu'a la visualisation et n'intervient jamais dans la
#       decision d'affectation
#   Une population peut donc apparaitre scindee sur l'UMAP tout en formant un
#   cluster unique, et reciproquement : c'est le clustering qui fait foi.
#
# ENTREES
#   --input_rds       Objet Seurat avec reduction ACP issu du module 2
#   --config          Fichier YAML de configuration
#
# SORTIES
#   --output_rds      Objet Seurat avec clusters et coordonnees UMAP (.rds)
#   --output_figures  Repertoire de destination des figures :
#                       10_umap_clusters.png          UMAP annote des clusters
#                       11_umap_clusters_clean.png    version sans etiquettes
#                       12_umap_canonical_markers.png marqueurs projetes
#   --log             Metriques cles pour le rapport final
#
# EXEMPLE D'APPEL
#   Rscript 03_clustering_umap.R \
#       --input_rds pbmc_pca.rds \
#       --config params/example_pbmc.yaml \
#       --output_rds pbmc_clustered.rds \
#       --output_figures . \
#       --log clustering_log.txt
#
# DEPENDANCES
#   optparse, Seurat (>= 5.0), yaml, ggplot2, patchwork
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
  library(patchwork)
})


# ------------------------------------------------------------------------------
# 1. INTERFACE EN LIGNE DE COMMANDE
# ------------------------------------------------------------------------------

option_list <- list(
  make_option("--input_rds", type = "character",
              help = "Objet Seurat avec ACP issu du module 2"),
  make_option("--config", type = "character",
              help = "Fichier YAML de configuration"),
  make_option("--output_rds", type = "character", default = "pbmc_clustered.rds",
              help = "Nom de l'objet Seurat en sortie"),
  make_option("--output_figures", type = "character", default = ".",
              help = "Repertoire de destination des figures"),
  make_option("--log", type = "character", default = "clustering_log.txt",
              help = "Fichier de log des metriques")
)
opt <- parse_args(OptionParser(option_list = option_list))


# ------------------------------------------------------------------------------
# 2. CHARGEMENT DE LA CONFIGURATION ET DES DONNEES
# ------------------------------------------------------------------------------
# Le nombre de composantes principales est relu depuis la section pca de la
# configuration plutot que depuis l'objet : la configuration reste l'unique
# source de verite pour ce parametre, ce qui evite toute divergence entre le
# nombre de composantes calculees et le nombre effectivement exploite.

cfg  <- yaml::read_yaml(opt$config)
clu  <- cfg$clustering
ump  <- cfg$umap
npcs <- cfg$pca$n_pcs_selected

log_lines <- c("=== Module Clustering + UMAP ===",
               paste("Date :", Sys.time()),
               paste("PC utilisees :", npcs))

obj <- readRDS(opt$input_rds)


# ------------------------------------------------------------------------------
# 3. CONSTRUCTION DU GRAPHE DE VOISINAGE
# ------------------------------------------------------------------------------
# Le clustering repose sur une representation en graphe : chaque cellule devient
# un noeud, relie a ses plus proches voisins dans l'espace des composantes
# principales. Le poids de chaque arete traduit la similarite transcriptomique
# entre les deux cellules qu'elle relie.
#
# Cette representation est mieux adaptee aux donnees single-cell que les
# methodes fondees sur des distances euclidiennes globales, car elle capture la
# structure locale des donnees sans presupposer que les populations aient une
# forme spherique ou une densite homogene.

obj <- FindNeighbors(obj, dims = 1:npcs, verbose = FALSE)


# ------------------------------------------------------------------------------
# 4. PARTITIONNEMENT PAR L'ALGORITHME DE LOUVAIN
# ------------------------------------------------------------------------------
# L'algorithme de Louvain identifie des communautes dans le graphe en optimisant
# la modularite, mesure qui compare la densite de connexions a l'interieur des
# communautes a celle attendue dans un graphe aleatoire de meme degre.
#
# Le parametre de resolution controle la granularite du partitionnement :
#   resolution faible  peu de clusters, larges, risque de fusionner des
#                      populations biologiquement distinctes
#   resolution elevee  nombreux clusters, fins, risque de scinder une population
#                      homogene en sous-groupes sans realite biologique
#
# Il n'existe pas de resolution optimale universelle : le choix se valide par
# la coherence biologique des clusters obtenus, verifiee a l'etape suivante par
# les marqueurs canoniques. La valeur retenue ici provient de la configuration.
#
# Note sur le parametre algorithm : Seurat attend un entier et non une chaine.
# La valeur 1 correspond a l'algorithme de Louvain original, 2 a Louvain avec
# raffinement multi-niveaux, 3 a SLM et 4 a Leiden.

obj <- FindClusters(obj,
                    resolution = clu$resolution,
                    algorithm  = clu$algorithm,
                    verbose    = FALSE)

n_clusters <- length(levels(Idents(obj)))
log_lines <- c(log_lines,
               paste("Resolution :", clu$resolution),
               paste("Nombre de clusters :", n_clusters))

# La distribution des effectifs est consignee par ordre decroissant. Elle
# constitue un controle qualite immediat : un cluster reduit a quelques cellules
# peut correspondre soit a une population rare authentique, comme les plaquettes
# ou les cellules dendritiques, soit a un artefact de sur-partitionnement. Seule
# l'inspection des marqueurs permet de trancher.
distrib <- as.data.frame(table(Idents(obj)))
colnames(distrib) <- c("cluster", "n_cellules")
distrib <- distrib[order(-distrib$n_cellules), ]
log_lines <- c(log_lines, "Distribution (decroissante) :",
               paste(distrib$cluster, ":", distrib$n_cellules, collapse = " | "))


# ------------------------------------------------------------------------------
# 5. PROJECTION UMAP
# ------------------------------------------------------------------------------
# L'UMAP produit une representation bidimensionnelle preservant a la fois la
# structure locale, c'est-a-dire le voisinage immediat de chaque cellule, et
# une partie de la structure globale, soit les distances relatives entre amas.
#
# Deux parametres gouvernent ce compromis :
#   n_neighbors  nombre de voisins consideres pour chaque cellule. Une valeur
#                faible privilegie la structure fine et fragmente la
#                representation ; une valeur elevee lisse et favorise la
#                structure d'ensemble.
#   min_dist     distance minimale imposee entre points dans la projection. Une
#                valeur faible produit des amas compacts et nettement separes ;
#                une valeur elevee etale les points et rend les continuums
#                progressifs plus lisibles.
#
# Avertissement d'interpretation : les distances sur une UMAP ne sont pas
# quantitativement interpretables. Deux amas eloignes ne sont pas
# necessairement plus dissemblables que deux amas proches, et la taille
# apparente d'un amas ne reflete ni son effectif ni son homogeneite.

obj <- RunUMAP(obj,
               dims        = 1:npcs,
               n.neighbors = ump$n_neighbors,
               min.dist    = ump$min_dist,
               verbose     = FALSE)


# ------------------------------------------------------------------------------
# 6. FIGURES DE LA PROJECTION
# ------------------------------------------------------------------------------
# Deux versions sont produites : une annotee, destinee a l'analyse et au
# diagnostic, et une epuree, adaptee a une figure de publication ou de
# presentation ou les etiquettes surchargeraient la lecture.

p_clust <- DimPlot(obj, reduction = "umap", label = TRUE, label.size = 6) +
  ggtitle(paste0("UMAP - Clustering Louvain (resolution ", clu$resolution, ")"),
          subtitle = paste(n_clusters, "clusters identifies sur", ncol(obj), "cellules"))
ggsave(file.path(opt$output_figures, "10_umap_clusters.png"),
       p_clust, width = 10, height = 8, dpi = 300)

p_clean <- DimPlot(obj, reduction = "umap", label = FALSE)
ggsave(file.path(opt$output_figures, "11_umap_clusters_clean.png"),
       p_clean, width = 9, height = 7, dpi = 300)


# ------------------------------------------------------------------------------
# 7. VALIDATION PAR MARQUEURS CANONIQUES
# ------------------------------------------------------------------------------
# Cette figure constitue le controle de coherence biologique du clustering,
# realise avant toute annotation formelle. Chaque panneau projette l'expression
# d'un marqueur etabli sur l'UMAP : si le partitionnement est pertinent, chaque
# marqueur doit se concentrer sur un ou plusieurs clusters bien delimites plutot
# que se disperser uniformement.
#
# Un marqueur diffus sur l'ensemble du nuage signale soit un choix de marqueur
# inadapte a l'espece ou au tissu, soit un clustering qui n'a pas capture la
# structure biologique reelle.
#
# Le filtrage sur les genes effectivement presents dans l'objet est
# indispensable : FeaturePlot echoue si un gene demande est absent de la
# matrice. Ce cas se produit couramment lorsqu'un gene a ete elimine au
# filtrage du module 1, ou lorsque la liste de marqueurs a ete etablie pour une
# autre espece.

markers <- ump$canonical_markers
markers <- markers[markers %in% rownames(obj)]

p_feat <- FeaturePlot(obj, features = markers, ncol = 4)
ggsave(file.path(opt$output_figures, "12_umap_canonical_markers.png"),
       p_feat, width = 16, height = 12, dpi = 300)

log_lines <- c(log_lines,
               paste("Marqueurs projetes :", paste(markers, collapse = ", ")))


# ------------------------------------------------------------------------------
# 8. SAUVEGARDE
# ------------------------------------------------------------------------------
# L'objet transmis au module suivant porte desormais l'identite de cluster de
# chaque cellule ainsi que ses coordonnees UMAP. Ces deux informations sont
# exploitees par tous les modules en aval : les marqueurs differentiels
# comparent les clusters entre eux, l'annotation les renomme, et l'analyse de
# trajectoire se deploie dans l'espace UMAP.

saveRDS(obj, file = opt$output_rds)
writeLines(log_lines, opt$log)

cat("Module Clustering + UMAP termine.\n")
cat(paste(log_lines, collapse = "\n"), "\n")
