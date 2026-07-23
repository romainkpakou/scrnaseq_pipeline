#!/usr/bin/env Rscript

# ==============================================================================
# MODULE 6 : ANALYSE DE TRAJECTOIRE CELLULAIRE ET GENES DYNAMIQUES
# ==============================================================================
#
# ROLE
#   Reconstruit une trajectoire de differentiation a partir des profils
#   transcriptomiques, puis identifie les genes dont l'expression evolue le long
#   de cette trajectoire.
#
# FONDEMENT BIOLOGIQUE
#   Les cellules capturees a un instant donne ne representent pas un etat
#   homogene mais un continuum d'etats de differentiation. Une analyse
#   pseudotemporelle exploite cette heterogeneite pour ordonner les cellules le
#   long d'un axe de progression, reconstruisant in silico une dynamique que
#   l'experience n'a observee qu'a un instant unique.
#
#   Le pseudotemps n'est pas une duree : c'est une coordonnee de progression le
#   long de la trajectoire. Deux cellules de pseudotemps proche sont a un stade
#   de differentiation comparable, sans qu'on puisse en deduire un intervalle
#   temporel reel.
#
# CHOIX METHODOLOGIQUE CENTRAL
#   L'analyse est volontairement restreinte a un sous-ensemble de types
#   cellulaires defini en configuration. Raison : les algorithmes de
#   reconstruction sont non supervises et ne connaissent pas l'ontogenie
#   hematopoietique. Appliques a l'ensemble des populations, ils produisent des
#   lignages qui violent les contraintes developpementales, par exemple une
#   trajectoire reliant des lymphocytes a des monocytes alors que ces lignees
#   divergent des les stades precoces de l'hematopoiese.
#
#   Restreindre l'analyse a un compartiment partageant une origine
#   developpementale commune est ce qui garantit l'interpretabilite du resultat.
#   Cette contrainte biologique est apportee par l'analyste via la
#   configuration : elle ne peut pas etre inferee des donnees.
#
# ENTREES
#   --input_rds        Objet Seurat annote issu du module 5
#   --config           Fichier YAML de configuration
#
# SORTIES
#   --output_rds       Objet SingleCellExperiment avec trajectoire (.rds)
#   --output_dynamic   Correlations gene/pseudotemps, tous genes testes (.csv)
#   --output_summary   Metriques de synthese de la trajectoire (.csv)
#   --output_figures   Repertoire de destination des figures :
#                        21_umap_subset_pseudotime.png    gradient de pseudotemps
#                        22_umap_subset_trajectories.png  courbe ajustee
#                        24_heatmap_dynamic_genes.png     genes dynamiques
#   --log              Metriques cles pour le rapport final
#
# EXEMPLE D'APPEL
#   Rscript 06_trajectories.R \
#       --input_rds pbmc_annotated.rds \
#       --config params/example_pbmc.yaml \
#       --output_rds sce_slingshot.rds \
#       --output_dynamic dynamic_genes.csv \
#       --output_summary trajectory_summary.csv \
#       --output_figures . \
#       --log trajectories_log.txt
#
# DEPENDANCES
#   optparse, Seurat (>= 5.0), yaml, dplyr, ggplot2,
#   slingshot (>= 2.0), SingleCellExperiment
#   Note : slingshot requiert DelayedMatrixStats, dependance parfois absente
#   d'une installation Bioconductor par defaut.
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
  library(slingshot)              # reconstruction de trajectoire
  library(SingleCellExperiment)   # structure de donnees requise par slingshot
})


# ------------------------------------------------------------------------------
# 1. INTERFACE EN LIGNE DE COMMANDE
# ------------------------------------------------------------------------------

option_list <- list(
  make_option("--input_rds", type = "character",
              help = "Objet Seurat annote issu du module 5"),
  make_option("--config", type = "character",
              help = "Fichier YAML de configuration"),
  make_option("--output_rds", type = "character", default = "sce_slingshot.rds",
              help = "Objet SingleCellExperiment avec trajectoire"),
  make_option("--output_dynamic", type = "character", default = "dynamic_genes.csv",
              help = "Table des correlations gene/pseudotemps"),
  make_option("--output_summary", type = "character",
              default = "trajectory_summary.csv",
              help = "Metriques de synthese de la trajectoire"),
  make_option("--output_figures", type = "character", default = ".",
              help = "Repertoire de destination des figures"),
  make_option("--log", type = "character", default = "trajectories_log.txt",
              help = "Fichier de log des metriques")
)
opt <- parse_args(OptionParser(option_list = option_list))


# ------------------------------------------------------------------------------
# 2. CHARGEMENT DE LA CONFIGURATION ET DES DONNEES
# ------------------------------------------------------------------------------

cfg <- yaml::read_yaml(opt$config)
tr  <- cfg$trajectories

log_lines <- c("=== Module Trajectoires cellulaires ===",
               paste("Date :", Sys.time()))

obj <- readRDS(opt$input_rds)


# ------------------------------------------------------------------------------
# 3. SELECTION DU COMPARTIMENT D'INTERET
# ------------------------------------------------------------------------------
# La selection s'appuie sur la colonne cell_type des metadonnees, renseignee par
# le module d'annotation. Utiliser les noms biologiques plutot que les numeros
# de cluster rend la configuration lisible et robuste : le numero attribue a une
# population peut varier d'une execution a l'autre si la resolution change,
# alors que son identite biologique reste stable.
#
# C'est aussi ce qui explique que ce module soit conditionne a l'execution
# prealable de l'annotation dans le workflow Nextflow.

types_keep <- tr$cell_types_to_include
obj_sub    <- subset(obj, subset = cell_type %in% types_keep)
n_sub      <- ncol(obj_sub)

log_lines <- c(log_lines,
               paste("Types inclus :", paste(types_keep, collapse = ", ")),
               paste("Cellules dans le subset :", n_sub))


# ------------------------------------------------------------------------------
# 4. CONVERSION VERS SINGLECELLEXPERIMENT
# ------------------------------------------------------------------------------
# Slingshot appartient a l'ecosysteme Bioconductor et opere sur des objets
# SingleCellExperiment, structure differente de celle de Seurat. La conversion
# preserve les matrices d'expression ainsi que les reductions de dimensionnalite
# calculees en amont, ACP et UMAP, qui sont indispensables a la suite.

sce <- as.SingleCellExperiment(obj_sub)


# ------------------------------------------------------------------------------
# 5. RECONSTRUCTION DE LA TRAJECTOIRE
# ------------------------------------------------------------------------------
# Slingshot procede en deux temps :
#
#   1. Topologie globale. Un arbre couvrant de poids minimal est construit entre
#      les centres des clusters dans l'espace de reduction choisi. Cet arbre
#      determine quels types cellulaires se succedent le long de chaque lignage.
#
#   2. Ajustement local. Des courbes principales sont ajustees le long de chaque
#      branche, et chaque cellule est projetee orthogonalement sur la courbe.
#      Sa position le long de celle-ci constitue son pseudotemps.
#
# Parametres et leur portee :
#
#   clusterLabels  Les identites biologiques servent de noeuds a l'arbre. C'est
#                  ce qui permet de lire directement le lignage sous forme de
#                  succession de types cellulaires.
#
#   reducedDim     L'espace UMAP est retenu ici pour la coherence visuelle entre
#                  la trajectoire et les figures produites en amont. L'espace ACP
#                  serait plus fidele aux distances transcriptomiques reelles,
#                  mais produirait une courbe impossible a superposer aux UMAP
#                  du reste du rapport.
#
#   start.clus     Point de depart biologique de la trajectoire. Ce parametre
#                  est determinant : Slingshot peut identifier la topologie sans
#                  lui, mais pas l'orienter. Sans ancrage, la trajectoire
#                  pourrait etre parcourue a l'envers, du differencie vers le
#                  naif, ce qui inverserait toute l'interpretation des genes
#                  dynamiques.

sce <- slingshot(sce,
                 clusterLabels = "cell_type",
                 reducedDim    = "UMAP",
                 start.clus    = tr$start_cluster)

# Le nombre de lignages est un indicateur a lire attentivement. Un lignage
# unique traduit une progression lineaire. Plusieurs lignages signalent des
# points de branchement, soit reels si le compartiment se differencie vers
# plusieurs destins, soit artefactuels si le sous-ensemble selectionne est
# trop heterogene.
n_lineages <- length(slingLineages(sce))
log_lines  <- c(log_lines, paste("Lignages identifies :", n_lineages))

for (ln in names(slingLineages(sce))) {
  path      <- paste(slingLineages(sce)[[ln]], collapse = " -> ")
  log_lines <- c(log_lines, paste0("  ", ln, " : ", path))
}

# Le pseudotemps de la premiere courbe sert de reference pour la suite. Les
# cellules situees a l'ecart de la courbe recoivent une valeur manquante :
# elles ne sont pas assignables a ce lignage.
pt <- slingPseudotime(sce)[, 1]

log_lines <- c(log_lines,
               paste("Pseudotemps : min", round(min(pt, na.rm = TRUE), 2),
                     "| max", round(max(pt, na.rm = TRUE), 2),
                     "| median", round(median(pt, na.rm = TRUE), 2)))


# ------------------------------------------------------------------------------
# 6. PREPARATION DES DONNEES DE FIGURES
# ------------------------------------------------------------------------------
# Les objets Slingshot ne se tracent pas directement avec ggplot2 : il faut
# extraire d'une part les coordonnees des cellules, d'autre part les points
# constituant la courbe ajustee.
#
# Le champ ord de la courbe donne l'ordre des points le long de celle-ci. Sans
# ce reordonnancement, geom_path relierait les points dans leur ordre de
# stockage et produirait un enchevetrement illisible au lieu d'une courbe.

umap_df <- as.data.frame(reducedDim(sce, "UMAP"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cell_type  <- colData(sce)$cell_type
umap_df$pseudotime <- pt

curve1   <- slingCurves(sce)[[1]]
curve_df <- as.data.frame(curve1$s[curve1$ord, ])
colnames(curve_df)[1:2] <- c("UMAP_1", "UMAP_2")


# ------------------------------------------------------------------------------
# 7. FIGURES DE LA TRAJECTOIRE
# ------------------------------------------------------------------------------
# Deux representations complementaires du meme resultat :
#   - coloration par pseudotemps : montre la progression continue et permet de
#     verifier que le gradient s'oriente bien du naif vers le differencie
#   - coloration par type cellulaire : montre la correspondance entre la
#     progression reconstruite et les identites biologiques annotees
#
# La lecture croisee de ces deux figures constitue la validation de la
# trajectoire : le gradient de pseudotemps doit suivre l'ordre biologique
# attendu des types cellulaires.

# Figure 21 : gradient de pseudotemps. Les cellules non assignees au lignage
# apparaissent en gris via na.value.
p_pt <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = pseudotime)) +
  geom_point(size = 1) +
  scale_color_viridis_c(option = "plasma", na.value = "grey80") +
  geom_path(data = curve_df, aes(UMAP_1, UMAP_2),
            color = "black", linewidth = 1, inherit.aes = FALSE) +
  ggtitle("Trajectoire - Pseudotemps Slingshot") +
  theme_minimal()
ggsave(file.path(opt$output_figures, "21_umap_subset_pseudotime.png"),
       p_pt, width = 10, height = 7, dpi = 300)

# Figure 22 : coloration par type cellulaire. inherit.aes = FALSE est
# indispensable sur la courbe : sans cela, geom_path heriterait de l'esthetique
# de couleur definie au niveau du graphique et chercherait une colonne
# cell_type inexistante dans le tableau de la courbe.
p_traj <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = cell_type)) +
  geom_point(size = 1) +
  geom_path(data = curve_df, aes(UMAP_1, UMAP_2),
            color = "black", linewidth = 1.2, inherit.aes = FALSE) +
  ggtitle("Trajectoire - Slingshot",
          subtitle = paste(n_lineages, "lignage(s) identifie(s)")) +
  labs(color = "Type cellulaire") +
  theme_minimal()
ggsave(file.path(opt$output_figures, "22_umap_subset_trajectories.png"),
       p_traj, width = 10, height = 7, dpi = 300)


# ------------------------------------------------------------------------------
# 8. IDENTIFICATION DES GENES DYNAMIQUES
# ------------------------------------------------------------------------------
# Une fois les cellules ordonnees, il devient possible de reperer les genes dont
# l'expression suit cette progression. Ces genes constituent les effecteurs
# moleculaires de la differentiation observee.
#
# Choix du test : la correlation de Spearman est retenue plutot que celle de
# Pearson car elle est fondee sur les rangs. Elle detecte donc toute relation
# monotone, meme non lineaire, ce qui correspond a la realite biologique ou
# l'induction d'un gene suit rarement une droite. Elle est de plus robuste aux
# valeurs extremes, frequentes dans les donnees d'expression.
#
# Filtre prealable : les genes exprimes dans une trop faible fraction de
# cellules sont ecartes. Sur des donnees creuses, un gene detecte dans quelques
# cellules seulement peut produire une correlation elevee par pur hasard, sans
# aucune signification biologique.

logcounts_mat <- as.matrix(logcounts(sce))

# Les cellules sans pseudotemps assigne sont exclues : elles n'ont pas de
# position sur l'axe de progression et ne peuvent contribuer a la correlation.
valid_cells <- !is.na(pt)
pt_valid    <- pt[valid_cells]
mat_valid   <- logcounts_mat[, valid_cells]

pct_expr    <- rowMeans(mat_valid > 0)
genes_keep  <- pct_expr > tr$min_pct_expressed
mat_valid   <- mat_valid[genes_keep, ]

log_lines <- c(log_lines,
               paste("Genes testes (>", tr$min_pct_expressed * 100, "% cellules) :",
                     nrow(mat_valid)))

correlations <- apply(mat_valid, 1, function(g) {
  cor(g, pt_valid, method = "spearman")
})

# Le tri decroissant place en tete les genes induits le long de la trajectoire
# et en queue ceux dont l'expression decroit. Cette double lecture est
# informative : la differentiation ne consiste pas seulement a activer un
# programme, mais aussi a en eteindre un autre.
cor_df <- data.frame(gene = names(correlations), spearman = correlations) %>%
  arrange(desc(spearman))

write.csv(cor_df, opt$output_dynamic, row.names = FALSE)

top_pos <- head(cor_df, 3)$gene
top_neg <- tail(cor_df, 3)$gene

log_lines <- c(log_lines,
               paste("Top 3 positifs :", paste(top_pos, collapse = ", ")),
               paste("Top 3 negatifs :", paste(rev(top_neg), collapse = ", ")))


# ------------------------------------------------------------------------------
# 9. CARTE THERMIQUE DES GENES DYNAMIQUES
# ------------------------------------------------------------------------------
# Les cellules sont ordonnees par pseudotemps croissant sur l'axe horizontal,
# les genes les plus dynamiques figurent en lignes. La transition
# transcriptomique se lit alors comme un basculement progressif d'intensite le
# long de l'axe des cellules.
#
# Parametres du trace :
#   Rowv, Colv = NA  desactive le regroupement hierarchique automatique. Il est
#                    imperatif de le desactiver : l'ordre des colonnes doit
#                    rester celui du pseudotemps, un reordonnancement
#                    detruirait l'information que la figure est censee montrer.
#   scale = "row"    centre et reduit chaque gene independamment, ce qui permet
#                    de comparer des genes de niveaux d'expression tres
#                    differents sur une meme echelle de couleur.

n_top      <- tr$n_top_dynamic_genes %/% 2
top_genes  <- c(head(cor_df$gene, n_top), tail(cor_df$gene, n_top))
cell_order <- order(pt_valid)
heat_mat   <- mat_valid[top_genes, cell_order]

png(file.path(opt$output_figures, "24_heatmap_dynamic_genes.png"),
    width = 2400, height = 1600, res = 200)
heatmap(heat_mat, Rowv = NA, Colv = NA, scale = "row",
        col    = colorRampPalette(c("navy", "white", "firebrick"))(100),
        labCol = "",
        main   = paste("Top", tr$n_top_dynamic_genes, "genes dynamiques"),
        ylab   = paste("Genes (top", n_top, "positifs +", n_top, "negatifs)"),
        xlab   = "Cellules (ordonnees par pseudotemps)")
dev.off()


# ------------------------------------------------------------------------------
# 10. SAUVEGARDE
# ------------------------------------------------------------------------------
# La table de synthese rassemble les metriques cles sous une forme exploitable
# par le rapport final, qui l'affiche telle quelle sans avoir a analyser le log.

summary_df <- data.frame(
  metric = c("n_cells_subset", "n_lineages", "pt_min", "pt_max", "pt_median",
             "n_genes_tested"),
  value  = c(n_sub, n_lineages,
             round(min(pt, na.rm = TRUE), 2),
             round(max(pt, na.rm = TRUE), 2),
             round(median(pt, na.rm = TRUE), 2),
             nrow(mat_valid))
)
write.csv(summary_df, opt$output_summary, row.names = FALSE)

# L'objet SingleCellExperiment conserve la trajectoire complete, incluant les
# courbes et les pseudotemps de tous les lignages. Il permet de reprendre
# l'analyse en aval sans recalculer, par exemple pour explorer un gene
# particulier.
saveRDS(sce, file = opt$output_rds)
writeLines(log_lines, opt$log)

cat("Module Trajectoires termine.\n")
cat(paste(log_lines, collapse = "\n"), "\n")
