#!/usr/bin/env Rscript

# ==============================================================================
# MODULE 1 : CONTROLE QUALITE ET NORMALISATION
# ==============================================================================
#
# ROLE
#   Premiere etape du pipeline scrnaseq_pipeline. Ce script transforme une
#   matrice de comptage brute 10X Genomics en un objet Seurat filtre, normalise
#   et pret pour la reduction de dimensionnalite.
#
#   Trois operations successives :
#     1. Filtrage des artefacts techniques (debris, doublets, cellules mortes)
#     2. Normalisation de la profondeur de sequencage entre cellules
#     3. Selection des genes portant le signal biologique discriminant
#
# ENTREES
#   --input           Dossier contenant la matrice 10X (matrix.mtx, barcodes.tsv
#                     et genes.tsv ou features.tsv selon la version du pipeline
#                     Cell Ranger)
#   --config          Fichier YAML de configuration de l'analyse
#
# SORTIES
#   --output_rds      Objet Seurat filtre et normalise (.rds)
#   --output_figures  Repertoire de destination des figures :
#                       01_qc_violin_before.png  distributions avant filtrage
#                       02_qc_scatter.png        correlations entre metriques
#                       03_qc_violin_after.png   distributions apres filtrage
#                       04_hvg_plot.png          genes hautement variables
#   --log             Fichier texte des metriques cles (lu par le rapport final)
#
# EXEMPLE D'APPEL
#   Rscript 01_qc_normalization.R \
#       --input data/raw/pbmc3k/hg19 \
#       --config params/example_pbmc.yaml \
#       --output_rds pbmc_qc_normalized.rds \
#       --output_figures . \
#       --log qc_log.txt
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
# suppressPackageStartupMessages evite que les bannieres de demarrage des
# packages ne polluent le flux de sortie capture par Nextflow. Les messages
# utiles sont ecrits explicitement dans le log en fin de script.

suppressPackageStartupMessages({
  library(optparse)    # interface en ligne de commande
  library(Seurat)      # ecosysteme single-cell
  library(yaml)        # lecture de la configuration
  library(ggplot2)     # figures
  library(patchwork)   # composition de figures multiples
})


# ------------------------------------------------------------------------------
# 1. INTERFACE EN LIGNE DE COMMANDE
# ------------------------------------------------------------------------------
# Le script est concu pour etre appele par un process Nextflow, jamais en dur
# depuis une session interactive. Tous les chemins sont donc parametrables :
# c'est ce qui permet au meme script de traiter n'importe quel jeu de donnees
# sans modification de son code source.

option_list <- list(
  make_option("--input",  type = "character",
              help = "Repertoire de la matrice 10X"),
  make_option("--config", type = "character",
              help = "Fichier YAML de configuration"),
  make_option("--output_rds", type = "character",
              default = "pbmc_qc_normalized.rds",
              help = "Nom de l'objet Seurat en sortie"),
  make_option("--output_figures", type = "character", default = ".",
              help = "Repertoire de destination des figures"),
  make_option("--log", type = "character", default = "qc_log.txt",
              help = "Fichier de log des metriques")
)
opt <- parse_args(OptionParser(option_list = option_list))


# ------------------------------------------------------------------------------
# 2. CHARGEMENT DE LA CONFIGURATION
# ------------------------------------------------------------------------------
# Aucun seuil n'est code en dur dans ce script : tous proviennent du YAML.
# Cette separation stricte entre code et parametres est le principe qui rend
# le pipeline reutilisable et son execution tracable, puisque le fichier de
# configuration constitue a lui seul l'enregistrement complet de l'analyse.

cfg  <- yaml::read_yaml(opt$config)
qc   <- cfg$qc              # seuils de filtrage
norm <- cfg$normalization   # parametres de normalisation

# Le log accumule les metriques au fil de l'execution. Il est ecrit en une
# seule fois a la fin, et sera relu par le module de rapport pour composer
# dynamiquement son texte.
log_lines <- c("=== Module QC + Normalisation ===",
               paste("Date :", Sys.time()),
               paste("Input :", opt$input))


# ------------------------------------------------------------------------------
# 3. LECTURE DE LA MATRICE ET CREATION DE L'OBJET SEURAT
# ------------------------------------------------------------------------------
# Read10X lit le triplet de fichiers produit par Cell Ranger et reconstruit
# une matrice creuse genes x cellules. Le format creux est indispensable :
# une matrice single-cell est composee a plus de 90 % de zeros, et sa version
# dense saturerait la memoire.

counts <- Read10X(data.dir = opt$input)

# Deux filtres preliminaires sont appliques des la creation de l'objet :
#   min.cells    : un gene detecte dans trop peu de cellules n'apporte aucun
#                  pouvoir discriminant et ajoute du bruit aux analyses
#   min.features : un code-barre avec trop peu de genes correspond
#                  vraisemblablement a une gouttelette vide ou a du debris
obj <- CreateSeuratObject(
  counts       = counts,
  project      = cfg$project_name,
  min.cells    = qc$min_cells_per_gene,
  min.features = qc$min_features
)

# --- Calcul de la fraction mitochondriale ---
# Rationnel biologique : une cellule dont la membrane est compromise libere
# son ARN cytoplasmique dans le milieu, mais conserve ses mitochondries plus
# longtemps. La proportion relative de transcrits mitochondriaux augmente
# donc mecaniquement, ce qui en fait un marqueur fiable de souffrance ou
# d'apoptose cellulaire.
#
# Le motif de reconnaissance des genes mitochondriaux depend de la convention
# de nomenclature de l'espece : majuscules chez l'humain (MT-), minuscules
# chez la souris (mt-). Cette bascule est ce qui rend le module utilisable
# au-dela du seul dataset humain.
mt_pattern <- if (cfg$species == "mouse") "^mt-" else "^MT-"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)

n_before <- ncol(obj)
log_lines <- c(log_lines, paste("Cellules avant filtrage :", n_before))


# ------------------------------------------------------------------------------
# 4. FIGURES DE DIAGNOSTIC AVANT FILTRAGE
# ------------------------------------------------------------------------------
# Ces figures ne servent pas a decider automatiquement des seuils : elles
# permettent a l'analyste de verifier a posteriori que les seuils du YAML
# etaient adaptes a la distribution reelle des donnees. Un pipeline qui
# filtre sans donner a voir ce qu'il filtre n'est pas auditable.

# Distributions des trois metriques. pt.size = 0 masque les points individuels
# qui, a plusieurs milliers de cellules, formeraient une masse illisible
# recouvrant les violons.
p1 <- VlnPlot(obj,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0)
ggsave(file.path(opt$output_figures, "01_qc_violin_before.png"),
       p1, width = 12, height = 5, dpi = 300)

# Nuages de correlation entre metriques :
#   s1 : nCount vs percent.mt  identifie les cellules stressees
#   s2 : nCount vs nFeature    doit montrer une correlation forte (> 0,95).
#        Une correlation degradee signalerait un probleme technique de
#        sequencage plutot qu'une variabilite biologique.
s1 <- FeatureScatter(obj, "nCount_RNA", "percent.mt")
s2 <- FeatureScatter(obj, "nCount_RNA", "nFeature_RNA")
ggsave(file.path(opt$output_figures, "02_qc_scatter.png"),
       s1 + s2, width = 12, height = 5, dpi = 300)


# ------------------------------------------------------------------------------
# 5. FILTRAGE DES CELLULES
# ------------------------------------------------------------------------------
# Trois criteres combines, chacun ciblant un artefact technique distinct :
#
#   nFeature_RNA > min_features   elimine les gouttelettes vides et les debris
#                                 cellulaires, pauvres en transcrits
#   nFeature_RNA < max_features   elimine les doublets, c'est-a-dire deux
#                                 cellules encapsulees dans une meme gouttelette,
#                                 qui presentent un nombre de genes anormalement
#                                 eleve et un profil transcriptomique chimerique
#   percent.mt   < max_mt_percent elimine les cellules apoptotiques
#
# Ces seuils sont volontairement conservateurs. Un filtrage trop agressif
# supprimerait des populations biologiquement reelles : les monocytes actives
# et les plaquettes, par exemple, presentent naturellement des profils
# atypiques qu'un seuil trop strict confondrait avec des artefacts.

obj <- subset(obj, subset = nFeature_RNA > qc$min_features &
                            nFeature_RNA < qc$max_features &
                            percent.mt   < qc$max_mt_percent)

n_after  <- ncol(obj)
pct_kept <- round(100 * n_after / n_before, 1)
med_mt   <- round(median(obj$percent.mt), 2)

# Le taux de conservation est l'indicateur de sante du filtrage. Un taux
# tres eleve (> 95 %) traduit un jeu de donnees de bonne qualite ; un taux
# faible (< 70 %) doit alerter sur la pertinence des seuils ou sur la qualite
# de la preparation cellulaire.
log_lines <- c(log_lines,
               paste("Cellules apres filtrage :", n_after,
                     paste0("(", pct_kept, "%)")),
               paste("Mediane percent.mt :", med_mt))


# ------------------------------------------------------------------------------
# 6. FIGURE DE CONTROLE APRES FILTRAGE
# ------------------------------------------------------------------------------
# Comparee a la figure 01, cette figure doit montrer des distributions
# resserrees et depourvues des queues extremes. C'est la verification
# visuelle que le filtrage a effectivement fait son travail.

p2 <- VlnPlot(obj,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0)
ggsave(file.path(opt$output_figures, "03_qc_violin_after.png"),
       p2, width = 12, height = 5, dpi = 300)


# ------------------------------------------------------------------------------
# 7. NORMALISATION
# ------------------------------------------------------------------------------
# Probleme resolu : la profondeur de sequencage varie fortement d'une cellule
# a l'autre pour des raisons purement techniques (efficacite de capture,
# amplification). Sans correction, deux cellules du meme type paraitraient
# transcriptomiquement differentes du seul fait de cette variation.
#
# LogNormalize procede en trois temps pour chaque cellule :
#   1. division de l'expression de chaque gene par le total des transcrits
#      de la cellule, ce qui ramene toutes les cellules a une echelle commune
#   2. multiplication par un facteur d'echelle, par convention 10 000, afin
#      de revenir a des valeurs d'ordre de grandeur interpretable
#   3. transformation log(x + 1), qui stabilise la variance des genes
#      fortement exprimes et rapproche la distribution de la normalite
#      attendue par les methodes lineaires en aval

obj <- NormalizeData(obj,
                     normalization.method = norm$method,
                     scale.factor         = norm$scale_factor)


# ------------------------------------------------------------------------------
# 8. SELECTION DES GENES HAUTEMENT VARIABLES
# ------------------------------------------------------------------------------
# Sur environ 20 000 genes detectes, la grande majorite est exprimee de facon
# homogene entre toutes les cellules : ces genes de menage ne portent aucune
# information permettant de distinguer les populations. Les conserver
# reviendrait a noyer le signal biologique dans du bruit.
#
# La methode VST (variance stabilizing transformation) modelise la relation
# attendue entre variance et expression moyenne, puis retient les genes dont
# la variance observee excede nettement cette attente. Elle corrige ainsi le
# biais qui ferait ressortir les genes fortement exprimes du seul fait de
# leur niveau d'expression.

obj <- FindVariableFeatures(obj,
                            selection.method = "vst",
                            nfeatures        = norm$n_hvg)

# Les dix genes les plus variables constituent un controle qualite implicite :
# ils devraient correspondre a des marqueurs de populations attendues dans
# l'echantillon. Leur consignation dans le log permet cette verification.
top10 <- head(VariableFeatures(obj), 10)
log_lines <- c(log_lines, paste("Top 10 HVG :", paste(top10, collapse = ", ")))

# LabelPoints annote les dix premiers genes ; repel = TRUE decale les
# etiquettes pour eviter leur chevauchement dans les zones denses du nuage.
vp <- VariableFeaturePlot(obj)
vp <- LabelPoints(plot = vp, points = top10, repel = TRUE)
ggsave(file.path(opt$output_figures, "04_hvg_plot.png"),
       vp, width = 9, height = 6, dpi = 300)


# ------------------------------------------------------------------------------
# 9. MISE A L'ECHELLE
# ------------------------------------------------------------------------------
# Derniere transformation avant l'ACP : chaque gene est centre (moyenne nulle)
# et reduit (ecart-type unitaire). Sans cette etape, l'ACP serait dominee par
# les quelques genes tres fortement exprimes, dont la variance absolue ecrase
# celle de tous les autres, alors meme qu'ils n'apportent pas necessairement
# le signal le plus discriminant.
#
# Note : la mise a l'echelle porte ici sur l'ensemble des genes plutot que sur
# les seuls HVG. C'est plus couteux en memoire, mais cela permet aux etapes
# ulterieures, notamment les cartes thermiques de marqueurs, de disposer des
# valeurs mises a l'echelle pour n'importe quel gene.

obj <- ScaleData(obj, features = rownames(obj))


# ------------------------------------------------------------------------------
# 10. SAUVEGARDE
# ------------------------------------------------------------------------------
# L'objet .rds est consomme par le module suivant (ACP) via un canal Nextflow.
# Le fichier de log est collecte par le module de rapport en fin de pipeline.

saveRDS(obj, file = opt$output_rds)
writeLines(log_lines, opt$log)

# Sortie console : reprise integrale du log. Nextflow la capture dans
# .command.out, ce qui rend le diagnostic possible sans avoir a fouiller
# le repertoire de travail en cas d'echec.
cat("Module QC termine.\n")
cat(paste(log_lines, collapse = "\n"), "\n")
