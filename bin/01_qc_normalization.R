#!/usr/bin/env Rscript

# ==========================================================================
# Module 1 : QC + Normalisation
# Pipeline scrnaseq_pipeline
# ==========================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(yaml)
  library(ggplot2)
  library(patchwork)
})

# --- Arguments CLI ---
option_list <- list(
  make_option("--input", type = "character", help = "Dossier matrice 10X"),
  make_option("--config", type = "character", help = "Fichier YAML de config"),
  make_option("--output_rds", type = "character", default = "pbmc_qc_normalized.rds"),
  make_option("--output_figures", type = "character", default = "."),
  make_option("--log", type = "character", default = "qc_log.txt")
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- Chargement config ---
cfg <- yaml::read_yaml(opt$config)
qc  <- cfg$qc
norm <- cfg$normalization

log_lines <- c("=== Module QC + Normalisation ===",
               paste("Date :", Sys.time()),
               paste("Input :", opt$input))

# --- Lecture 10X + objet Seurat ---
counts <- Read10X(data.dir = opt$input)
obj <- CreateSeuratObject(counts = counts,
                          project = cfg$project_name,
                          min.cells = qc$min_cells_per_gene,
                          min.features = qc$min_features)

# Motif mitochondrial selon espèce
mt_pattern <- if (cfg$species == "mouse") "^mt-" else "^MT-"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)

n_before <- ncol(obj)
log_lines <- c(log_lines, paste("Cellules avant filtrage :", n_before))

# --- Figures QC avant ---
p1 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0)
ggsave(file.path(opt$output_figures, "01_qc_violin_before.png"),
       p1, width = 12, height = 5, dpi = 300)

s1 <- FeatureScatter(obj, "nCount_RNA", "percent.mt")
s2 <- FeatureScatter(obj, "nCount_RNA", "nFeature_RNA")
ggsave(file.path(opt$output_figures, "02_qc_scatter.png"),
       s1 + s2, width = 12, height = 5, dpi = 300)

# --- Filtrage (seuils depuis le YAML) ---
obj <- subset(obj, subset = nFeature_RNA > qc$min_features &
                            nFeature_RNA < qc$max_features &
                            percent.mt < qc$max_mt_percent)

n_after <- ncol(obj)
pct_kept <- round(100 * n_after / n_before, 1)
med_mt <- round(median(obj$percent.mt), 2)
log_lines <- c(log_lines,
               paste("Cellules apres filtrage :", n_after,
                     paste0("(", pct_kept, "%)")),
               paste("Mediane percent.mt :", med_mt))

# --- Figure QC apres ---
p2 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0)
ggsave(file.path(opt$output_figures, "03_qc_violin_after.png"),
       p2, width = 12, height = 5, dpi = 300)

# --- Normalisation + HVG + scaling ---
obj <- NormalizeData(obj, normalization.method = norm$method,
                     scale.factor = norm$scale_factor)
obj <- FindVariableFeatures(obj, selection.method = "vst",
                            nfeatures = norm$n_hvg)

top10 <- head(VariableFeatures(obj), 10)
log_lines <- c(log_lines, paste("Top 10 HVG :", paste(top10, collapse = ", ")))

vp <- VariableFeaturePlot(obj)
vp <- LabelPoints(plot = vp, points = top10, repel = TRUE)
ggsave(file.path(opt$output_figures, "04_hvg_plot.png"),
       vp, width = 9, height = 6, dpi = 300)

obj <- ScaleData(obj, features = rownames(obj))

# --- Sauvegarde ---
saveRDS(obj, file = opt$output_rds)
writeLines(log_lines, opt$log)

cat("Module QC termine.\n")
cat(paste(log_lines, collapse = "\n"), "\n")
