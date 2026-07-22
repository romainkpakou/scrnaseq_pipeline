#!/usr/bin/env Rscript

# ==========================================================================
# Module 6 : Trajectoires cellulaires (Slingshot) + genes dynamiques
# ==========================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(yaml)
  library(dplyr)
  library(ggplot2)
  library(slingshot)
  library(SingleCellExperiment)
})

option_list <- list(
  make_option("--input_rds", type = "character", help = "RDS annote"),
  make_option("--config", type = "character"),
  make_option("--output_rds", type = "character", default = "sce_slingshot.rds"),
  make_option("--output_dynamic", type = "character", default = "dynamic_genes.csv"),
  make_option("--output_summary", type = "character", default = "trajectory_summary.csv"),
  make_option("--output_figures", type = "character", default = "."),
  make_option("--log", type = "character", default = "trajectories_log.txt")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg <- yaml::read_yaml(opt$config)
tr  <- cfg$trajectories

log_lines <- c("=== Module Trajectoires cellulaires ===",
               paste("Date :", Sys.time()))

obj <- readRDS(opt$input_rds)

# --- Subset des types cellulaires d'interet ---
types_keep <- tr$cell_types_to_include
obj_sub <- subset(obj, subset = cell_type %in% types_keep)
n_sub <- ncol(obj_sub)
log_lines <- c(log_lines,
               paste("Types inclus :", paste(types_keep, collapse = ", ")),
               paste("Cellules dans le subset :", n_sub))

# --- Conversion en SingleCellExperiment (preserve PCA + UMAP) ---
sce <- as.SingleCellExperiment(obj_sub)

# --- Slingshot sur l'espace UMAP, ancre sur le type de depart ---
sce <- slingshot(sce,
                 clusterLabels = "cell_type",
                 reducedDim = "UMAP",
                 start.clus = tr$start_cluster)

n_lineages <- length(slingLineages(sce))
log_lines <- c(log_lines, paste("Lignages identifies :", n_lineages))
for (ln in names(slingLineages(sce))) {
  path <- paste(slingLineages(sce)[[ln]], collapse = " -> ")
  log_lines <- c(log_lines, paste0("  ", ln, " : ", path))
}

# Pseudotemps (curve 1)
pt <- slingPseudotime(sce)[, 1]
log_lines <- c(log_lines,
               paste("Pseudotemps : min", round(min(pt, na.rm = TRUE), 2),
                     "| max", round(max(pt, na.rm = TRUE), 2),
                     "| median", round(median(pt, na.rm = TRUE), 2)))

# --- Figures trajectoire ---
umap_df <- as.data.frame(reducedDim(sce, "UMAP"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cell_type <- colData(sce)$cell_type
umap_df$pseudotime <- pt

# Extraire la courbe principale
curve1 <- slingCurves(sce)[[1]]
curve_df <- as.data.frame(curve1$s[curve1$ord, ])
colnames(curve_df)[1:2] <- c("UMAP_1", "UMAP_2")

# 21 : pseudotemps
p_pt <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = pseudotime)) +
  geom_point(size = 1) +
  scale_color_viridis_c(option = "plasma", na.value = "grey80") +
  geom_path(data = curve_df, aes(UMAP_1, UMAP_2),
            color = "black", linewidth = 1, inherit.aes = FALSE) +
  ggtitle("Trajectoire T - Pseudotemps Slingshot") +
  theme_minimal()
ggsave(file.path(opt$output_figures, "21_umap_subset_pseudotime.png"),
       p_pt, width = 10, height = 7, dpi = 300)

# 22 : trajectoire coloree par type
p_traj <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = cell_type)) +
  geom_point(size = 1) +
  geom_path(data = curve_df, aes(UMAP_1, UMAP_2),
            color = "black", linewidth = 1.2, inherit.aes = FALSE) +
  ggtitle("Trajectoire T pure - Slingshot",
          subtitle = paste(n_lineages, "lignage(s) T identifie(s)")) +
  labs(color = "Type cellulaire") +
  theme_minimal()
ggsave(file.path(opt$output_figures, "22_umap_subset_trajectories.png"),
       p_traj, width = 10, height = 7, dpi = 300)

# --- Genes dynamiques : correlation de Spearman avec le pseudotemps ---
logcounts_mat <- as.matrix(logcounts(sce))
valid_cells <- !is.na(pt)
pt_valid <- pt[valid_cells]
mat_valid <- logcounts_mat[, valid_cells]

# Filtrer genes exprimes dans > min_pct des cellules
pct_expr <- rowMeans(mat_valid > 0)
genes_keep <- pct_expr > tr$min_pct_expressed
mat_valid <- mat_valid[genes_keep, ]
log_lines <- c(log_lines,
               paste("Genes testes (>", tr$min_pct_expressed * 100, "% cellules) :",
                     nrow(mat_valid)))

# Correlation de Spearman
correlations <- apply(mat_valid, 1, function(g) {
  cor(g, pt_valid, method = "spearman")
})
cor_df <- data.frame(gene = names(correlations),
                     spearman = correlations) %>%
  arrange(desc(spearman))
write.csv(cor_df, opt$output_dynamic, row.names = FALSE)

top_pos <- head(cor_df, 3)$gene
top_neg <- tail(cor_df, 3)$gene
log_lines <- c(log_lines,
               paste("Top 3 positifs :", paste(top_pos, collapse = ", ")),
               paste("Top 3 negatifs :", paste(rev(top_neg), collapse = ", ")))

# --- Heatmap des top genes dynamiques ---
n_top <- tr$n_top_dynamic_genes %/% 2
top_genes <- c(head(cor_df$gene, n_top), tail(cor_df$gene, n_top))
cell_order <- order(pt_valid)
heat_mat <- mat_valid[top_genes, cell_order]

png(file.path(opt$output_figures, "24_heatmap_dynamic_genes.png"),
    width = 2400, height = 1600, res = 200)
heatmap(heat_mat, Rowv = NA, Colv = NA, scale = "row",
        col = colorRampPalette(c("navy", "white", "firebrick"))(100),
        labCol = "", main = "Top 30 genes dynamiques - Trajectoire T",
        ylab = "Genes (top 15 positifs + top 15 negatifs)",
        xlab = "Cellules (ordonnees par pseudotemps)")
dev.off()

# --- Sauvegardes ---
summary_df <- data.frame(
  metric = c("n_cells_subset", "n_lineages", "pt_min", "pt_max", "pt_median",
             "n_genes_tested"),
  value = c(n_sub, n_lineages, round(min(pt, na.rm = TRUE), 2),
            round(max(pt, na.rm = TRUE), 2), round(median(pt, na.rm = TRUE), 2),
            nrow(mat_valid))
)
write.csv(summary_df, opt$output_summary, row.names = FALSE)

saveRDS(sce, file = opt$output_rds)
writeLines(log_lines, opt$log)

cat("Module Trajectoires termine.\n")
cat(paste(log_lines, collapse = "\n"), "\n")
