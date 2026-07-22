#!/usr/bin/env Rscript

# ==========================================================================
# Module 2 : ACP + selection des composantes principales
# ==========================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(yaml)
  library(ggplot2)
})

option_list <- list(
  make_option("--input_rds", type = "character", help = "RDS issu du module QC"),
  make_option("--config", type = "character"),
  make_option("--output_rds", type = "character", default = "pbmc_pca.rds"),
  make_option("--output_figures", type = "character", default = "."),
  make_option("--log", type = "character", default = "pca_log.txt")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg <- yaml::read_yaml(opt$config)
pca <- cfg$pca

log_lines <- c("=== Module ACP + selection PC ===",
               paste("Date :", Sys.time()))

obj <- readRDS(opt$input_rds)

# --- ACP ---
obj <- RunPCA(obj, npcs = pca$n_dims, verbose = FALSE)
log_lines <- c(log_lines, paste("ACP calculee :", pca$n_dims, "composantes"))

# --- Figures ---
p_dim <- DimPlot(obj, reduction = "pca")
ggsave(file.path(opt$output_figures, "05_pca_dimplot.png"),
       p_dim, width = 7, height = 6, dpi = 300)

p_load <- VizDimLoadings(obj, dims = 1:2, reduction = "pca")
ggsave(file.path(opt$output_figures, "06_pca_loadings.png"),
       p_load, width = 9, height = 7, dpi = 300)

png(file.path(opt$output_figures, "07_pca_heatmap.png"),
    width = 2400, height = 2400, res = 300)
DimHeatmap(obj, dims = 1:9, cells = 500, balanced = TRUE)
dev.off()

p_elbow <- ElbowPlot(obj, ndims = pca$n_dims) +
  geom_vline(xintercept = pca$n_pcs_selected, linetype = "dashed", color = "red")
ggsave(file.path(opt$output_figures, "08_elbow_plot.png"),
       p_elbow, width = 8, height = 5, dpi = 300)

# --- JackStraw (optionnel) ---
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

saveRDS(obj, file = opt$output_rds)
writeLines(log_lines, opt$log)

cat("Module ACP termine.\n")
cat(paste(log_lines, collapse = "\n"), "\n")
