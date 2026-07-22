#!/usr/bin/env Rscript

# ==========================================================================
# Module 4 : Identification des marqueurs differentiels (FindAllMarkers)
# ==========================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(yaml)
  library(dplyr)
  library(ggplot2)
})

option_list <- list(
  make_option("--input_rds", type = "character", help = "RDS issu du clustering"),
  make_option("--config", type = "character"),
  make_option("--output_markers", type = "character", default = "all_markers.csv"),
  make_option("--output_top", type = "character", default = "top_markers_by_cluster.csv"),
  make_option("--output_figures", type = "character", default = "."),
  make_option("--log", type = "character", default = "markers_log.txt")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg <- yaml::read_yaml(opt$config)
mk  <- cfg$markers

log_lines <- c("=== Module Marqueurs differentiels ===",
               paste("Date :", Sys.time()))

obj <- readRDS(opt$input_rds)

# --- FindAllMarkers ---
t0 <- Sys.time()
all_markers <- FindAllMarkers(obj,
                              test.use = mk$test_use,
                              min.pct = mk$min_pct,
                              logfc.threshold = mk$logfc_threshold,
                              only.pos = mk$only_positive,
                              verbose = FALSE)
dt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)

log_lines <- c(log_lines,
               paste("Test :", mk$test_use),
               paste("Temps de calcul :", dt, "min"),
               paste("Marqueurs identifies :", nrow(all_markers)))

write.csv(all_markers, opt$output_markers, row.names = FALSE)

# --- Top N par cluster ---
top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = mk$top_n_per_cluster) %>%
  ungroup()
write.csv(top_markers, opt$output_top, row.names = FALSE)

# --- Figures ---
# DotPlot des top marqueurs (top 3 par cluster pour lisibilite)
top3 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3) %>%
  ungroup()
feats <- unique(top3$gene)

p_dot <- DotPlot(obj, features = feats) + RotatedAxis() +
  theme(axis.text.x = element_text(size = 8))
ggsave(file.path(opt$output_figures, "13_dotplot_canonical.png"),
       p_dot, width = 16, height = 7, dpi = 300)

# VlnPlot du top 1 marqueur de chaque cluster
top1 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 1) %>%
  ungroup()
p_vln <- VlnPlot(obj, features = unique(top1$gene), stack = TRUE, flip = TRUE) +
  NoLegend()
ggsave(file.path(opt$output_figures, "14_vlnplot_top_markers.png"),
       p_vln, width = 10, height = 8, dpi = 300)

# Heatmap des top 10 par cluster
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

writeLines(log_lines, opt$log)

cat("Module Marqueurs termine.\n")
cat(paste(log_lines, collapse = "\n"), "\n")
