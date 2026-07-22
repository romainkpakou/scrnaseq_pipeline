#!/usr/bin/env Rscript

# ==========================================================================
# Module 3 : Clustering (Louvain) + UMAP
# ==========================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(yaml)
  library(ggplot2)
  library(patchwork)
})

option_list <- list(
  make_option("--input_rds", type = "character", help = "RDS issu du module ACP"),
  make_option("--config", type = "character"),
  make_option("--output_rds", type = "character", default = "pbmc_clustered.rds"),
  make_option("--output_figures", type = "character", default = "."),
  make_option("--log", type = "character", default = "clustering_log.txt")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg  <- yaml::read_yaml(opt$config)
clu  <- cfg$clustering
ump  <- cfg$umap
npcs <- cfg$pca$n_pcs_selected

log_lines <- c("=== Module Clustering + UMAP ===",
               paste("Date :", Sys.time()),
               paste("PC utilisees :", npcs))

obj <- readRDS(opt$input_rds)

# --- Graphe de voisinage + clustering ---
obj <- FindNeighbors(obj, dims = 1:npcs, verbose = FALSE)
obj <- FindClusters(obj, resolution = clu$resolution,
                    algorithm = clu$algorithm, verbose = FALSE)

n_clusters <- length(levels(Idents(obj)))
log_lines <- c(log_lines,
               paste("Resolution :", clu$resolution),
               paste("Nombre de clusters :", n_clusters))

# Distribution des clusters
distrib <- as.data.frame(table(Idents(obj)))
colnames(distrib) <- c("cluster", "n_cellules")
distrib <- distrib[order(-distrib$n_cellules), ]
log_lines <- c(log_lines, "Distribution (decroissante) :",
               paste(distrib$cluster, ":", distrib$n_cellules, collapse = " | "))

# --- UMAP ---
obj <- RunUMAP(obj, dims = 1:npcs,
               n.neighbors = ump$n_neighbors,
               min.dist = ump$min_dist, verbose = FALSE)

# --- Figures UMAP ---
p_clust <- DimPlot(obj, reduction = "umap", label = TRUE, label.size = 6) +
  ggtitle(paste0("UMAP - Clustering Louvain (resolution ", clu$resolution, ")"),
          subtitle = paste(n_clusters, "clusters identifies sur", ncol(obj), "cellules"))
ggsave(file.path(opt$output_figures, "10_umap_clusters.png"),
       p_clust, width = 10, height = 8, dpi = 300)

p_clean <- DimPlot(obj, reduction = "umap", label = FALSE)
ggsave(file.path(opt$output_figures, "11_umap_clusters_clean.png"),
       p_clean, width = 9, height = 7, dpi = 300)

# --- UMAP marqueurs canoniques ---
markers <- ump$canonical_markers
markers <- markers[markers %in% rownames(obj)]
p_feat <- FeaturePlot(obj, features = markers, ncol = 4)
ggsave(file.path(opt$output_figures, "12_umap_canonical_markers.png"),
       p_feat, width = 16, height = 12, dpi = 300)
log_lines <- c(log_lines,
               paste("Marqueurs projetes :", paste(markers, collapse = ", ")))

saveRDS(obj, file = opt$output_rds)
writeLines(log_lines, opt$log)

cat("Module Clustering + UMAP termine.\n")
cat(paste(log_lines, collapse = "\n"), "\n")
