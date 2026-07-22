#!/usr/bin/env Rscript

# ==========================================================================
# Module 5 : Annotation cellulaire automatisee (AddModuleScore)
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
  make_option("--markers_db", type = "character", help = "YAML base marqueurs canoniques"),
  make_option("--output_rds", type = "character", default = "pbmc_annotated.rds"),
  make_option("--output_summary", type = "character", default = "annotation_summary.csv"),
  make_option("--output_figures", type = "character", default = "."),
  make_option("--log", type = "character", default = "annotation_log.txt")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg  <- yaml::read_yaml(opt$config)
ann  <- cfg$annotation
db   <- yaml::read_yaml(opt$markers_db)

log_lines <- c("=== Module Annotation cellulaire ===",
               paste("Date :", Sys.time()),
               paste("Types dans la base :", length(db)))

obj <- readRDS(opt$input_rds)

# --- Calcul d'un module score par type cellulaire ---
cell_types <- names(db)
score_cols <- c()

for (i in seq_along(cell_types)) {
  ct <- cell_types[i]
  markers_present <- db[[ct]]$markers[db[[ct]]$markers %in% rownames(obj)]
  if (length(markers_present) == 0) {
    log_lines <- c(log_lines, paste("ATTENTION : aucun marqueur trouve pour", ct))
    next
  }
  col_name <- paste0("score_", i)
  obj <- AddModuleScore(obj, features = list(markers_present),
                        name = col_name, ctrl = 20)
  # AddModuleScore ajoute une colonne suffixee "1"
  score_cols[ct] <- paste0(col_name, "1")
}

# --- Score moyen par cluster ---
meta <- obj@meta.data
clusters <- levels(Idents(obj))
score_matrix <- sapply(score_cols, function(col) {
  tapply(meta[[col]], Idents(obj), mean)
})
# score_matrix : lignes = clusters, colonnes = types cellulaires
colnames(score_matrix) <- names(score_cols)

# --- Attribution : argmax par cluster, avec garde-fou seuil ---
assignments <- data.frame(cluster = rownames(score_matrix),
                          stringsAsFactors = FALSE)
assignments$best_type <- NA
assignments$best_score <- NA

for (cl in rownames(score_matrix)) {
  scores <- score_matrix[cl, ]
  best <- which.max(scores)
  best_score <- scores[best]
  if (best_score >= ann$score_threshold) {
    assignments$best_type[assignments$cluster == cl] <- names(scores)[best]
  } else {
    assignments$best_type[assignments$cluster == cl] <- "Unknown"
  }
  assignments$best_score[assignments$cluster == cl] <- round(best_score, 3)
}

# --- Log detaille des attributions ---
log_lines <- c(log_lines, "Attribution par cluster :")
for (cl in rownames(score_matrix)) {
  bt <- assignments$best_type[assignments$cluster == cl]
  bs <- assignments$best_score[assignments$cluster == cl]
  n_cells <- sum(Idents(obj) == cl)
  log_lines <- c(log_lines,
                 paste0("  Cluster ", cl, " -> ", bt,
                        " (score ", bs, ", ", n_cells, " cellules)"))
}

# --- Renommage des clusters ---
new_ids <- assignments$best_type
names(new_ids) <- assignments$cluster
obj <- RenameIdents(obj, new_ids)
obj$cell_type <- Idents(obj)

# --- Tableau de synthese ---
summary_df <- as.data.frame(table(obj$cell_type))
colnames(summary_df) <- c("cell_type", "n_cells")
summary_df$pct <- round(100 * summary_df$n_cells / sum(summary_df$n_cells), 2)
summary_df <- summary_df[order(-summary_df$n_cells), ]
write.csv(summary_df, opt$output_summary, row.names = FALSE)

# --- Figures ---
p_annot <- DimPlot(obj, reduction = "umap", label = TRUE, label.size = 4,
                   repel = TRUE) + NoLegend() +
  ggtitle("UMAP - PBMC 3k (types cellulaires annotes)")
ggsave(file.path(opt$output_figures, "16_umap_annotated.png"),
       p_annot, width = 10, height = 8, dpi = 300)

p_pub <- DimPlot(obj, reduction = "umap", label = FALSE) +
  ggtitle("UMAP - PBMC 3k (types cellulaires annotes)")
ggsave(file.path(opt$output_figures, "17_umap_publication.png"),
       p_pub, width = 11, height = 8, dpi = 300)

# DotPlot annote : 1 marqueur canonique phare par type
key_markers <- sapply(db, function(x) x$markers[1])
key_markers <- key_markers[key_markers %in% rownames(obj)]
p_dot <- DotPlot(obj, features = unique(key_markers)) + RotatedAxis()
ggsave(file.path(opt$output_figures, "18_dotplot_annotated.png"),
       p_dot, width = 12, height = 7, dpi = 300)

saveRDS(obj, file = opt$output_rds)
writeLines(log_lines, opt$log)

cat("Module Annotation termine.\n")
cat(paste(log_lines, collapse = "\n"), "\n")
