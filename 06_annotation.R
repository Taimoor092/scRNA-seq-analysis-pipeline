# ============================================================
# 06_annotation.R
# Annotate clusters with cell type labels.
# Two modes:
#   "singler" — automated reference-based annotation (SingleR)
#   "manual"  — user-defined marker gene list
# ============================================================

source("R/utils.R")
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

step_banner("06", "Cell Type Annotation")

config     <- load_config()
seurat_obj <- load_seurat("05_clustered", config)

method <- config$annotation$method %||% "singler"

# ============================================================
# METHOD 1: Automated annotation with SingleR
# ============================================================
if (method == "singler") {

  if (!requireNamespace("SingleR", quietly = TRUE) ||
      !requireNamespace("celldex", quietly = TRUE)) {
    stop("Please install SingleR and celldex:\n",
         "BiocManager::install(c('SingleR', 'celldex'))")
  }

  library(SingleR)
  library(celldex)

  ref_name <- config$annotation$reference %||%
    "HumanPrimaryCellAtlasData"

  message("Loading reference: ", ref_name)
  ref <- do.call(ref_name, list())

  message("Running SingleR automated annotation...")

  # Extract log-normalised counts for SingleR
  expr_mat <- GetAssayData(seurat_obj, layer = "data")

  singler_results <- SingleR(
    test      = expr_mat,
    ref       = ref,
    labels    = ref$label.main
  )

  # Assign per-cell labels to Seurat object
  seurat_obj$singler_label <- singler_results$labels
  seurat_obj$singler_score <- apply(singler_results$scores, 1, max)

  # Assign majority label per cluster
  cluster_labels <- seurat_obj@meta.data %>%
    group_by(seurat_clusters) %>%
    count(singler_label) %>%
    slice_max(n, n = 1) %>%
    select(seurat_clusters, cell_type = singler_label)

  label_map <- setNames(cluster_labels$cell_type,
                        cluster_labels$seurat_clusters)

  seurat_obj$cell_type <- label_map[as.character(seurat_obj$seurat_clusters)]

  # Save annotation summary
  anno_summary <- cluster_labels
  anno_summary$n_cells <- table(seurat_obj$seurat_clusters)[
    as.character(anno_summary$seurat_clusters)]
  save_table(as.data.frame(anno_summary), "cluster_annotation", config)

  # SingleR score heatmap
  score_heatmap <- plotScoreHeatmap(singler_results,
                                    show.pruned = TRUE)
  pdf(file.path(config$output$figures_dir, "singler_score_heatmap.pdf"),
      width = 12, height = 8)
  print(score_heatmap)
  dev.off()
  message("SingleR score heatmap saved.")

# ============================================================
# METHOD 2: Manual annotation
# Edit the marker list below to match your dataset and tissue.
# ============================================================
} else if (method == "manual") {

  message("Using manual cell type annotation...")

  # ── EDIT THIS SECTION for your dataset ──────────────────────
  # PBMC example (adjust cluster numbers from your elbow/UMAP)
  manual_labels <- c(
    "0" = "CD14+ Monocytes",
    "1" = "CD4+ T cells",
    "2" = "CD8+ T cells",
    "3" = "NK cells",
    "4" = "B cells",
    "5" = "CD16+ Monocytes",
    "6" = "Dendritic cells",
    "7" = "Megakaryocytes",
    "8" = "Plasmacytoid DCs"
  )

  seurat_obj$cell_type <- manual_labels[
    as.character(seurat_obj$seurat_clusters)]

  # Replace NAs (unmapped clusters) with cluster number
  seurat_obj$cell_type[is.na(seurat_obj$cell_type)] <- paste0(
    "Cluster_",
    seurat_obj$seurat_clusters[is.na(seurat_obj$cell_type)]
  )

} else {
  stop("Unknown annotation method: ", method,
       ". Use 'singler' or 'manual'.")
}

# ── Set cell_type as the active identity ──────────────────────
Idents(seurat_obj) <- "cell_type"

n_types <- length(unique(seurat_obj$cell_type))
message("Annotated ", n_types, " cell types:")
print(table(seurat_obj$cell_type))

# ── UMAP coloured by cell type ────────────────────────────────
umap_celltypes <- DimPlot(seurat_obj, reduction = "umap",
                           group.by = "cell_type",
                           label = TRUE, label.size = 4,
                           pt.size = 0.5, repel = TRUE) +
  theme_scrna() +
  labs(title    = "UMAP — Cell Types",
       subtitle = paste0(n_types, " cell populations identified"))

save_plot(umap_celltypes, "umap_celltypes", config, width = 10, height = 8)

# ── Side-by-side: clusters vs cell types ─────────────────────
p_clust <- DimPlot(seurat_obj, reduction = "umap",
                   group.by = "seurat_clusters",
                   label = TRUE, pt.size = 0.3) +
  theme_scrna() + labs(title = "Clusters") + NoLegend()

p_types <- DimPlot(seurat_obj, reduction = "umap",
                   group.by = "cell_type",
                   label = TRUE, label.size = 3,
                   pt.size = 0.3, repel = TRUE) +
  theme_scrna() + labs(title = "Cell Types") + NoLegend()

combined <- p_clust | p_types
save_plot(combined, "umap_cluster_vs_celltype", config,
          width = 16, height = 7)

# ── Save annotated object ─────────────────────────────────────
save_seurat(seurat_obj, "06_annotated", config)

message("✓ Step 06 complete — cell types annotated.")
