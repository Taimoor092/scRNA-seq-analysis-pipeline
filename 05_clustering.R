# ============================================================
# 05_clustering.R
# Build KNN graph, run Louvain/Leiden clustering, find
# marker genes for each cluster with Wilcoxon test.
# ============================================================

source("R/utils.R")
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

step_banner("05", "Clustering & Marker Genes")

config     <- load_config()
seurat_obj <- load_seurat("04_dim_reduced", config)

pcs_use    <- config$dim_reduction$pcs_use
resolution <- config$clustering$resolution
algorithm  <- config$clustering$algorithm

# ── Build KNN graph ───────────────────────────────────────────
message("Building KNN graph (dims 1:", pcs_use, ")...")

seurat_obj <- FindNeighbors(seurat_obj,
                            dims      = 1:pcs_use,
                            verbose   = FALSE)

# ── Cluster ───────────────────────────────────────────────────
message("Clustering (resolution = ", resolution, ", algorithm = ",
        algorithm, ")...")

seurat_obj <- FindClusters(seurat_obj,
                           resolution = resolution,
                           algorithm  = algorithm,
                           random.seed = config$project$seed,
                           verbose    = FALSE)

n_clusters <- length(unique(Idents(seurat_obj)))
message("Found ", n_clusters, " clusters.")

# ── UMAP coloured by cluster ──────────────────────────────────
umap_clusters <- DimPlot(seurat_obj, reduction = "umap",
                          label = TRUE, label.size = 5,
                          pt.size = 0.5, repel = TRUE) +
  theme_scrna() +
  labs(title    = "UMAP — Clusters",
       subtitle = paste0(n_clusters, " clusters at resolution ", resolution))

save_plot(umap_clusters, "umap_clusters", config, width = 9, height = 8)

# ── Cluster size bar chart ────────────────────────────────────
cluster_sizes <- as.data.frame(table(Idents(seurat_obj)))
colnames(cluster_sizes) <- c("Cluster", "Cells")

bar_plot <- ggplot(cluster_sizes, aes(x = Cluster, y = Cells,
                                       fill = Cluster)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = Cells), vjust = -0.3, size = 3) +
  scale_fill_viridis_d(option = "turbo") +
  theme_scrna() +
  labs(title = "Cells per Cluster", x = "Cluster", y = "Number of cells")

save_plot(bar_plot, "cluster_sizes", config, width = 8, height = 5)
save_table(cluster_sizes, "cluster_sizes", config)

# ── Find marker genes ─────────────────────────────────────────
message("Finding marker genes for each cluster (this may take a minute)...")

all_markers <- FindAllMarkers(
  seurat_obj,
  only.pos         = TRUE,
  min.pct          = config$differential_expression$min_pct,
  logfc.threshold  = config$differential_expression$logfc_threshold,
  test.use         = config$differential_expression$test,
  verbose          = FALSE
)

# Filter significant markers
sig_markers <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC))

message("Found ", nrow(sig_markers), " significant marker genes.")

save_table(all_markers, "all_cluster_markers", config)
save_table(sig_markers, "significant_markers", config)

# ── Top N markers per cluster ──────────────────────────────────
n_top <- config$differential_expression$n_top_markers

top_markers <- sig_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = n_top) %>%
  ungroup()

save_table(top_markers, "top_markers_per_cluster", config)

# ── Heatmap of top 5 markers per cluster ─────────────────────
top5 <- sig_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  pull(gene)

heatmap_plot <- DoHeatmap(seurat_obj,
                           features = unique(top5),
                           size = 3) +
  theme(axis.text.y = element_text(size = 6)) +
  labs(title = "Top 5 Marker Genes per Cluster")

save_plot(heatmap_plot, "heatmap_top_markers", config, width = 12, height = 10)

# ── Save ─────────────────────────────────────────────────────
save_seurat(seurat_obj, "05_clustered", config)

message("✓ Step 05 complete — ", n_clusters, " clusters found, markers identified.")
