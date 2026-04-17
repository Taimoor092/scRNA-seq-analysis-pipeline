# ============================================================
# 04_dim_reduction.R
# Run PCA, generate elbow plot to choose dimensionality,
# then compute UMAP (and optionally tSNE).
# ============================================================

source("R/utils.R")
library(Seurat)
library(ggplot2)
library(patchwork)

step_banner("04", "Dimensionality Reduction")

config <- load_config()
seurat_obj <- load_seurat("03_normalized", config)

n_pcs    <- config$dim_reduction$n_pcs
pcs_use  <- config$dim_reduction$pcs_use
run_tsne <- isTRUE(config$dim_reduction$run_tsne)

# ── PCA ──────────────────────────────────────────────────────
message("Running PCA (", n_pcs, " components)...")

seurat_obj <- RunPCA(seurat_obj,
                     npcs    = n_pcs,
                     verbose = FALSE)

# ── Elbow plot — use to choose pcs_use in config ─────────────
elbow <- ElbowPlot(seurat_obj, ndims = n_pcs) +
  geom_vline(xintercept = pcs_use,
             linetype = "dashed", colour = "red", linewidth = 0.8) +
  annotate("text", x = pcs_use + 0.5, y = Inf,
           label = paste0("Using PC 1:", pcs_use),
           hjust = 0, vjust = 1.5, colour = "red", size = 3.5) +
  theme_scrna() +
  labs(title    = "PCA Elbow Plot",
       subtitle = "Red line = number of PCs used for clustering & UMAP")

save_plot(elbow, "pca_elbow_plot", config, width = 8, height = 5)

# ── PCA loadings plot ─────────────────────────────────────────
pca_dim_plot <- DimPlot(seurat_obj, reduction = "pca",
                        dims = c(1, 2)) +
  theme_scrna() +
  labs(title = "PCA — PC1 vs PC2")

pca_heat <- DimHeatmap(seurat_obj, dims = 1:6,
                        cells = 500, balanced = TRUE)

save_plot(pca_dim_plot, "pca_scatter", config, width = 8, height = 7)

# ── UMAP ─────────────────────────────────────────────────────
message("Running UMAP (using ", pcs_use, " PCs)...")

seurat_obj <- RunUMAP(
  seurat_obj,
  dims         = 1:pcs_use,
  min.dist     = config$dim_reduction$umap_min_dist     %||% 0.3,
  n.neighbors  = config$dim_reduction$umap_n_neighbors  %||% 30,
  seed.use     = config$project$seed
)

umap_plot <- DimPlot(seurat_obj, reduction = "umap",
                     pt.size = 0.5) +
  theme_scrna() +
  labs(title    = "UMAP — before clustering",
       subtitle = paste0("n = ", ncol(seurat_obj), " cells, ",
                         pcs_use, " PCs"))

save_plot(umap_plot, "umap_preclustering", config, width = 8, height = 7)

# ── tSNE (optional) ───────────────────────────────────────────
if (run_tsne) {
  message("Running tSNE...")
  seurat_obj <- RunTSNE(seurat_obj, dims = 1:pcs_use,
                        seed.use = config$project$seed)

  tsne_plot <- DimPlot(seurat_obj, reduction = "tsne") +
    theme_scrna() +
    labs(title = "tSNE")

  save_plot(tsne_plot, "tsne_preclustering", config, width = 8, height = 7)
}

# ── Save ─────────────────────────────────────────────────────
save_seurat(seurat_obj, "04_dim_reduced", config)

message("✓ Step 04 complete — PCA and UMAP computed.")
