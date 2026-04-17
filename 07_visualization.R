# ============================================================
# 07_visualization.R
# Generate a full suite of publication-ready plots from the
# final annotated Seurat object.
# ============================================================

source("R/utils.R")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)
library(viridis)

step_banner("07", "Publication-Ready Visualization")

config     <- load_config()
seurat_obj <- load_seurat("06_annotated", config)

Idents(seurat_obj) <- "cell_type"
cell_types <- sort(unique(seurat_obj$cell_type))
n_types    <- length(cell_types)

# Colour palette — reproducible across all plots
pal <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(n_types),
  cell_types
)

# ── 1. Final UMAP (cell types, no labels) ────────────────────
umap_final <- DimPlot(seurat_obj, reduction = "umap",
                       group.by = "cell_type",
                       cols = pal, pt.size = 0.6,
                       label = FALSE) +
  theme_scrna() +
  labs(title    = "scRNA-seq UMAP",
       subtitle = paste0(ncol(seurat_obj), " cells · ",
                         n_types, " cell types"),
       colour   = "Cell type")

save_plot(umap_final, "figure1_umap_final", config, width = 10, height = 8)

# ── 2. UMAP with labels ───────────────────────────────────────
umap_labelled <- DimPlot(seurat_obj, reduction = "umap",
                          group.by = "cell_type",
                          cols = pal, pt.size = 0.5,
                          label = TRUE, label.size = 4,
                          repel = TRUE) +
  theme_scrna() + NoLegend() +
  labs(title = "UMAP — labelled cell types")

save_plot(umap_labelled, "figure1b_umap_labelled", config, width = 9, height = 8)

# ── 3. Cell type proportion bar chart ────────────────────────
prop_df <- as.data.frame(prop.table(table(seurat_obj$cell_type)) * 100)
colnames(prop_df) <- c("CellType", "Proportion")
prop_df <- prop_df %>% arrange(desc(Proportion))

prop_bar <- ggplot(prop_df, aes(x = reorder(CellType, Proportion),
                                  y = Proportion, fill = CellType)) +
  geom_col(show.legend = FALSE, width = 0.75) +
  scale_fill_manual(values = pal) +
  geom_text(aes(label = paste0(round(Proportion, 1), "%")),
            hjust = -0.1, size = 3.5) +
  coord_flip() +
  expand_limits(y = max(prop_df$Proportion) * 1.12) +
  theme_scrna() +
  labs(title = "Cell type proportions",
       x = NULL, y = "% of total cells")

save_plot(prop_bar, "figure2_cell_proportions", config, width = 9, height = 6)
save_table(prop_df, "cell_type_proportions", config)

# ── 4. Dot plot — canonical marker genes ─────────────────────
# PBMC canonical markers — update if using a different tissue
canonical_markers <- list(
  "CD14+ Monocytes"  = c("CD14", "LYZ", "CST3"),
  "CD16+ Monocytes"  = c("FCGR3A", "MS4A7"),
  "CD4+ T cells"     = c("CD3D", "IL7R", "CCR7"),
  "CD8+ T cells"     = c("CD3D", "CD8A", "GZMB"),
  "NK cells"         = c("GNLY", "NKG7", "KLRD1"),
  "B cells"          = c("MS4A1", "CD79A", "CD19"),
  "Dendritic cells"  = c("FCER1A", "CST7"),
  "Megakaryocytes"   = c("PPBP", "PF4"),
  "Plasmacytoid DCs" = c("LILRA4", "CLEC4C")
)

marker_genes <- unique(unlist(canonical_markers))
present_markers <- marker_genes[marker_genes %in% rownames(seurat_obj)]

if (length(present_markers) > 0) {
  dot_plot <- DotPlot(seurat_obj, features = present_markers,
                       cols = c("lightgrey", "#2166AC"),
                       dot.scale = 6) +
    RotatedAxis() +
    theme_scrna() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
    labs(title    = "Canonical marker gene expression",
         subtitle = "Dot size = % expressing; colour = mean expression",
         x = NULL, y = NULL)

  save_plot(dot_plot, "figure3_dotplot_markers", config, width = 12, height = 7)
}

# ── 5. Feature plots — key genes ─────────────────────────────
feature_genes <- c("CD3D", "CD14", "MS4A1", "GNLY",
                   "FCGR3A", "NKG7", "CD8A", "IL7R")
present_features <- feature_genes[feature_genes %in% rownames(seurat_obj)]

if (length(present_features) > 0) {
  feat_plot <- FeaturePlot(seurat_obj,
                            features  = present_features,
                            reduction = "umap",
                            ncol      = 4,
                            pt.size   = 0.3,
                            cols      = c("lightgrey", "#B2182B")) &
    theme_scrna() &
    theme(axis.text = element_blank(), axis.ticks = element_blank())

  save_plot(feat_plot, "figure4_feature_plots", config,
            width = 16, height = 8)
}

# ── 6. Violin plots — QC metrics by cell type ─────────────────
vln_qc <- VlnPlot(seurat_obj,
                   features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   group.by = "cell_type",
                   cols     = pal,
                   pt.size  = 0,
                   ncol     = 1) &
  theme_scrna() &
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

save_plot(vln_qc, "figure5_qc_by_celltype", config, width = 12, height = 14)

# ── 7. Top marker genes violin (top 2 per type) ───────────────
markers_tbl <- read.csv(file.path(config$output$tables_dir,
                                   "top_markers_per_cluster.csv"),
                         row.names = 1)
top2_genes <- markers_tbl %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 2) %>%
  pull(gene) %>% unique()

top2_genes <- top2_genes[top2_genes %in% rownames(seurat_obj)]
top2_genes <- head(top2_genes, 16)   # cap at 16 for readability

if (length(top2_genes) > 0) {
  vln_markers <- VlnPlot(seurat_obj,
                          features = top2_genes,
                          group.by = "cell_type",
                          cols     = pal,
                          pt.size  = 0,
                          ncol     = 4) &
    theme_scrna() &
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.title.x = element_blank())

  save_plot(vln_markers, "figure6_marker_violins", config,
            width = 16, height = 12)
}

# ── Summary: all figures saved ────────────────────────────────
figs <- list.files(config$output$figures_dir)
message("\n── Figures saved (", length(figs), " files) ──────────────────")
for (f in figs) message("  ", f)
message("─────────────────────────────────────────────────\n")

message("✓ Step 07 complete — all figures generated.")
message("\n🎉 Pipeline complete! Find your results in: results/")
