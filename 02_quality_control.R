# ============================================================
# 02_quality_control.R
# Calculate QC metrics, visualize distributions, and filter
# low-quality cells based on config thresholds.
# ============================================================

source("R/utils.R")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

step_banner("02", "Quality Control")

config <- load_config()
seurat_obj <- load_seurat("01_raw", config)

# ── Calculate QC metrics ─────────────────────────────────────
mt_pattern <- config$qc$mt_pattern %||% "^MT-"

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj,
                                                   pattern = mt_pattern)

seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj,
                                                   pattern = "^RP[SL]")

message("QC metrics calculated.")
message("  Median nFeature_RNA : ", median(seurat_obj$nFeature_RNA))
message("  Median nCount_RNA   : ", median(seurat_obj$nCount_RNA))
message("  Median percent.mt   : ", round(median(seurat_obj$percent.mt), 2), "%")

# ── Plot QC metrics BEFORE filtering ────────────────────────
p1 <- VlnPlot(seurat_obj,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0.1) &
  theme_scrna() &
  theme(axis.text.x = element_blank())

p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA",
                     feature2 = "nFeature_RNA") +
  geom_hline(yintercept = c(config$qc$min_features,
                             config$qc$max_features),
             linetype = "dashed", colour = "red") +
  theme_scrna() +
  labs(title = "UMI counts vs Genes per cell")

p3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA",
                     feature2 = "percent.mt") +
  geom_hline(yintercept = config$qc$max_mt_percent,
             linetype = "dashed", colour = "red") +
  theme_scrna() +
  labs(title = "UMI counts vs Mitochondrial %")

qc_before <- p1 / (p2 | p3) +
  plot_annotation(title   = "QC Metrics — Before Filtering",
                  subtitle = paste0("n = ", ncol(seurat_obj), " cells"),
                  theme    = theme_scrna())

save_plot(qc_before, "qc_before_filtering", config, width = 12, height = 10)

# ── Apply QC filters ─────────────────────────────────────────
cells_before <- ncol(seurat_obj)

seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA >= config$qc$min_features &
           nFeature_RNA <= config$qc$max_features &
           nCount_RNA   >= config$qc$min_counts   &
           nCount_RNA   <= config$qc$max_counts   &
           percent.mt   <= config$qc$max_mt_percent
)

cells_after  <- ncol(seurat_obj)
cells_removed <- cells_before - cells_after

message("\n── QC Filtering results ────────────────────────")
message("  Cells before : ", cells_before)
message("  Cells removed: ", cells_removed,
        " (", round(cells_removed / cells_before * 100, 1), "%)")
message("  Cells after  : ", cells_after)
message("────────────────────────────────────────────────\n")

# ── Plot QC metrics AFTER filtering ─────────────────────────
qc_after <- VlnPlot(seurat_obj,
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    ncol = 3, pt.size = 0.1) &
  theme_scrna() &
  theme(axis.text.x = element_blank())
qc_after <- qc_after +
  plot_annotation(title    = "QC Metrics — After Filtering",
                  subtitle = paste0("n = ", cells_after, " cells retained"),
                  theme    = theme_scrna())

save_plot(qc_after, "qc_after_filtering", config, width = 10, height = 5)

# ── Save QC summary table ────────────────────────────────────
qc_summary <- data.frame(
  Metric         = c("Cells before QC", "Cells removed", "Cells after QC",
                     "Removal %", "Median nFeature", "Median nCount",
                     "Median MT%"),
  Value          = c(cells_before, cells_removed, cells_after,
                     round(cells_removed / cells_before * 100, 1),
                     round(median(seurat_obj$nFeature_RNA), 0),
                     round(median(seurat_obj$nCount_RNA), 0),
                     round(median(seurat_obj$percent.mt), 2))
)
save_table(qc_summary, "qc_summary", config)

# ── Save filtered Seurat object ──────────────────────────────
save_seurat(seurat_obj, "02_qc_filtered", config)

message("✓ Step 02 complete — QC filtering done.")
