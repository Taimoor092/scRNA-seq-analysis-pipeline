# ============================================================
# 03_normalization.R
# Normalize expression data, identify highly variable genes,
# and scale data for downstream analysis.
# Supports LogNormalize and SCTransform methods.
# ============================================================

source("R/utils.R")
library(Seurat)
library(ggplot2)
library(patchwork)

step_banner("03", "Normalization & Feature Selection")

config <- load_config()
seurat_obj <- load_seurat("02_qc_filtered", config)

method      <- config$normalization$method
n_hvg       <- config$normalization$n_variable_features
regress_out <- config$normalization$vars_to_regress

# ── Normalize ────────────────────────────────────────────────
if (method == "SCTransform") {
  message("Using SCTransform normalization...")

  seurat_obj <- SCTransform(
    seurat_obj,
    variable.features.n = n_hvg,
    vars.to.regress     = if (length(regress_out) > 0) regress_out else NULL,
    verbose             = FALSE
  )

  message("SCTransform complete. Default assay set to 'SCT'.")

} else {
  message("Using LogNormalize normalization...")

  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor         = config$normalization$scale_factor
  )

  message("Identifying ", n_hvg, " highly variable genes...")

  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = "vst",
    nfeatures        = n_hvg
  )

  # Scale data (optionally regressing out covariates)
  vars_to_scale <- if (length(regress_out) > 0) regress_out else NULL
  message("Scaling data",
          if (!is.null(vars_to_scale))
            paste0(" (regressing: ", paste(vars_to_scale, collapse = ", "), ")")
          else "...")

  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = vars_to_scale)
}

# ── Plot top variable genes ───────────────────────────────────
top_genes <- head(VariableFeatures(seurat_obj), 20)

hvg_plot <- VariableFeaturePlot(seurat_obj)
hvg_plot <- LabelPoints(plot = hvg_plot, points = top_genes,
                        repel = TRUE, xnudge = 0, ynudge = 0) +
  theme_scrna() +
  labs(title    = "Highly Variable Genes",
       subtitle = paste0("Top ", n_hvg, " HVGs selected (labelling top 20)"))

save_plot(hvg_plot, "highly_variable_genes", config, width = 10, height = 7)

# ── Save HVG list ────────────────────────────────────────────
hvg_df <- data.frame(
  gene = VariableFeatures(seurat_obj),
  rank = seq_along(VariableFeatures(seurat_obj))
)
save_table(hvg_df, "highly_variable_genes", config)

message("Top 10 HVGs: ", paste(head(top_genes, 10), collapse = ", "))

# ── Save normalized Seurat object ────────────────────────────
save_seurat(seurat_obj, "03_normalized", config)

message("✓ Step 03 complete — normalization and scaling done.")
