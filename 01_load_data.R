# ============================================================
# 01_load_data.R
# Load raw count matrix into a Seurat object.
# Supports 10x Genomics MEX format and .rds input.
# Demo mode: automatically downloads PBMC 3k dataset.
# ============================================================

source("R/utils.R")
library(Seurat)
library(Matrix)

step_banner("01", "Load Data")

config <- load_config()
setup_dirs(config)
set.seed(config$project$seed)

# ── Demo mode: download PBMC 3k from 10x Genomics ───────────
if (isTRUE(config$project$use_demo_data)) {
  message("Demo mode: downloading PBMC 3k dataset...")

  url <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
  tar_file <- file.path(config$data$input_dir, "pbmc3k.tar.gz")
  raw_dir  <- file.path(config$data$input_dir, "pbmc3k")

  if (!dir.exists(raw_dir)) {
    dir.create(raw_dir, recursive = TRUE)
    download.file(url, destfile = tar_file)
    untar(tar_file, exdir = config$data$input_dir)
    message("PBMC 3k data downloaded and extracted to: ", raw_dir)
  } else {
    message("PBMC 3k data already present — skipping download.")
  }

  data_dir <- file.path(config$data$input_dir,
                        "filtered_gene_bc_matrices", "hg19")

} else {
  data_dir <- config$data$input_dir
}

# ── Load count matrix ────────────────────────────────────────
if (config$data$input_format == "rds") {
  message("Loading Seurat object from .rds file...")
  rds_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0) stop("No .rds file found in: ", data_dir)
  seurat_obj <- readRDS(rds_files[1])

} else {
  message("Loading 10x MEX format from: ", data_dir)
  counts <- Read10X(data.dir = data_dir)

  seurat_obj <- CreateSeuratObject(
    counts       = counts,
    project      = config$project$name,
    min.cells    = 3,    # Gene must be in at least 3 cells
    min.features = config$qc$min_features
  )
}

# ── Quick summary ────────────────────────────────────────────
message("\n── Dataset summary ─────────────────────────────")
message("  Cells  : ", ncol(seurat_obj))
message("  Genes  : ", nrow(seurat_obj))
message("  Assays : ", paste(Assays(seurat_obj), collapse = ", "))
message("────────────────────────────────────────────────\n")

# ── Save raw Seurat object ───────────────────────────────────
save_seurat(seurat_obj, "01_raw", config)

message("✓ Step 01 complete — raw Seurat object saved.")
