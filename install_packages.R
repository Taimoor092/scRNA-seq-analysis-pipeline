# ============================================================
# install_packages.R
# Run this once to install all required packages.
# ============================================================

cat("Installing scRNA-seq pipeline dependencies...\n\n")

# --- CRAN packages ---
cran_packages <- c(
  "Seurat",       # Core scRNA-seq framework
  "SeuratObject", # Seurat data structures
  "ggplot2",      # Visualization
  "dplyr",        # Data manipulation
  "patchwork",    # Combine ggplot2 panels
  "yaml",         # Parse config file
  "Matrix",       # Sparse matrix support
  "scales",       # Plot scales/colours
  "viridis",      # Colour palettes
  "RColorBrewer", # Colour palettes
  "ggrepel",      # Non-overlapping labels
  "cowplot",      # Plot utilities
  "future"        # Parallelisation
)

installed <- rownames(installed.packages())
to_install <- cran_packages[!cran_packages %in% installed]

if (length(to_install) > 0) {
  cat("Installing CRAN packages:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  cat("All CRAN packages already installed.\n")
}

# --- Bioconductor packages ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "SingleR",   # Automated cell type annotation
  "celldex",   # Reference datasets for SingleR
  "scran",     # Additional normalization methods
  "scater"     # QC utilities
)

to_install_bioc <- bioc_packages[!bioc_packages %in% installed]

if (length(to_install_bioc) > 0) {
  cat("Installing Bioconductor packages:", paste(to_install_bioc, collapse = ", "), "\n")
  BiocManager::install(to_install_bioc, ask = FALSE)
} else {
  cat("All Bioconductor packages already installed.\n")
}

cat("\n✓ All packages installed. You're ready to run the pipeline!\n")
cat("Next step: source('scripts/run_pipeline.R')\n")
