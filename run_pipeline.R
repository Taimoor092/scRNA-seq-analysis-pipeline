# ============================================================
# run_pipeline.R
# Master script — runs all 7 steps in order.
# Run this from the project root directory:
#   setwd("/path/to/scRNA-seq-pipeline")
#   source("scripts/run_pipeline.R")
#
# To run individual steps, source them directly:
#   source("R/02_quality_control.R")
# ============================================================

cat("
╔══════════════════════════════════════════════════╗
║        scRNA-seq Analysis Pipeline               ║
║        github.com/YOUR_USERNAME/scRNA-seq        ║
╚══════════════════════════════════════════════════╝
")

start_time <- proc.time()

# ── Verify working directory ──────────────────────────────────
if (!file.exists("config/params.yaml")) {
  stop(
    "Config file not found.\n",
    "Please run this script from the project root directory:\n",
    "  setwd('/path/to/scRNA-seq-pipeline')\n",
    "  source('scripts/run_pipeline.R')"
  )
}

# ── Load config ───────────────────────────────────────────────
source("R/utils.R")
config <- load_config()

cat("Project  :", config$project$name, "\n")
cat("Demo mode:", config$project$use_demo_data, "\n")
cat("Norm method:", config$normalization$method, "\n")
cat("Annotation:", config$annotation$method, "\n\n")

# ── Run pipeline steps ────────────────────────────────────────
steps <- list(
  list(file = "R/01_load_data.R",     name = "Load Data"),
  list(file = "R/02_quality_control.R", name = "Quality Control"),
  list(file = "R/03_normalization.R",  name = "Normalization"),
  list(file = "R/04_dim_reduction.R",  name = "Dimensionality Reduction"),
  list(file = "R/05_clustering.R",     name = "Clustering"),
  list(file = "R/06_annotation.R",     name = "Cell Type Annotation"),
  list(file = "R/07_visualization.R",  name = "Visualization")
)

for (i in seq_along(steps)) {
  step <- steps[[i]]
  cat(sprintf("\n[%d/%d] Running: %s\n", i, length(steps), step$name))
  tryCatch(
    source(step$file),
    error = function(e) {
      cat("✗ ERROR in step", i, "(", step$name, "):\n")
      cat(conditionMessage(e), "\n")
      cat("Fix the error and re-run from step", i, ":\n")
      cat("  source('", step$file, "')\n", sep = "")
      stop("Pipeline halted at step ", i)
    }
  )
}

# ── Final summary ─────────────────────────────────────────────
elapsed <- proc.time() - start_time
mins    <- floor(elapsed["elapsed"] / 60)
secs    <- round(elapsed["elapsed"] %% 60)

cat("\n╔══════════════════════════════════════════════════╗\n")
cat("║  ✓ Pipeline complete!                            ║\n")
cat(sprintf("║  Time: %d min %d sec                             ║\n", mins, secs))
cat("╚══════════════════════════════════════════════════╝\n\n")

cat("Output files:\n")
cat("  Figures   →", config$output$figures_dir, "\n")
cat("  Tables    →", config$output$tables_dir, "\n")
cat("  Objects   →", config$data$output_dir, "\n\n")

cat("Final annotated Seurat object:\n")
cat("  data/processed/06_annotated.rds\n\n")

cat("To explore results in R:\n")
cat("  obj <- readRDS('data/processed/06_annotated.rds')\n")
cat("  DimPlot(obj, group.by = 'cell_type', label = TRUE)\n")
