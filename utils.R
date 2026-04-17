# ============================================================
# utils.R — Shared helper functions
# Sourced automatically by each pipeline step.
# ============================================================

library(yaml)
library(ggplot2)

# ── Load config ─────────────────────────────────────────────
load_config <- function(config_path = "config/params.yaml") {
  if (!file.exists(config_path)) {
    stop("Config file not found: ", config_path,
         "\nMake sure you're running from the project root directory.")
  }
  config <- yaml::read_yaml(config_path)
  message("Config loaded from: ", config_path)
  return(config)
}

# ── Ensure output directories exist ─────────────────────────
setup_dirs <- function(config) {
  dirs <- c(
    config$data$output_dir,
    config$output$figures_dir,
    config$output$tables_dir,
    config$output$reports_dir
  )
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
      message("Created directory: ", d)
    }
  }
}

# ── Consistent ggplot2 theme ────────────────────────────────
theme_scrna <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(colour = "grey40", size = base_size - 1),
      axis.title    = element_text(size = base_size),
      legend.title  = element_text(face = "bold", size = base_size - 1),
      legend.text   = element_text(size = base_size - 2),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text    = element_text(face = "bold")
    )
}

# ── Save a plot in multiple formats ─────────────────────────
save_plot <- function(plot, filename, config, width = NULL, height = NULL) {
  w <- width  %||% config$output$figure_width
  h <- height %||% config$output$figure_height
  dpi <- config$output$dpi %||% 300

  for (fmt in config$output$save_formats) {
    out_path <- file.path(config$output$figures_dir,
                          paste0(filename, ".", fmt))
    ggsave(out_path, plot = plot, width = w, height = h,
           dpi = dpi, device = fmt)
    message("Saved: ", out_path)
  }
}

# ── Save a data frame as CSV ────────────────────────────────
save_table <- function(df, filename, config) {
  out_path <- file.path(config$output$tables_dir,
                        paste0(filename, ".csv"))
  write.csv(df, out_path, row.names = TRUE)
  message("Saved table: ", out_path)
}

# ── Load or save Seurat object ───────────────────────────────
save_seurat <- function(obj, step_name, config) {
  out_path <- file.path(config$data$output_dir,
                        paste0(step_name, ".rds"))
  saveRDS(obj, out_path)
  message("Seurat object saved: ", out_path)
}

load_seurat <- function(step_name, config) {
  path <- file.path(config$data$output_dir,
                    paste0(step_name, ".rds"))
  if (!file.exists(path)) {
    stop("Seurat object not found: ", path,
         "\nHave you run the previous pipeline step?")
  }
  obj <- readRDS(path)
  message("Loaded Seurat object: ", path,
          " (", ncol(obj), " cells, ", nrow(obj), " genes)")
  return(obj)
}

# ── Null-coalescing operator ─────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Print pipeline step banner ───────────────────────────────
step_banner <- function(step_num, step_name) {
  line <- strrep("─", 55)
  cat("\n", line, "\n", sep = "")
  cat(sprintf("  Step %s: %s\n", step_num, step_name))
  cat(line, "\n\n", sep = "")
}

message("utils.R loaded.")
