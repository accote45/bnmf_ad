# assemble_figure.R
# Collation script: loads data, generates all panels, and assembles the
# multi-panel figure using patchwork.
#
# Usage:
#   Rscript scripts/analysis/assemble_figure.R [--results-dir results] [--ancestries EUR,AFR,META]
#
# Or source interactively:
#   source("scripts/analysis/assemble_figure.R")
#   fig <- assemble_figure()

library(tidyverse)
library(data.table)
library(patchwork)

# Source panel scripts (now in a1_analysis/)
script_dir <- "scripts/a1_analysis"
source(file.path(script_dir, "figure_utils.R"))
source(file.path("scripts/analysis", "figure1_panel_a_heatmap.R"))
source(file.path("scripts/analysis", "figure1_panel_c_barplots.R"))
source(file.path("scripts/analysis", "figure1_panel_f_variantselection.R"))
source(file.path("scripts/analysis", "figure1_panel_g_circos.R"))


#' Check whether an H matrix has any active (non-zero) clusters
#'
#' @param h data.table H-matrix with Cluster column + trait weight columns
#' @return TRUE if at least one cluster has non-zero weights
has_active_clusters <- function(h) {
  cluster_cols <- setdiff(colnames(h), "Cluster")
  if (length(cluster_cols) == 0) return(FALSE)
  any(as.matrix(h[, ..cluster_cols]) != 0)
}


#' Assemble all panels into the final multi-panel figure
#'
#' @param results_dir Path to results directory (default: "results")
#' @param ancestries Character vector of ancestry codes (default: c("EUR", "AFR", "META"))
#' @param gtf_path Path to GRCh37 GTF file for SNP-to-gene mapping
#' @param output_dir Path for figure output (default: "results/figures")
#' @param width Figure width in inches
#' @param height Figure height in inches
#' @return Invisible patchwork object
assemble_figure <- function(results_dir = "results",
                            ancestries = c("EUR", "AFR", "META"),
                            gtf_path = "/sc/arion/projects/paul_oreilly/lab/kestir01/Homo_sapiens.GRCh37.75.gtf",
                            output_dir = "results/figures",
                            width = 20,
                            height = 20) {

  cat("Loading H matrices...\n")
  h_list <- load_h_matrices(results_dir, ancestries)
  if (length(h_list) == 0) stop("No H matrices found. Check results_dir and ancestries.")

  cat(sprintf("Loaded %d ancestries: %s\n", length(h_list), paste(names(h_list), collapse = ", ")))

  # Filter out ancestries with no active clusters (K=0)
  h_list <- Filter(has_active_clusters, h_list)
  if (length(h_list) == 0) stop("No ancestries have active clusters (all K=0).")

  active_ancestries <- names(h_list)
  cat(sprintf("Active ancestries (K>0): %s\n", paste(active_ancestries, collapse = ", ")))

  # Build panels dynamically
  panels <- list()

  cat("Generating panel: Heatmaps...\n")
  panels[["heatmap"]] <- plot_panel_a(h_list)

  cat("Generating panel: Weight bar plots...\n")
  panels[["barplots"]] <- plot_panel_c(h_list)

  cat("Generating panel: Variant selection...\n")
  panels[["variant_sel"]] <- plot_panel_f(results_dir, active_ancestries)

  cat("Generating panel: Gene bar plots...\n")
  panels[["gene_bars"]] <- plot_panel_g(results_dir, active_ancestries, gtf_path)

  # Assemble with patchwork — dynamic 2-column grid
  n_panels <- length(panels)
  cat(sprintf("Assembling %d panels in 2-column layout...\n", n_panels))

  final <- wrap_plots(panels, ncol = 2) +
    plot_annotation(
      tag_levels = "A",
      theme = theme(
        plot.tag = element_text(size = 22, face = "bold")
      )
    )

  # Save
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  pdf_path <- file.path(output_dir, "figure_main.pdf")
  ggsave(pdf_path, final, width = width, height = height, device = cairo_pdf)
  cat(sprintf("Saved PDF: %s\n", pdf_path))

  png_path <- file.path(output_dir, "figure_main.png")
  ggsave(png_path, final, width = width, height = height, dpi = 300)
  cat(sprintf("Saved PNG: %s\n", png_path))

  invisible(final)
}


# --- Command-line execution ---
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  results_dir <- "results"
  ancestries  <- c("EUR", "AFR", "META")

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--results-dir" && i < length(args)) {
      results_dir <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--ancestries" && i < length(args)) {
      ancestries <- strsplit(args[i + 1], ",")[[1]]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }

  assemble_figure(results_dir = results_dir, ancestries = ancestries)
}
