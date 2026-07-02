# assemble_figure2.R
# Collation script: loads phenotype + PRS data, generates all panels, and
# assembles the multi-panel PRS figure using patchwork.
#
# Usage:
#   Rscript scripts/analysis/assemble_figure2.R
#
# Or with custom paths:
#   Rscript scripts/analysis/assemble_figure2.R \
#     --pheno-path phenotypes/aou_t2d_cad_phenotypes.csv \
#     --prs-path results/prs/output/prs_all_clusters.tsv \
#     --output-dir results/figures
#
# Or source interactively:
#   source("scripts/analysis/assemble_figure2.R")
#   fig <- assemble_figure2()

library(tidyverse)
library(data.table)
library(patchwork)

# Source panel scripts (relative to project root)
script_dir <- "scripts/analysis"
source(file.path(script_dir, "figure_utils.R"))
source(file.path(script_dir, "figure2_utils.R"))
source(file.path(script_dir, "figure2_panel_a_density.R"))
source(file.path(script_dir, "figure2_panel_b_corr.R"))
source(file.path(script_dir, "figure2_panel_c_forest.R"))
source(file.path(script_dir, "figure2_panel_d_liability.R"))
source(file.path(script_dir, "figure2_panel_e_auc.R"))


#' Assemble all panels into the final multi-panel PRS figure
#'
#' @param pheno_path  Path to phenotype CSV
#' @param prs_path    Path to PRS TSV (output of compute_prs.R)
#' @param output_dir  Path for figure output
#' @param width       Figure width in inches
#' @param height      Figure height in inches
#' @return Invisible patchwork object
assemble_figure2 <- function(
    pheno_path = "phenotypes/aou_t2d_cad_phenotypes.csv",
    prs_path   = "results/prs/output/prs_all_clusters.tsv",
    output_dir = "results/figures",
    width      = 18,
    height     = 24) {

  # Load and merge data
  cat("Loading and merging phenotype + PRS data...\n")
  prs_data <- load_prs_data(pheno_path, prs_path)
  dat      <- prs_data$dat
  prs_meta <- prs_data$prs_meta

  # Generate panels
  cat("\nGenerating Panel A (PRS density curves)...\n")
  p_a <- plot_panel_a_density(dat, prs_meta)

  cat("\nGenerating Panel B (Pearson correlation)...\n")
  p_b <- plot_panel_b_corr(dat, prs_meta)

  cat("\nGenerating Panel C (Forest plots)...\n")
  p_c <- plot_panel_c_forest(dat, prs_meta)

  cat("\nGenerating Panel D (Liability R-squared)...\n")
  p_d <- plot_panel_d_liability(dat, prs_meta)

  cat("\nGenerating Panel E (AUC barplot)...\n")
  p_e <- plot_panel_e_auc(dat, prs_meta)

  # Assemble with patchwork
  # Layout:
  #   A A B B
  #   A A B B
  #   C C D D
  #   C C D D
  #   C C D D
  #   E E E E
  #   E E E E
  layout <- "
AABB
AABB
CCDD
CCDD
CCDD
EEEE
EEEE
"

  final <- p_a + p_b + p_c + p_d + p_e +
    plot_layout(design = layout) +
    plot_annotation(
      tag_levels = "A",
      theme = theme(
        plot.tag = element_text(size = 22, face = "bold")
      )
    )

  # Save main figure (PNG only)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  png_path <- file.path(output_dir, "figure2_prs.png")
  ggsave(png_path, final, width = width, height = height, dpi = 300)
  cat(sprintf("\nSaved main figure: %s\n", png_path))

  # --- Supplementary figures ---
  cat("\nGenerating supplementary figures...\n")

  # Panel B: ancestry-specific correlation matrices
  for (anc in c("eur", "afr")) {
    p <- plot_panel_b_corr(dat, prs_meta, ancestry_filter = anc)
    out <- file.path(output_dir, sprintf("figure2_corr_%s.png", anc))
    ggsave(out, p, width = 10, height = 9, dpi = 300)
    cat(sprintf("  Saved: %s\n", out))
  }
  # Also save overall correlation as standalone
  out <- file.path(output_dir, "figure2_corr_overall.png")
  ggsave(out, p_b, width = 10, height = 9, dpi = 300)
  cat(sprintf("  Saved: %s\n", out))

  # Panel C: ancestry-specific forest plots
  for (anc in c("eur", "afr")) {
    p <- plot_panel_c_forest(dat, prs_meta, ancestry_filter = anc)
    out <- file.path(output_dir, sprintf("figure2_forest_%s.png", anc))
    ggsave(out, p, width = 10, height = 14, dpi = 300)
    cat(sprintf("  Saved: %s\n", out))
  }
  out <- file.path(output_dir, "figure2_forest_overall.png")
  ggsave(out, p_c, width = 10, height = 14, dpi = 300)
  cat(sprintf("  Saved: %s\n", out))

  # Panel D: Nagelkerke R² version
  cat("\nGenerating Panel D (Nagelkerke R-squared)...\n")
  p_d_nag <- plot_panel_d_liability(dat, prs_meta, metric = "nagelkerke")
  out <- file.path(output_dir, "figure2_nagelkerke_r2.png")
  ggsave(out, p_d_nag, width = 10, height = 10, dpi = 300)
  cat(sprintf("  Saved: %s\n", out))

  invisible(final)
}


# --- Command-line execution ---
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  pheno_path <- "phenotypes/aou_t2d_cad_phenotypes.csv"
  prs_path   <- "results/prs/output/prs_all_clusters.tsv"
  output_dir <- "results/figures"

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--pheno-path" && i < length(args)) {
      pheno_path <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--prs-path" && i < length(args)) {
      prs_path <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--output-dir" && i < length(args)) {
      output_dir <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }

  assemble_figure2(pheno_path = pheno_path, prs_path = prs_path,
                   output_dir = output_dir)
}
