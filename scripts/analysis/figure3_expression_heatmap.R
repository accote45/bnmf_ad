# figure3_expression_heatmap.R
# Gene expression heatmap: tests whether genes within bNMF clusters show
# elevated absolute expression in Tabula Sapiens cell types vs background.
#
# Methodology (following GeneExpressionLandscape):
#   - One-sided t-test (alternative = "greater") on log2(expr + 1)
#   - Empirical permutation test (10K permutations)
#   - Heatmap of -log10(permutation p-value)
#
# Usage (CLI):
#   Rscript scripts/analysis/figure3_expression_heatmap.R \
#     --gmt-dir /path/to/expression_specificity_gmt/ \
#     --expr-dir /path/to/tabula_sapiens_cell_types/ \
#     --gtf-path /path/to/Homo_sapiens.GRCh37.75.gtf \
#     --results-dir results \
#     --ancestries EUR,AFR \
#     --top-n 10 \
#     --n-perm 10000 \
#     --output-dir results/figures
#
# Or interactive:
#   source("scripts/analysis/figure3_expression_heatmap.R")
#   generate_figure3(gmt_dir = ..., expr_dir = ..., gtf_path = ...)

library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# Source shared utilities
script_dir <- "scripts/analysis"
source(file.path(script_dir, "figure_utils.R"))
source(file.path(script_dir, "figure1_panel_g_circos.R"))


# =============================================================================
# 1. Extended GTF parser (adds Ensembl gene_id)
# =============================================================================

#' Parse GTF for protein-coding genes with both gene_name and Ensembl gene_id
#'
#' Extends parse_gtf_genes() from figure1_panel_g_circos.R to also extract the
#' Ensembl gene_id, needed to bridge SNP-mapped gene symbols to expression data.
#'
#' @param gtf_path Path to Ensembl GTF file (GRCh37), plain or .gz
#' @return tibble with columns: chr, start, end, gene_name, gene_id
parse_gtf_genes_extended <- function(gtf_path) {
  gtf <- read_tsv(
    gtf_path, comment = "#",
    col_names = c("chr", "source", "feature", "start", "end",
                  "score", "strand", "frame", "attributes"),
    col_types = cols(
      chr = col_character(), start = col_double(), end = col_double(),
      .default = col_character()
    )
  )

  gtf %>%
    filter(feature == "gene", str_detect(attributes, 'gene_biotype "protein_coding"')) %>%
    mutate(
      gene_name = str_replace(attributes, '.*gene_name "([^"]+)".*', "\\1"),
      gene_id   = str_replace(attributes, '.*gene_id "([^"]+)".*', "\\1"),
      gene_id   = str_replace(gene_id, "\\.\\d+$", "")   # strip version suffix
    ) %>%
    filter(chr %in% c(as.character(1:22), "X", "Y")) %>%
    select(chr, start, end, gene_name, gene_id)
}


# =============================================================================
# 2. Cluster gene extraction
# =============================================================================

#' Get top N genes per bNMF cluster with Ensembl IDs
#'
#' Reads W_matrix for one ancestry, maps SNPs to genes via GTF overlap,
#' ranks by mean weight, and returns top N genes per cluster.
#'
#' @param results_dir Path to results directory
#' @param ancestry Character ancestry code (e.g., "EUR")
#' @param gene_df tibble from parse_gtf_genes_extended()
#' @param top_n Number of top genes per cluster (default: 10)
#' @return tibble with columns: cluster, gene_name, gene_id, mean_weight, rank
get_top_cluster_genes <- function(results_dir, ancestry, gene_df, top_n = 10) {
  w_path <- file.path(results_dir, ancestry, "W_matrix.tsv")
  if (!file.exists(w_path)) {
    # try alternate naming convention
    w_path <- file.path(results_dir, ancestry, sprintf("W_matrix_%s.tsv", ancestry))
  }
  if (!file.exists(w_path)) {
    warning(sprintf("W_matrix not found for %s", ancestry))
    return(NULL)
  }

  w_df <- read_tsv(w_path, show_col_types = FALSE) %>%
    mutate(
      chr = str_replace(VAR_ID, "^([^_]+)_.*", "\\1"),
      pos = as.numeric(str_replace(VAR_ID, "^[^_]+_([^_]+)_.*", "\\1"))
    )
  cluster_cols <- setdiff(colnames(w_df), c("VAR_ID", "chr", "pos"))

  mapped <- map_snps_to_genes(w_df, gene_df)
  cat(sprintf("  %s: %d SNPs mapped to genes (of %d total)\n",
              ancestry, n_distinct(mapped$VAR_ID), nrow(w_df)))

  # For each cluster, get top N genes by mean weight
  map_dfr(cluster_cols, function(kcol) {
    mapped %>%
      group_by(gene_name) %>%
      summarise(
        mean_weight = mean(.data[[kcol]], na.rm = TRUE),
        gene_id = first(gene_id),
        .groups = "drop"
      ) %>%
      filter(mean_weight > 0) %>%
      slice_max(mean_weight, n = top_n) %>%
      mutate(
        cluster = kcol,
        rank = row_number()
      ) %>%
      select(cluster, gene_name, gene_id, mean_weight, rank)
  })
}


# =============================================================================
# 3. GMT loading (top-decile gene extraction)
# =============================================================================

#' Parse GMT files and extract top-decile (:10) gene sets per cell type
#'
#' @param gmt_dir Directory containing GMT files
#' @return tibble with columns: cell_type, gene_id (Ensembl IDs from :10 rows)
load_gmt_top_decile <- function(gmt_dir) {
  gmt_files <- list.files(gmt_dir, pattern = "\\.gmt$", full.names = TRUE)
  if (length(gmt_files) == 0) {
    stop(sprintf("No GMT files found in %s", gmt_dir))
  }

  cat(sprintf("  Reading %d GMT file(s) from %s\n", length(gmt_files), gmt_dir))

  map_dfr(gmt_files, function(gmt_file) {
    lines <- readLines(gmt_file)
    # Filter to :10 rows (top expression decile)
    top_lines <- lines[str_detect(lines, ":10\\t")]

    map_dfr(top_lines, function(line) {
      fields <- str_split(line, "\t")[[1]]
      cell_type <- str_replace(fields[1], ":10$", "")
      # GMT format: name, description, gene1, gene2, ...
      gene_ids <- fields[3:length(fields)]
      gene_ids <- gene_ids[gene_ids != "" & !is.na(gene_ids)]
      # Strip version suffix from Ensembl IDs
      gene_ids <- str_replace(gene_ids, "\\.\\d+$", "")
      tibble(cell_type = cell_type, gene_id = gene_ids)
    })
  })
}


# =============================================================================
# 4. Expression data loading
# =============================================================================

#' Load expression data from a single CSV file
#'
#' Reads a per-tissue expression CSV and computes median expression per gene.
#' Handles both gene x sample matrices and pre-aggregated formats.
#'
#' @param csv_path Path to expression CSV
#' @return tibble with columns: gene_id, median_expr
load_tissue_expression <- function(csv_path) {
  expr <- read_csv(csv_path, show_col_types = FALSE)

  # First column is gene identifier
  id_col <- colnames(expr)[1]
  data_cols <- setdiff(colnames(expr), id_col)

  expr <- expr %>% rename(gene_id = !!sym(id_col))
  expr$gene_id <- str_replace(expr$gene_id, "\\.\\d+$", "")

  if (length(data_cols) == 1) {
    # Already aggregated: single value column
    expr %>%
      rename(median_expr = !!sym(data_cols[1])) %>%
      select(gene_id, median_expr)
  } else {
    # Multiple samples/cells: compute median per gene
    expr %>%
      rowwise() %>%
      mutate(median_expr = median(c_across(all_of(data_cols)), na.rm = TRUE)) %>%
      ungroup() %>%
      select(gene_id, median_expr)
  }
}


#' Load all expression CSVs from a directory
#'
#' @param expr_dir Directory containing *_expression.csv files
#' @return tibble with columns: cell_type, gene_id, median_expr
load_all_expression <- function(expr_dir) {
  csv_files <- list.files(expr_dir, pattern = "_expression\\.csv$",
                          full.names = TRUE, recursive = TRUE)
  if (length(csv_files) == 0) {
    # Fall back to any CSV files
    csv_files <- list.files(expr_dir, pattern = "\\.csv$",
                            full.names = TRUE, recursive = TRUE)
  }
  if (length(csv_files) == 0) {
    stop(sprintf("No expression CSV files found in %s", expr_dir))
  }

  cat(sprintf("  Loading %d expression file(s) from %s\n",
              length(csv_files), expr_dir))

  map_dfr(csv_files, function(csv_path) {
    # Derive cell type name from filename
    cell_type <- basename(csv_path) %>%
      str_replace("_expression\\.csv$", "") %>%
      str_replace("\\.csv$", "")

    tryCatch({
      te <- load_tissue_expression(csv_path)
      te %>% mutate(cell_type = cell_type)
    }, error = function(e) {
      warning(sprintf("Failed to load %s: %s", basename(csv_path), e$message))
      tibble(gene_id = character(), median_expr = numeric(), cell_type = character())
    })
  })
}


#' Load expression data from a single MAGMA-format matrix file
#'
#' Fallback for when per-tissue CSVs are unavailable (e.g., local testing).
#'
#' @param magma_path Path to MAGMA-format file (TSV: Name + tissue columns)
#' @return tibble with columns: cell_type, gene_id, median_expr
load_expression_magma <- function(magma_path) {
  cat(sprintf("  Loading MAGMA-format expression from %s\n", magma_path))
  expr <- read_tsv(magma_path, show_col_types = FALSE)

  id_col <- colnames(expr)[1]
  tissue_cols <- setdiff(colnames(expr), id_col)

  expr %>%
    rename(gene_id = !!sym(id_col)) %>%
    mutate(gene_id = str_replace(gene_id, "\\.\\d+$", "")) %>%
    pivot_longer(
      cols = all_of(tissue_cols),
      names_to = "cell_type",
      values_to = "median_expr"
    )
}


# =============================================================================
# 5. Statistical testing
# =============================================================================

#' Test if cluster genes have higher expression than background in one cell type
#'
#' Performs a one-sided t-test and empirical permutation test.
#'
#' @param cluster_gene_ids Character vector of Ensembl IDs for cluster genes
#' @param expr_vec Named numeric vector: gene_id -> expression value (one cell type)
#' @param n_perm Number of permutations (default: 10000)
#' @param seed Random seed
#' @return list with t_pvalue, perm_pvalue, n_cluster_found, n_background, mean_cluster, mean_bg
run_expression_test <- function(cluster_gene_ids, expr_vec, n_perm = 10000, seed = 42) {
  # Partition into cluster and background
  found_ids <- intersect(cluster_gene_ids, names(expr_vec))
  bg_ids <- setdiff(names(expr_vec), cluster_gene_ids)

  n_found <- length(found_ids)
  n_bg <- length(bg_ids)

  if (n_found < 2) {
    return(list(t_pvalue = NA_real_, perm_pvalue = NA_real_,
                n_cluster_found = n_found, n_background = n_bg,
                mean_cluster = NA_real_, mean_bg = NA_real_))
  }

  cluster_expr <- log2(expr_vec[found_ids] + 1)
  bg_expr <- log2(expr_vec[bg_ids] + 1)

  # One-sided t-test
  t_result <- tryCatch(
    t.test(cluster_expr, bg_expr, alternative = "greater"),
    error = function(e) list(p.value = NA_real_)
  )

  # Permutation test
  all_expr <- c(cluster_expr, bg_expr)
  n_total <- length(all_expr)
  obs_diff <- mean(cluster_expr) - mean(bg_expr)

  set.seed(seed)
  perm_diffs <- replicate(n_perm, {
    idx <- sample.int(n_total, n_found)
    mean(all_expr[idx]) - mean(all_expr[-idx])
  })
  perm_pvalue <- (sum(perm_diffs >= obs_diff) + 1) / (n_perm + 1)

  list(
    t_pvalue        = t_result$p.value,
    perm_pvalue     = perm_pvalue,
    n_cluster_found = n_found,
    n_background    = n_bg,
    mean_cluster    = mean(cluster_expr),
    mean_bg         = mean(bg_expr)
  )
}


#' Run expression tests across all cluster x cell type combinations
#'
#' @param cluster_genes_df tibble with cluster, gene_id columns
#' @param expr_long tibble with cell_type, gene_id, median_expr columns
#' @param n_perm Number of permutations
#' @param seed Random seed
#' @return tibble with cluster, cell_type, t_pvalue, perm_pvalue, etc.
run_all_expression_tests <- function(cluster_genes_df, expr_long,
                                     n_perm = 10000, seed = 42) {
  clusters <- unique(cluster_genes_df$cluster)
  cell_types <- unique(expr_long$cell_type)
  n_tests <- length(clusters) * length(cell_types)

  cat(sprintf("  Running %d tests (%d clusters x %d cell types, %d permutations each)\n",
              n_tests, length(clusters), length(cell_types), n_perm))

  counter <- 0
  results <- list()

  for (cl in clusters) {
    cl_gene_ids <- cluster_genes_df %>%
      filter(cluster == cl) %>%
      pull(gene_id)

    for (ct in cell_types) {
      counter <- counter + 1
      if (counter %% 20 == 0 || counter == n_tests) {
        cat(sprintf("    %d / %d\r", counter, n_tests))
      }

      # Build named expression vector for this cell type
      ct_data <- expr_long %>% filter(cell_type == ct)
      expr_vec <- setNames(ct_data$median_expr, ct_data$gene_id)

      res <- run_expression_test(cl_gene_ids, expr_vec, n_perm, seed)

      results[[counter]] <- tibble(
        cluster         = cl,
        cell_type       = ct,
        t_pvalue        = res$t_pvalue,
        perm_pvalue     = res$perm_pvalue,
        n_cluster_found = res$n_cluster_found,
        n_background    = res$n_background,
        mean_cluster    = res$mean_cluster,
        mean_bg         = res$mean_bg
      )
    }
  }
  cat("\n")

  bind_rows(results)
}


# =============================================================================
# 6. Heatmap visualization
# =============================================================================

#' Plot expression association heatmap using pheatmap
#'
#' @param test_results tibble from run_all_expression_tests()
#' @param ancestry Character ancestry code for title
#' @param output_dir Directory for output files
#' @param pvalue_col Which p-value column: "perm_pvalue" or "t_pvalue"
#' @param sig_threshold P-value threshold for significance markers (default: 0.05)
#' @return Invisible pheatmap object
plot_expression_heatmap <- function(test_results, ancestry, output_dir,
                                    pvalue_col = "perm_pvalue",
                                    sig_threshold = 0.05) {
  # Pivot to matrix: rows = cell_type, columns = cluster
  mat_df <- test_results %>%
    mutate(neg_log10_p = -log10(pmax(.data[[pvalue_col]], 1e-10))) %>%
    select(cell_type, cluster, neg_log10_p) %>%
    pivot_wider(names_from = cluster, values_from = neg_log10_p)

  cell_types <- mat_df$cell_type
  mat <- mat_df %>% select(-cell_type) %>% as.matrix()
  rownames(mat) <- cell_types

  mat[is.na(mat)] <- 0

  # Clean up cell type names for display
  rownames(mat) <- str_replace_all(rownames(mat), "_", " ")

  # Build significance label matrix (* for p < threshold)
  sig_df <- test_results %>%
    mutate(label = ifelse(.data[[pvalue_col]] < sig_threshold, "*", "")) %>%
    select(cell_type, cluster, label) %>%
    pivot_wider(names_from = cluster, values_from = label)
  sig_labels <- sig_df %>% select(-cell_type) %>% as.matrix()
  rownames(sig_labels) <- str_replace_all(sig_df$cell_type, "_", " ")
  sig_labels[is.na(sig_labels)] <- ""

  # Color palette: white -> orange -> dark red
  n_colors <- 100
  color_palette <- colorRampPalette(c("white", "#FDAE61", "#D73027"))(n_colors)

  # Determine figure dimensions
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  fig_height <- max(6, n_rows * 0.22 + 3)
  fig_width <- max(5, n_cols * 1.5 + 4)
  fontsize_row <- max(5, min(10, 200 / n_rows))

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Save PDF
  pdf_path <- file.path(output_dir, sprintf("figure3_expression_heatmap_%s.pdf", ancestry))
  ht <- pheatmap(
    mat,
    color = color_palette,
    cluster_rows = n_rows > 2,
    cluster_cols = FALSE,
    display_numbers = sig_labels,
    fontsize_number = 12,
    fontsize_row = fontsize_row,
    fontsize_col = 12,
    main = sprintf("%s \u2014 Absolute Expression\n-log10(permutation p-value)", ancestry),
    angle_col = 0,
    filename = pdf_path,
    width = fig_width,
    height = fig_height
  )
  cat(sprintf("  Saved %s\n", pdf_path))

  # Save PNG
  png_path <- file.path(output_dir, sprintf("figure3_expression_heatmap_%s.png", ancestry))
  pheatmap(
    mat,
    color = color_palette,
    cluster_rows = n_rows > 2,
    cluster_cols = FALSE,
    display_numbers = sig_labels,
    fontsize_number = 12,
    fontsize_row = fontsize_row,
    fontsize_col = 12,
    main = sprintf("%s \u2014 Absolute Expression\n-log10(permutation p-value)", ancestry),
    angle_col = 0,
    filename = png_path,
    width = fig_width,
    height = fig_height
  )
  cat(sprintf("  Saved %s\n", png_path))

  invisible(ht)
}


# =============================================================================
# 7. Main orchestrator
# =============================================================================

#' Generate Figure 3: Expression association heatmaps for bNMF clusters
#'
#' @param gmt_dir Directory with Tabula Sapiens GMT files (expression_specificity_gmt/)
#' @param expr_dir Directory with per-tissue *_expression.csv files
#' @param expr_magma_path Alternative: single MAGMA-format file (for local testing)
#' @param gtf_path Path to GRCh37 GTF file
#' @param results_dir Path to bNMF results directory
#' @param ancestries Character vector of ancestry codes
#' @param top_n Number of top genes per cluster (default: 10)
#' @param n_perm Number of permutations (default: 10000)
#' @param output_dir Output directory
#' @param seed Random seed
#' @return Invisible list of Heatmap objects
generate_figure3 <- function(
    gmt_dir         = NULL,
    expr_dir        = NULL,
    expr_magma_path = NULL,
    gtf_path        = "/sc/arion/projects/paul_oreilly/lab/kestir01/Homo_sapiens.GRCh37.75.gtf",
    results_dir     = "results",
    ancestries      = c("EUR", "AFR"),
    top_n           = 10,
    n_perm          = 10000,
    output_dir      = "results/figures",
    seed            = 42) {

  # --- Validate inputs ---
  if (is.null(expr_dir) && is.null(expr_magma_path)) {
    stop("Must provide either --expr-dir (per-tissue CSVs) or --expr-magma (single matrix)")
  }

  # --- Step 1: Parse GTF ---
  cat("Step 1: Parsing GTF for protein-coding genes...\n")
  gene_df <- parse_gtf_genes_extended(gtf_path)
  cat(sprintf("  Found %d protein-coding genes\n", nrow(gene_df)))

  # Ensembl ID <-> gene_name mapping
  id_map <- gene_df %>% select(gene_name, gene_id) %>% distinct()

  # --- Step 2: Load expression data ---
  cat("Step 2: Loading expression data...\n")
  if (!is.null(expr_dir) && dir.exists(expr_dir)) {
    expr_long <- load_all_expression(expr_dir)
  } else if (!is.null(expr_magma_path) && file.exists(expr_magma_path)) {
    expr_long <- load_expression_magma(expr_magma_path)
  } else {
    stop("Expression data not found at provided path(s)")
  }

  n_cell_types <- n_distinct(expr_long$cell_type)
  n_genes_expr <- n_distinct(expr_long$gene_id)
  cat(sprintf("  Loaded expression for %d genes across %d cell types\n",
              n_genes_expr, n_cell_types))

  # --- Step 3: Optionally load GMT top-decile genes for context ---
  gmt_top10 <- NULL
  if (!is.null(gmt_dir) && dir.exists(gmt_dir)) {
    cat("Step 3: Loading GMT top-decile genes...\n")
    gmt_top10 <- load_gmt_top_decile(gmt_dir)
    cat(sprintf("  Found top-decile gene sets for %d cell types\n",
                n_distinct(gmt_top10$cell_type)))
  } else {
    cat("Step 3: No GMT directory provided, skipping top-decile loading\n")
  }

  # --- Step 4: Process each ancestry ---
  heatmaps <- list()

  for (anc in ancestries) {
    cat(sprintf("\n=== Processing %s ===\n", anc))

    # Get cluster genes
    cat("  Mapping SNPs to genes...\n")
    cluster_genes <- get_top_cluster_genes(results_dir, anc, gene_df, top_n)
    if (is.null(cluster_genes) || nrow(cluster_genes) == 0) {
      warning(sprintf("No cluster genes found for %s, skipping", anc))
      next
    }

    # Report coverage
    n_in_expr <- sum(cluster_genes$gene_id %in% unique(expr_long$gene_id))
    cat(sprintf("  %d / %d cluster genes found in expression data\n",
                n_in_expr, nrow(cluster_genes)))

    # Save cluster genes
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    genes_path <- file.path(output_dir, sprintf("figure3_cluster_genes_%s.tsv", anc))
    write_tsv(cluster_genes, genes_path)
    cat(sprintf("  Saved cluster genes to %s\n", genes_path))

    # Run statistical tests
    cat("  Running statistical tests...\n")
    test_results <- run_all_expression_tests(cluster_genes, expr_long, n_perm, seed)

    # Save test results
    results_path <- file.path(output_dir, sprintf("figure3_test_results_%s.tsv", anc))
    write_tsv(test_results, results_path)
    cat(sprintf("  Saved test results to %s\n", results_path))

    # Report summary
    n_sig <- sum(test_results$perm_pvalue < 0.05, na.rm = TRUE)
    cat(sprintf("  %d / %d tests significant at p < 0.05\n",
                n_sig, nrow(test_results)))

    # Generate heatmap
    cat("  Generating heatmap...\n")
    ht <- plot_expression_heatmap(test_results, anc, output_dir)
    heatmaps[[anc]] <- ht
  }

  cat("\nDone.\n")
  invisible(heatmaps)
}


# =============================================================================
# 8. CLI argument parsing
# =============================================================================

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  # Defaults
  gmt_dir         <- NULL
  expr_dir        <- NULL
  expr_magma_path <- NULL
  gtf_path        <- "/sc/arion/projects/paul_oreilly/lab/kestir01/Homo_sapiens.GRCh37.75.gtf"
  results_dir     <- "results"
  ancestries      <- c("EUR", "AFR")
  top_n           <- 10
  n_perm          <- 10000
  output_dir      <- "results/figures"
  seed            <- 42

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--gmt-dir" && i < length(args)) {
      gmt_dir <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--expr-dir" && i < length(args)) {
      expr_dir <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--expr-magma" && i < length(args)) {
      expr_magma_path <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--gtf-path" && i < length(args)) {
      gtf_path <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--results-dir" && i < length(args)) {
      results_dir <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--ancestries" && i < length(args)) {
      ancestries <- str_split(args[i + 1], ",")[[1]]; i <- i + 2
    } else if (args[i] == "--top-n" && i < length(args)) {
      top_n <- as.integer(args[i + 1]); i <- i + 2
    } else if (args[i] == "--n-perm" && i < length(args)) {
      n_perm <- as.integer(args[i + 1]); i <- i + 2
    } else if (args[i] == "--output-dir" && i < length(args)) {
      output_dir <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--seed" && i < length(args)) {
      seed <- as.integer(args[i + 1]); i <- i + 2
    } else {
      warning(sprintf("Unknown argument: %s", args[i]))
      i <- i + 1
    }
  }

  generate_figure3(
    gmt_dir         = gmt_dir,
    expr_dir        = expr_dir,
    expr_magma_path = expr_magma_path,
    gtf_path        = gtf_path,
    results_dir     = results_dir,
    ancestries      = ancestries,
    top_n           = top_n,
    n_perm          = n_perm,
    output_dir      = output_dir,
    seed            = seed
  )
}
