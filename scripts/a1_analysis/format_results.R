# format_results.R
# Summarize bNMF results, generate output files and heatmaps.

library(data.table)


#' Summarize bNMF results across replicates
#'
#' Determines optimal K, extracts best W/H matrices, writes output files.
#'
#' @param results_list List of bNMF results from run_bnmf()
#' @param trait_names Character vector of original trait names (before non-neg expansion)
#' @param variant_ids Character vector of variant IDs
#' @param output_dir Directory to write output files
#' @return List with optimal K, W matrix, H matrix
summarize_bnmf <- function(results_list, trait_names, variant_ids, output_dir, ancestry) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Determine optimal K (most frequent converged K)
  k_values <- sapply(results_list, function(r) r$K_converged)
  k_table <- table(k_values)
  optimal_k <- as.integer(names(k_table)[which.max(k_table)])
  cat(sprintf("Optimal K = %d (appeared in %d/%d replicates)\n",
              optimal_k, max(k_table), length(results_list)))

  # Select best run (among those with optimal K, pick lowest final error)
  optimal_runs <- which(k_values == optimal_k)
  final_errors <- sapply(optimal_runs, function(i) {
    errs <- results_list[[i]]$n.error
    errs[[length(errs)]]
  })
  best_idx <- optimal_runs[which.min(final_errors)]
  best_result <- results_list[[best_idx]]
  cat(sprintf("Best run: replicate %d (error = %.2f)\n", best_idx, min(final_errors)))

  # Extract W and H, keeping only active clusters (lambda >= lambda.cut)
  W_full <- best_result$W
  H_full <- best_result$H
  lambda <- best_result$lambda
  lambda_cut <- best_result$lambda.cut
  active_clusters <- which(lambda >= lambda_cut)

  # Handle edge case: no active clusters survived
  if (length(active_clusters) == 0) {
    cat("WARNING: No active clusters found (K=0). Using top cluster by lambda.\n")
    active_clusters <- which.max(lambda)
    optimal_k <- 1L
  }

  W <- W_full[, active_clusters, drop = FALSE]
  H <- H_full[active_clusters, , drop = FALSE]

  # Label columns/rows
  cluster_names <- paste0("K", seq_along(active_clusters))
  colnames(W) <- cluster_names
  rownames(W) <- variant_ids

  # Map non-negative column names back to trait labels
  nonneg_names <- c(rbind(paste0(trait_names, "_pos"), paste0(trait_names, "_neg")))
  # Filter to active nodes
  active_nodes <- best_result$active_nodes
  nonneg_names_active <- nonneg_names[active_nodes]
  rownames(H) <- cluster_names
  colnames(H) <- nonneg_names_active

  # Write W matrix
  w_dt <- data.table(VAR_ID = variant_ids, W)
  w_file <- file.path(output_dir, sprintf("W_matrix_%s.tsv", ancestry))
  fwrite(w_dt, w_file, sep = "\t")
  cat(sprintf("Wrote W matrix: %s (%d variants x %d clusters)\n",
              w_file, nrow(W), ncol(W)))

  # Write H matrix
  h_dt <- data.table(Cluster = cluster_names, H)
  h_file <- file.path(output_dir, sprintf("H_matrix_%s.tsv", ancestry))
  fwrite(h_dt, h_file, sep = "\t")
  cat(sprintf("Wrote H matrix: %s (%d clusters x %d trait columns)\n",
              h_file, nrow(H), ncol(H)))

  # Write run summary
  run_summary <- data.table(
    replicate = seq_along(results_list),
    K_converged = k_values,
    final_error = sapply(seq_along(results_list), function(i) {
      errs <- results_list[[i]]$n.error
      errs[[length(errs)]]
    }),
    final_evidence = sapply(seq_along(results_list), function(i) {
      evid <- results_list[[i]]$n.evid
      evid[[length(evid)]]
    }),
    iterations = sapply(results_list, function(r) {
      if (!is.null(r$iterations)) r$iterations else NA_integer_
    }),
    converged = sapply(results_list, function(r) {
      if (!is.null(r$converged)) r$converged else NA
    }),
    is_optimal_k = k_values == optimal_k
  )
  summary_file <- file.path(output_dir, sprintf("run_summary_%s.csv", ancestry))
  fwrite(run_summary, summary_file)
  cat(sprintf("Wrote run summary: %s\n", summary_file))

  list(
    optimal_k = optimal_k,
    W = W,
    H = H,
    best_replicate = best_idx,
    run_summary = run_summary
  )
}


#' Generate heatmaps for W and H matrices
#'
#' @param W Variant weight matrix (variants x clusters)
#' @param H Trait weight matrix (clusters x trait columns)
#' @param output_dir Directory to save PDFs
#' @param top_n_variants Number of top variants per cluster to show in W heatmap
plot_heatmaps <- function(W, H, output_dir, ancestry, top_n_variants = 50) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    cat("pheatmap not installed — skipping heatmap generation.\n")
    cat("Install with: install.packages('pheatmap')\n")
    return(invisible(NULL))
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # H matrix heatmap (trait-cluster weights)
  h_file <- file.path(output_dir, sprintf("heatmap_H_%s.png", ancestry))
  tryCatch({
    png(h_file, width = max(8, ncol(H) * 0.5), height = max(6, nrow(H) * 0.8), units = "in", res = 300)
    pheatmap::pheatmap(
      H,
      main = "Trait-Cluster Weights (H matrix)",
      cluster_rows = nrow(H) > 1,
      cluster_cols = ncol(H) > 1,
      color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
      fontsize_row = 10,
      fontsize_col = 8,
      border_color = NA
    )
    dev.off()
    cat(sprintf("Wrote heatmap: %s\n", h_file))
  }, error = function(e) {
    cat(sprintf("Error generating H heatmap: %s\n", e$message))
    if (dev.cur() > 1) dev.off()
  })

  # W matrix heatmap (top variants per cluster)
  w_file <- file.path(output_dir, sprintf("heatmap_W_%s.png", ancestry))
  tryCatch({
    # Select top variants by max weight across any cluster
    max_weights <- apply(W, 1, max)
    top_idx <- head(order(max_weights, decreasing = TRUE), top_n_variants)
    W_top <- W[top_idx, , drop = FALSE]

    png(w_file, width = max(8, ncol(W_top) * 0.8), height = max(8, nrow(W_top) * 0.15), units = "in", res = 300)
    pheatmap::pheatmap(
      W_top,
      main = sprintf("Top %d Variant-Cluster Weights (W matrix)", nrow(W_top)),
      cluster_rows = nrow(W_top) > 1,
      cluster_cols = ncol(W_top) > 1,
      color = colorRampPalette(c("white", "orange", "red3"))(100),
      fontsize_row = 6,
      fontsize_col = 10,
      border_color = NA
    )
    dev.off()
    cat(sprintf("Wrote heatmap: %s\n", w_file))
  }, error = function(e) {
    cat(sprintf("Error generating W heatmap: %s\n", e$message))
    if (dev.cur() > 1) dev.off()
  })

  invisible(NULL)
}


#' Generate text summary report
#'
#' @param summary_result Output from summarize_bnmf()
#' @param output_dir Directory to write summary
generate_summary_report <- function(summary_result, output_dir, ancestry) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  report_file <- file.path(output_dir, sprintf("summary_%s.txt", ancestry))

  lines <- c(
    "=== bNMF Analysis Summary ===",
    sprintf("Date: %s", Sys.time()),
    "",
    sprintf("Optimal K: %d", summary_result$optimal_k),
    sprintf("Best replicate: %d", summary_result$best_replicate),
    sprintf("W matrix dimensions: %d variants x %d clusters",
            nrow(summary_result$W), ncol(summary_result$W)),
    sprintf("H matrix dimensions: %d clusters x %d trait columns",
            nrow(summary_result$H), ncol(summary_result$H)),
    "",
    "K values across replicates:",
    paste("  ", capture.output(print(summary_result$run_summary)), collapse = "\n"),
    "",
    "Top traits per cluster (highest H weights):"
  )

  H <- summary_result$H
  for (k in seq_len(nrow(H))) {
    top_traits <- head(sort(H[k, ], decreasing = TRUE), 5)
    trait_str <- paste(sprintf("%s (%.3f)", names(top_traits), top_traits), collapse = ", ")
    lines <- c(lines, sprintf("  %s: %s", rownames(H)[k], trait_str))
  }

  writeLines(lines, report_file)
  cat(sprintf("Wrote summary report: %s\n", report_file))
}
