#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# 04_compare_arms.R
# ------------------------------------------------------------------
# Compares bNMF results from two arms:
#   Arm 1 (published): results/a1_analysis/EUR/
#   Arm 2 (split-UKB): results/a1_analysis_splitukb/EUR/
#
# Outputs comparison metrics, cluster matching, and figures to
# results/a1_comparison/.
#
# Usage:
#   Rscript scripts/a1_analysis/04_compare_arms.R
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

# ---- CLI arguments ----
option_list <- list(
  make_option("--arm1-dir", type = "character",
              default = "results/a1_analysis/EUR",
              help = "Directory with Arm 1 (published) results"),
  make_option("--arm2-dir", type = "character",
              default = "results/a1_analysis_splitukb/EUR",
              help = "Directory with Arm 2 (split-UKB) results"),
  make_option("--output-dir", type = "character",
              default = "results/a1_comparison",
              help = "Output directory for comparison results"),
  make_option("--arm1-label", type = "character", default = "Published",
              help = "Label for Arm 1"),
  make_option("--arm2-label", type = "character", default = "SplitUKB",
              help = "Label for Arm 2"),
  make_option("--trait-map", type = "character", default = NULL,
              help = "CSV with arm1_name,arm2_name columns to harmonize trait names across arms")
)
opt <- parse_args(OptionParser(option_list = option_list))

arm1_dir    <- opt$`arm1-dir`
arm2_dir    <- opt$`arm2-dir`
output_dir  <- opt$`output-dir`
arm1_label  <- opt$`arm1-label`
arm2_label  <- opt$`arm2-label`
trait_map_file <- opt$`trait-map`

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ---- Helper: collapse non-negative H to signed loadings ----
# H matrix has columns like Trait_pos, Trait_neg. Signed loading = pos - neg.
# strip_source: if TRUE, remove _AuthorYear suffixes and deduplicate by keeping
#   the trait version with highest mean absolute loading across clusters.
collapse_h_to_signed <- function(h_dt, strip_source = FALSE) {
  cluster_col <- names(h_dt)[1]  # "Cluster"
  clusters <- h_dt[[cluster_col]]
  h_mat <- as.matrix(h_dt[, -1, with = FALSE])

  # Extract trait names from column pairs
  all_cols <- colnames(h_mat)
  pos_cols <- grep("_pos$", all_cols, value = TRUE)
  trait_names <- sub("_pos$", "", pos_cols)

  signed <- matrix(NA, nrow = length(clusters), ncol = length(trait_names),
                   dimnames = list(clusters, trait_names))
  for (tr in trait_names) {
    pc <- paste0(tr, "_pos")
    nc <- paste0(tr, "_neg")
    if (pc %in% all_cols && nc %in% all_cols) {
      signed[, tr] <- h_mat[, pc] - h_mat[, nc]
    }
  }

  if (strip_source) {
    clean_names <- sub("_[A-Za-z]+[0-9]{4}$", "", colnames(signed))
    if (any(duplicated(clean_names))) {
      # For each duplicate base trait, keep the version with highest mean |loading|
      mean_abs <- colMeans(abs(signed))
      keep <- tapply(seq_along(clean_names), clean_names, function(idx) idx[which.max(mean_abs[idx])])
      signed <- signed[, unlist(keep), drop = FALSE]
      colnames(signed) <- clean_names[unlist(keep)]
    } else {
      colnames(signed) <- clean_names
    }
  }

  signed
}

# ---- Helper: cosine similarity ----
cosine_sim <- function(a, b) sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)) + 1e-300)

# ---- Semantic cluster labels ----
dk_labels <- c(
  K1 = "Glucose-",                   K2 = "Lpa",
  K3 = "Metabolic Syndrome",         K4 = "Body Composition-",
  K5 = "CRP ApoB",                   K6 = "Liver",
  K7 = "Erythrocyte",                K8 = "Glucose+",
  K9 = "Triglycerides-"
)
tails_labels <- c(
  K1 = "Glucose+",                   K2 = "Erythrocyte",
  K3 = "Metabolic Syndrome",         K4 = "Body Composition-",
  K5 = "Triglycerides ApoA",         K6 = "Body Composition+",
  K7 = "Reticulocyte",               K8 = "Triglycerides-"
)
core_published_labels <- c(
  K1 = "Height+",                    K2 = "Height-",
  K3 = "Body Composition-",          K4 = "Erythrocyte",
  K5 = "Glucose-",                   K6 = "Metabolic Syndrome",
  K7 = "Glucose+",                   K8 = "Triglycerides-"
)
core_splitukb_labels <- c(
  K1 = "Glucose+",                   K2 = "Body Composition-",
  K3 = "Reticulocytes",              K4 = "Metabolic Syndrome",
  K5 = "Triglycerides ApoA",         K6 = "Triglycerides-"
)

# Map arm labels to their semantic cluster names
label_registry <- list(
  Published      = dk_labels,
  SplitUKB       = tails_labels,
  Core_Published = core_published_labels,
  Core_SplitUKB  = core_splitukb_labels
)


# ---- Load H matrices ----
cat("Loading H matrices...\n")
h1_file <- file.path(arm1_dir, "H_matrix_EUR.tsv")
h2_file <- file.path(arm2_dir, "H_matrix_EUR.tsv")
stopifnot(file.exists(h1_file), file.exists(h2_file))

h1_dt <- fread(h1_file)
h2_dt <- fread(h2_file)

# Strip GWAS source suffixes (_AuthorYear) from both arms so trait names align.
# The regex is a no-op on names that lack _AuthorYear patterns (e.g., split-UKB names).
h1_signed <- collapse_h_to_signed(h1_dt, strip_source = TRUE)
h2_signed <- collapse_h_to_signed(h2_dt, strip_source = TRUE)

# Apply trait name mapping if provided (harmonize names across arms)
if (!is.null(trait_map_file)) {
  trait_map <- fread(trait_map_file)
  h1_cols <- colnames(h1_signed)
  mapped <- match(h1_cols, trait_map$arm1_name)
  h1_cols[!is.na(mapped)] <- trait_map$arm2_name[mapped[!is.na(mapped)]]
  colnames(h1_signed) <- h1_cols
  cat(sprintf("Applied trait name mapping: %d traits renamed\n", sum(!is.na(mapped))))
}

k1 <- nrow(h1_signed)
k2 <- nrow(h2_signed)
cat(sprintf("Arm 1 (%s): K = %d, traits = %d\n", arm1_label, k1, ncol(h1_signed)))
cat(sprintf("Arm 2 (%s): K = %d, traits = %d\n", arm2_label, k2, ncol(h2_signed)))


# ---- Load W matrices ----
cat("Loading W matrices...\n")
w1_file <- file.path(arm1_dir, "W_matrix_EUR.tsv")
w2_file <- file.path(arm2_dir, "W_matrix_EUR.tsv")
stopifnot(file.exists(w1_file), file.exists(w2_file))

w1_dt <- fread(w1_file)
w2_dt <- fread(w2_file)


# ---- 1. Cluster count comparison ----
cat("\n=== Cluster Count ===\n")
cat(sprintf("  %s: K = %d\n", arm1_label, k1))
cat(sprintf("  %s: K = %d\n", arm2_label, k2))


# ---- 2. Shared traits analysis ----
traits1 <- colnames(h1_signed)
traits2 <- colnames(h2_signed)
shared_traits <- intersect(traits1, traits2)
cat(sprintf("\n=== Shared Traits ===\n"))
cat(sprintf("  %s traits: %d\n", arm1_label, length(traits1)))
cat(sprintf("  %s traits: %d\n", arm2_label, length(traits2)))
cat(sprintf("  Shared: %d\n", length(shared_traits)))
if (length(shared_traits) > 0) {
  cat(sprintf("  Shared trait names: %s\n", paste(shared_traits, collapse = ", ")))
}


# ---- 3. Cluster matching (on shared traits) ----
cat("\n=== Cluster Matching ===\n")

if (length(shared_traits) >= 2) {
  h1_shared <- h1_signed[, shared_traits, drop = FALSE]
  h2_shared <- h2_signed[, shared_traits, drop = FALSE]

  # Correlation matrix: rows = Arm1 clusters, cols = Arm2 clusters
  cor_mat <- cor(t(h1_shared), t(h2_shared))
  rownames(cor_mat) <- paste0(arm1_label, "_", rownames(h1_signed))
  colnames(cor_mat) <- paste0(arm2_label, "_", rownames(h2_signed))

  cat("Correlation matrix (Arm1 clusters x Arm2 clusters):\n")
  print(round(cor_mat, 3))

  # Greedy matching: for each Arm1 cluster, find best Arm2 match
  matches <- data.table(
    Arm1_Cluster = rownames(cor_mat),
    Best_Arm2_Cluster = colnames(cor_mat)[apply(abs(cor_mat), 1, which.max)],
    Correlation = apply(cor_mat, 1, function(x) x[which.max(abs(x))])
  )
  cat("\nBest matches:\n")
  print(matches)

  # Write cluster matching
  fwrite(matches, file.path(output_dir, "cluster_matching.tsv"), sep = "\t")

  # ---- Correlation heatmap ----
  cor_long <- as.data.table(expand.grid(
    Arm1 = rownames(cor_mat), Arm2 = colnames(cor_mat),
    stringsAsFactors = FALSE
  ))
  cor_long[, Correlation := as.vector(cor_mat)]

  p_heatmap <- ggplot(cor_long, aes(x = Arm2, y = Arm1, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = "Cluster Correlation (Shared Traits)",
         x = arm2_label, y = arm1_label) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(output_dir, "H_correlation_heatmap.png"), p_heatmap,
         width = max(6, k2 * 1.2), height = max(5, k1 * 1.0), dpi = 150)
  cat("Wrote H_correlation_heatmap.png\n")

  # ---- Styled heatmaps: Pearson correlation & cosine similarity ----
  # Helper to build a heatmap from a matrix with semantic labels
  make_styled_heatmap <- function(mat, arm1_labs, arm2_labs, fill_label, title) {
    # Apply semantic labels to rows (Arm1) and cols (Arm2)
    row_ids <- sub(paste0("^", arm1_label, "_"), "", rownames(mat))
    col_ids <- sub(paste0("^", arm2_label, "_"), "", colnames(mat))
    row_labels <- ifelse(row_ids %in% names(arm1_labs), arm1_labs[row_ids], row_ids)
    col_labels <- ifelse(col_ids %in% names(arm2_labs), arm2_labs[col_ids], col_ids)

    long <- as.data.table(expand.grid(
      Row = rownames(mat), Col = colnames(mat), stringsAsFactors = FALSE
    ))
    long[, Value := as.vector(mat)]
    long[, Row_label := factor(row_labels[match(Row, rownames(mat))],
                               levels = rev(row_labels))]
    long[, Col_label := factor(col_labels[match(Col, colnames(mat))],
                               levels = col_labels)]

    ggplot(long, aes(x = Col_label, y = Row_label, fill = Value)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.2f", Value)), size = 3.2) +
      scale_fill_gradient2(low = "steelblue4", mid = "white", high = "firebrick3",
                           midpoint = 0, limits = c(-1, 1), name = fill_label) +
      labs(title = title, x = arm2_label, y = arm1_label) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white"))
  }

  # Resolve semantic labels for each arm (fall back to raw K labels)
  arm1_labs <- label_registry[[arm1_label]]
  arm2_labs <- label_registry[[arm2_label]]
  if (is.null(arm1_labs)) arm1_labs <- setNames(paste0("K", seq_len(k1)), paste0("K", seq_len(k1)))
  if (is.null(arm2_labs)) arm2_labs <- setNames(paste0("K", seq_len(k2)), paste0("K", seq_len(k2)))

  # Pearson correlation heatmap (already computed as cor_mat)
  p_pearson <- make_styled_heatmap(cor_mat, arm1_labs, arm2_labs,
                                   "Pearson r",
                                   sprintf("Pearson Correlation: %s vs %s Clusters", arm1_label, arm2_label))
  ggsave(file.path(output_dir, "heatmap_pearson_correlation.png"), p_pearson,
         width = max(8, k2 * 1.4), height = max(7, k1 * 1.2), dpi = 300)
  cat("Wrote heatmap_pearson_correlation.png\n")

  # H-matrix cosine similarity
  h_cosine_mat <- matrix(0, nrow = nrow(h1_shared), ncol = nrow(h2_shared),
                         dimnames = list(rownames(cor_mat), colnames(cor_mat)))
  for (i in seq_len(nrow(h1_shared))) {
    for (j in seq_len(nrow(h2_shared))) {
      h_cosine_mat[i, j] <- cosine_sim(h1_shared[i, ], h2_shared[j, ])
    }
  }

  p_cosine <- make_styled_heatmap(h_cosine_mat, arm1_labs, arm2_labs,
                                  "Cosine Sim",
                                  sprintf("Cosine Similarity: %s vs %s Clusters", arm1_label, arm2_label))
  ggsave(file.path(output_dir, "heatmap_cosine_similarity.png"), p_cosine,
         width = max(8, k2 * 1.4), height = max(7, k1 * 1.2), dpi = 300)
  cat("Wrote heatmap_cosine_similarity.png\n")

} else {
  cat("Fewer than 2 shared traits — skipping cluster matching on H.\n")
  matches <- NULL
}


# ---- 4. W matrix comparison (shared variants) ----
cat("\n=== Variant Overlap (W matrix) ===\n")
vars1 <- w1_dt$VAR_ID
vars2 <- w2_dt$VAR_ID
shared_vars <- intersect(vars1, vars2)
cat(sprintf("  %s variants: %d\n", arm1_label, length(vars1)))
cat(sprintf("  %s variants: %d\n", arm2_label, length(vars2)))
cat(sprintf("  Shared: %d (%.1f%% of %s, %.1f%% of %s)\n",
            length(shared_vars),
            100 * length(shared_vars) / length(vars1), arm1_label,
            100 * length(shared_vars) / length(vars2), arm2_label))

# W-matrix cosine similarity between all cluster pairs
w1_kcols <- grep("^K", colnames(w1_dt), value = TRUE)
w2_kcols <- grep("^K", colnames(w2_dt), value = TRUE)
w1_mat <- as.matrix(w1_dt[, ..w1_kcols])
w2_mat <- as.matrix(w2_dt[, ..w2_kcols])

cosine_matrix <- matrix(0, nrow = length(w1_kcols), ncol = length(w2_kcols),
                         dimnames = list(w1_kcols, w2_kcols))
for (i in seq_along(w1_kcols)) {
  for (j in seq_along(w2_kcols)) {
    cosine_matrix[i, j] <- cosine_sim(w1_mat[, i], w2_mat[, j])
  }
}

cat("\nW-matrix cosine similarity (Arm1 rows x Arm2 cols):\n")
print(round(cosine_matrix, 3))

# Write cosine similarity matrix
cos_dt <- as.data.table(cosine_matrix, keep.rownames = paste0(arm1_label, "_Cluster"))
fwrite(cos_dt, file.path(output_dir, "cosine_similarity_matrix.csv"))
cat("Wrote cosine_similarity_matrix.csv\n")

if (length(shared_vars) >= 10 && !is.null(matches)) {
  # For matched cluster pairs, correlate W vectors on shared variants
  w1_shared <- as.matrix(w1_dt[match(shared_vars, VAR_ID), -1, with = FALSE])
  w2_shared <- as.matrix(w2_dt[match(shared_vars, VAR_ID), -1, with = FALSE])
  rownames(w1_shared) <- shared_vars
  rownames(w2_shared) <- shared_vars

  w_cors <- data.table(
    Arm1_Cluster = matches$Arm1_Cluster,
    Arm2_Cluster = matches$Best_Arm2_Cluster,
    H_Correlation = matches$Correlation,
    W_Correlation = NA_real_
  )

  for (i in seq_len(nrow(w_cors))) {
    c1 <- sub(paste0("^", arm1_label, "_"), "", w_cors$Arm1_Cluster[i])
    c2 <- sub(paste0("^", arm2_label, "_"), "", w_cors$Arm2_Cluster[i])
    if (c1 %in% colnames(w1_shared) && c2 %in% colnames(w2_shared)) {
      w_cors$W_Correlation[i] <- cor(w1_shared[, c1], w2_shared[, c2],
                                     use = "complete.obs")
    }
  }

  cat("\nMatched cluster W correlations:\n")
  print(w_cors)
  fwrite(w_cors, file.path(output_dir, "cluster_W_correlations.tsv"), sep = "\t")
}


# ---- 5. Top traits per cluster (both arms) ----
cat("\n=== Top Traits Per Cluster ===\n")

get_top_traits <- function(h_signed, arm_name, n_top = 5) {
  results <- list()
  for (k in rownames(h_signed)) {
    loadings <- sort(abs(h_signed[k, ]), decreasing = TRUE)
    top <- head(loadings, n_top)
    signs <- ifelse(h_signed[k, names(top)] > 0, "+", "-")
    results[[length(results) + 1]] <- data.table(
      Arm = arm_name,
      Cluster = k,
      Rank = seq_along(top),
      Trait = names(top),
      Loading = as.numeric(top),
      Direction = signs
    )
  }
  rbindlist(results)
}

top1 <- get_top_traits(h1_signed, arm1_label)
top2 <- get_top_traits(h2_signed, arm2_label)
top_all <- rbind(top1, top2)

cat(sprintf("\n--- %s ---\n", arm1_label))
for (k in unique(top1$Cluster)) {
  sub <- top1[Cluster == k]
  cat(sprintf("  %s: %s\n", k,
              paste(sprintf("%s%s(%.2f)", sub$Direction, sub$Trait, sub$Loading),
                    collapse = ", ")))
}

cat(sprintf("\n--- %s ---\n", arm2_label))
for (k in unique(top2$Cluster)) {
  sub <- top2[Cluster == k]
  cat(sprintf("  %s: %s\n", k,
              paste(sprintf("%s%s(%.2f)", sub$Direction, sub$Trait, sub$Loading),
                    collapse = ", ")))
}

fwrite(top_all, file.path(output_dir, "top_traits_per_cluster.tsv"), sep = "\t")


# ---- 6. Shared trait loading scatter (matched clusters) ----
if (length(shared_traits) >= 2 && !is.null(matches)) {
  scatter_data <- list()
  for (i in seq_len(nrow(matches))) {
    c1_name <- sub(paste0("^", arm1_label, "_"), "", matches$Arm1_Cluster[i])
    c2_name <- sub(paste0("^", arm2_label, "_"), "", matches$Best_Arm2_Cluster[i])
    c1_label <- if (c1_name %in% names(arm1_labs)) arm1_labs[c1_name] else c1_name
    c2_label <- if (c2_name %in% names(arm2_labs)) arm2_labs[c2_name] else c2_name
    scatter_data[[i]] <- data.table(
      Trait = shared_traits,
      Arm1_Loading = h1_shared[c1_name, ],
      Arm2_Loading = h2_shared[c2_name, ],
      Pair = sprintf("%s vs %s", c1_label, c2_label)
    )
  }
  scatter_dt <- rbindlist(scatter_data)

  p_scatter <- ggplot(scatter_dt, aes(x = Arm1_Loading, y = Arm2_Loading)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_text(aes(label = Trait), size = 2, hjust = -0.1, vjust = -0.3,
              check_overlap = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~Pair, scales = "free") +
    labs(title = "Shared Trait Loadings: Matched Clusters",
         x = paste(arm1_label, "signed loading"),
         y = paste(arm2_label, "signed loading")) +
    theme_minimal(base_size = 11) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white"))

  n_pairs <- nrow(matches)
  ggsave(file.path(output_dir, "shared_trait_loadings.png"), p_scatter,
         width = min(16, max(6, n_pairs * 4)), height = max(5, ceiling(n_pairs / 3) * 4),
         dpi = 150)
  cat("\nWrote shared_trait_loadings.png\n")
}


# ---- 7. Write summary ----
sink(file.path(output_dir, "comparison_summary.txt"))
cat("=== bNMF Arm Comparison Summary ===\n\n")
cat(sprintf("Arm 1 (%s): %s\n", arm1_label, arm1_dir))
cat(sprintf("Arm 2 (%s): %s\n\n", arm2_label, arm2_dir))

cat(sprintf("Cluster count: %s K=%d, %s K=%d\n\n", arm1_label, k1, arm2_label, k2))

cat(sprintf("Trait count: %s=%d, %s=%d, shared=%d\n",
            arm1_label, length(traits1), arm2_label, length(traits2),
            length(shared_traits)))
if (length(shared_traits) > 0) {
  cat(sprintf("Shared traits: %s\n\n", paste(shared_traits, collapse = ", ")))
}

cat(sprintf("Variant count: %s=%d, %s=%d, shared=%d\n\n",
            arm1_label, length(vars1), arm2_label, length(vars2),
            length(shared_vars)))

if (!is.null(matches)) {
  cat("Cluster matching (by H correlation on shared traits):\n")
  print(matches)
}

cat("\nTop traits per cluster:\n")
print(top_all)
sink()

cat(sprintf("\nAll outputs written to %s/\n", output_dir))
cat("Done.\n")
