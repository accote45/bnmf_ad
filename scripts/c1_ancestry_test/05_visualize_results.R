#!/usr/bin/env Rscript
# 05_visualize_results.R
# Generate all 4 C1.1 output visualizations:
#   1. Per-cluster significance table
#   2. Null distribution histogram with cross-ancestry overlay
#   3. Matched cluster correspondence table
#   4. Effect size scatter plots (EUR vs AFR trait weights)
#
# Usage:
#   Rscript scripts/c1_ancestry_test/05_visualize_results.R
#   Rscript scripts/c1_ancestry_test/05_visualize_results.R --config config/c1_config.yaml

library(data.table)
library(yaml)
library(ggplot2)
library(patchwork)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/c1_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}

# --- Load config ---
cfg <- yaml::read_yaml(config_path)

project_root <- getwd()
results_dir  <- file.path(project_root, cfg$results_dir)
jaccard_dir  <- file.path(results_dir, "jaccard")
figures_dir  <- file.path(results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Source existing figure utilities
source(file.path(project_root, "scripts", "a1_analysis", "figure_utils.R"))
source(file.path(project_root, "scripts", "c1_ancestry_test", "jaccard_utils.R"))

traits <- cfg$gwas$traits
gtf_path <- cfg$gtf_path
gene_flank_kb <- if (is.null(cfg$jaccard$gene_flank_kb)) 50 else cfg$jaccard$gene_flank_kb

cat("=== C1.1 Visualization ===\n\n")


# =====================================================================
# Helper function (defined before use)
# =====================================================================

#' Extract SNP set from W data.table for gene mapping
get_cluster_snp_set_from_W <- function(W_dt, cluster_col, top_n = 100) {
  if (!cluster_col %in% colnames(W_dt)) return(character(0))
  weights <- W_dt[[cluster_col]]
  var_ids <- W_dt$VAR_ID
  names(weights) <- var_ids
  weights <- weights[weights > 0]
  if (length(weights) == 0) return(character(0))
  sorted <- sort(weights, decreasing = TRUE)
  names(head(sorted, min(top_n, length(sorted))))
}


# =====================================================================
# Load data
# =====================================================================

cat("--- Loading results ---\n")

null_dist <- fread(file.path(jaccard_dir, "null_distribution.tsv"))
cross_results <- fread(file.path(jaccard_dir, "cross_ancestry_results.tsv"))
cross_jaccards <- fread(file.path(jaccard_dir, "cross_ancestry_jaccards.tsv"))
per_cluster <- fread(file.path(jaccard_dir, "per_cluster_summary.tsv"))
k_values <- fread(file.path(jaccard_dir, "k_values.tsv"))

# Load precomputed SNP-to-gene mapping (from 04_compute_jaccard.R)
snp_to_gene_file <- file.path(jaccard_dir, "snp_to_gene_map.tsv")
snp_to_gene_map <- if (file.exists(snp_to_gene_file)) fread(snp_to_gene_file) else NULL

cat(sprintf("  Null Jaccard values: %d\n", nrow(null_dist)))
cat(sprintf("  Cross-ancestry results: %d rows\n", nrow(cross_results)))
cat(sprintf("  Per-cluster summary: %d clusters\n", nrow(per_cluster)))
if (!is.null(snp_to_gene_map)) {
  cat(sprintf("  SNP-to-gene map: %d mappings\n", nrow(snp_to_gene_map)))
}

# Load H and W matrices for trait comparisons
# EUR representative is chosen later (after Output 3 selects the median subsample)
# to ensure the scatter plot and correspondence table use the same subsample.
afr_h_file <- file.path(results_dir, "afr", "H_matrix_AFR.tsv")
afr_w_file <- file.path(results_dir, "afr", "W_matrix_AFR.tsv")
afr_H <- if (file.exists(afr_h_file)) fread(afr_h_file) else NULL
afr_W <- if (file.exists(afr_w_file)) fread(afr_w_file) else NULL

# Placeholder — will be set after Output 3 determines the representative subsample
eur_H <- NULL
eur_W <- NULL
eur_rep_label <- NULL
correspondence <- NULL
n_subs <- cfg$subsampling$n_subsamples


# =====================================================================
# Output 1: Per-cluster significance table
# =====================================================================

cat("\n--- Output 1: Significance table ---\n")

# The cross_results table already has per-comparison p-values
# Write it directly and also the aggregate summary
sig_table_file <- file.path(figures_dir, "output1_significance_table.tsv")
fwrite(cross_results, sig_table_file, sep = "\t")
cat(sprintf("  Written: %s\n", sig_table_file))

agg_file <- file.path(figures_dir, "output1_aggregate_pvalues.tsv")
fwrite(per_cluster, agg_file, sep = "\t")
cat(sprintf("  Written: %s\n", agg_file))


# =====================================================================
# Output 2: Null distribution histogram (2x3 grid: Jaccard + Spearman)
# =====================================================================

cat("\n--- Output 2: Null distribution histogram ---\n")

# Color palette for levels (shared between Jaccard and Spearman rows)
level_colors <- c(
  "Null"          = "grey50",
  "SNP-level"     = "steelblue4",
  "Gene-level"    = "seagreen4",
  "Pathway-level" = "mediumpurple4"
)

# Helper for dual histograms with legend
make_null_histogram <- function(null_dt, cross_dt, metric_col, plot_title, x_label,
                                cross_color = "steelblue4",
                                null_label = "Null EUR-EUR",
                                cross_label = "EUR-AFR") {
  # Prepare stacked data with a source column for legend
  null_plot <- copy(null_dt)[, source := null_label]
  cross_plot <- copy(cross_dt)[, source := cross_label]
  combined_data <- rbind(
    null_plot[, .(value = .SD[[metric_col]], source)],
    cross_plot[, .(value = .SD[[metric_col]], source)]
  )
  combined_data[, source := factor(source, levels = c(null_label, cross_label))]
  fill_vals <- setNames(c(level_colors["Null"], cross_color),
                        c(null_label, cross_label))
  color_vals <- setNames(c("grey40", cross_color),
                         c(null_label, cross_label))

  ggplot(combined_data, aes(x = value, fill = source, color = source)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 30, alpha = 0.6, position = "identity") +
    geom_rug(
      data = combined_data[source == cross_label],
      aes(x = value),
      color = cross_color, linewidth = 1, sides = "b", inherit.aes = FALSE
    ) +
    scale_fill_manual(values = fill_vals) +
    scale_color_manual(values = color_vals) +
    labs(x = x_label, y = "Density", title = plot_title, fill = NULL, color = NULL) +
    theme_big_text(base_size = 12) +
    theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
      legend.key.size = unit(0.4, "cm")
    )
}

# --- Top row: Jaccard distributions ---
p_snp_j <- make_null_histogram(null_dist, cross_jaccards, "snp_jaccard",
                                "SNP-level", "Jaccard Similarity",
                                cross_color = level_colors["SNP-level"],
                                null_label = "Null EUR-EUR Jaccard",
                                cross_label = "EUR-AFR Jaccard")
p_gene_j <- make_null_histogram(null_dist, cross_jaccards, "gene_jaccard",
                                 "Gene-level", "Jaccard Similarity",
                                 cross_color = level_colors["Gene-level"],
                                 null_label = "Null EUR-EUR Jaccard",
                                 cross_label = "EUR-AFR Jaccard")

has_pathway_data <- "pathway_jaccard" %in% colnames(null_dist) &&
                    "pathway_jaccard" %in% colnames(cross_jaccards)

if (has_pathway_data) {
  p_pathway_j <- make_null_histogram(null_dist, cross_jaccards, "pathway_jaccard",
                                      "Pathway-level", "Jaccard Similarity",
                                      cross_color = level_colors["Pathway-level"],
                                      null_label = "Null EUR-EUR Jaccard",
                                      cross_label = "EUR-AFR Jaccard")
} else {
  p_pathway_j <- plot_spacer()
}

# --- Bottom row: Spearman rho distributions ---
has_snp_spearman <- "snp_spearman" %in% colnames(null_dist) &&
                    "snp_spearman" %in% colnames(cross_jaccards)
has_gene_spearman <- "gene_spearman" %in% colnames(null_dist) &&
                     "gene_spearman" %in% colnames(cross_jaccards)
has_pathway_spearman <- "pathway_spearman" %in% colnames(null_dist) &&
                        "pathway_spearman" %in% colnames(cross_jaccards)

if (has_snp_spearman) {
  p_snp_s <- make_null_histogram(null_dist, cross_jaccards, "snp_spearman",
                                  "SNP-level", "Spearman rho",
                                  cross_color = level_colors["SNP-level"],
                                  null_label = "Null EUR-EUR Spearman",
                                  cross_label = "EUR-AFR Spearman")
} else {
  p_snp_s <- plot_spacer()
}

if (has_gene_spearman) {
  p_gene_s <- make_null_histogram(null_dist, cross_jaccards, "gene_spearman",
                                   "Gene-level", "Spearman rho",
                                   cross_color = level_colors["Gene-level"],
                                   null_label = "Null EUR-EUR Spearman",
                                   cross_label = "EUR-AFR Spearman")
} else {
  p_gene_s <- plot_spacer()
}

if (has_pathway_spearman) {
  p_pathway_s <- make_null_histogram(null_dist, cross_jaccards, "pathway_spearman",
                                      "Pathway-level", "Spearman rho",
                                      cross_color = level_colors["Pathway-level"],
                                      null_label = "Null EUR-EUR Spearman",
                                      cross_label = "EUR-AFR Spearman")
} else {
  p_pathway_s <- plot_spacer()
}

# Assemble 2x3 grid: row 1 = Jaccard, row 2 = Spearman
p_combined <- wrap_plots(list(p_snp_j, p_gene_j, p_pathway_j,
                              p_snp_s, p_gene_s, p_pathway_s),
                         nrow = 2, ncol = 3)

hist_file <- file.path(figures_dir, "output2_null_distribution")
ggsave(paste0(hist_file, ".png"), p_combined, width = 24, height = 14, dpi = 300)
cat(sprintf("  Written: %s.png\n", hist_file))


# =====================================================================
# Output 2b: Top-N pathway overlap (1 row x 3 col)
# =====================================================================

cat("\n--- Output 2b: Top-N pathway overlap ---\n")

has_topn_data <- all(c("pathway_top10_overlap", "pathway_top25_overlap", "pathway_top50_overlap") %in%
                       colnames(null_dist)) &&
                 all(c("pathway_top10_overlap", "pathway_top25_overlap", "pathway_top50_overlap") %in%
                       colnames(cross_jaccards))

if (has_topn_data) {
  p_top10 <- make_null_histogram(null_dist, cross_jaccards, "pathway_top10_overlap",
                                  "Pathway Top-10 Overlap", "Top-10 Overlap Fraction",
                                  cross_color = "darkorange3")
  p_top25b <- make_null_histogram(null_dist, cross_jaccards, "pathway_top25_overlap",
                                   "Pathway Top-25 Overlap", "Top-25 Overlap Fraction",
                                   cross_color = "goldenrod3")
  p_top50 <- make_null_histogram(null_dist, cross_jaccards, "pathway_top50_overlap",
                                  "Pathway Top-50 Overlap", "Top-50 Overlap Fraction",
                                  cross_color = "sienna3")
  p_topn <- wrap_plots(list(p_top10, p_top25b, p_top50), nrow = 1, ncol = 3)
  topn_file <- file.path(figures_dir, "output2b_topn_pathway_overlap")
  ggsave(paste0(topn_file, ".png"), p_topn, width = 24, height = 7, dpi = 300)
  cat(sprintf("  Written: %s.png\n", topn_file))
} else {
  cat("  Skipping: top-N pathway overlap columns not found in null/cross data.\n")
}


# =====================================================================
# Output 3: Matched cluster correspondence table
# =====================================================================

cat("\n--- Output 3: Correspondence table ---\n")

if (!is.null(afr_H) && nrow(cross_results) > 0) {

  # Use the most common matching across all EUR subsamples
  # Take the matching from the median subsample (by mean Jaccard)
  sub_mean_jaccard <- cross_results[, .(mean_j = mean(jaccard_sim)), by = eur_subsample]
  median_sub <- sub_mean_jaccard[which.min(abs(mean_j - median(mean_j))), eur_subsample]
  representative_matches <- cross_results[eur_subsample == median_sub]

  # Load EUR H and W matrices from the representative (median) subsample
  # so that Output 3 (correspondence), Output 4 (scatter), and Output 5
  # (heatmap) all use the same subsample whose matching is displayed.
  eur_rep_label <- median_sub
  eur_h_file <- file.path(results_dir, "subsamples", median_sub,
                           sprintf("H_matrix_%s.tsv", median_sub))
  eur_w_file <- file.path(results_dir, "subsamples", median_sub,
                           sprintf("W_matrix_%s.tsv", median_sub))
  eur_H <- if (file.exists(eur_h_file)) fread(eur_h_file) else NULL
  eur_W <- if (file.exists(eur_w_file)) fread(eur_w_file) else NULL
  cat(sprintf("  EUR representative subsample: %s (K=%d)\n",
              median_sub, if (!is.null(eur_H)) nrow(eur_H) else 0))

  correspondence <- data.table(
    EUR_cluster       = representative_matches$cluster_1,
    AFR_cluster       = representative_matches$cluster_2,
    snp_jaccard       = representative_matches$jaccard_sim,
    snp_p_value       = representative_matches$p_value,
    gene_jaccard      = representative_matches$gene_jaccard,
    gene_p_value      = representative_matches$gene_p_value,
    ld_ratio          = representative_matches$ld_ratio,
    trait_correlation  = NA_real_,
    shared_top_genes  = NA_character_,
    n_genes_eur       = NA_integer_,
    n_genes_afr       = NA_integer_
  )

  # Add pathway columns if available
  if ("pathway_jaccard" %in% colnames(representative_matches)) {
    correspondence[, pathway_jaccard := representative_matches$pathway_jaccard]
    correspondence[, pathway_p_value := representative_matches$pathway_p_value]
    correspondence[, pathway_ratio   := representative_matches$pathway_ratio]
  }

  # Add ranked pathway columns if available
  if ("pathway_spearman" %in% colnames(representative_matches)) {
    correspondence[, pathway_spearman      := representative_matches$pathway_spearman]
    correspondence[, pathway_top10_overlap  := representative_matches$pathway_top10_overlap]
    correspondence[, pathway_top25_overlap  := representative_matches$pathway_top25_overlap]
    correspondence[, pathway_top50_overlap  := representative_matches$pathway_top50_overlap]
  }

  # Compute trait weight correlations (from H matrices)
  for (row_i in seq_len(nrow(correspondence))) {
    eur_k <- correspondence$EUR_cluster[row_i]
    afr_k <- correspondence$AFR_cluster[row_i]

    if (eur_k %in% eur_H$Cluster && afr_k %in% afr_H$Cluster) {
      # Get common trait columns
      eur_cols <- setdiff(colnames(eur_H), "Cluster")
      afr_cols <- setdiff(colnames(afr_H), "Cluster")
      common_cols <- intersect(eur_cols, afr_cols)

      if (length(common_cols) >= 2) {
        eur_vec <- as.numeric(eur_H[Cluster == eur_k, ..common_cols])
        afr_vec <- as.numeric(afr_H[Cluster == afr_k, ..common_cols])
        correspondence$trait_correlation[row_i] <- cor(eur_vec, afr_vec)
      }
    }
  }

  # Map SNPs to genes using precomputed mapping from 04_compute_jaccard.R
  if (!is.null(eur_W) && !is.null(afr_W) && !is.null(snp_to_gene_map)) {
    for (row_i in seq_len(nrow(correspondence))) {
      eur_k <- correspondence$EUR_cluster[row_i]
      afr_k <- correspondence$AFR_cluster[row_i]

      # Get top SNPs per cluster
      eur_snps <- get_cluster_snp_set_from_W(eur_W, eur_k, top_n = cfg$jaccard$top_n)
      afr_snps <- get_cluster_snp_set_from_W(afr_W, afr_k, top_n = cfg$jaccard$top_n)

      # Gene-level intersection (captures same-locus-different-tag-SNP cases)
      eur_genes <- unique(snp_to_gene_map[VAR_ID %in% eur_snps, gene_name])
      afr_genes <- unique(snp_to_gene_map[VAR_ID %in% afr_snps, gene_name])
      shared_genes <- intersect(eur_genes, afr_genes)

      correspondence$n_genes_eur[row_i] <- length(eur_genes)
      correspondence$n_genes_afr[row_i] <- length(afr_genes)
      if (length(shared_genes) > 0) {
        correspondence$shared_top_genes[row_i] <- paste(
          head(shared_genes, 10), collapse = ", "
        )
      }
    }
  } else if (!is.null(eur_W) && !is.null(afr_W) && file.exists(gtf_path)) {
    # Fallback: use GTF directly if precomputed map not available
    tryCatch({
      gene_df <- parse_gtf_genes(gtf_path)
      for (row_i in seq_len(nrow(correspondence))) {
        eur_k <- correspondence$EUR_cluster[row_i]
        afr_k <- correspondence$AFR_cluster[row_i]
        eur_snps <- get_cluster_snp_set_from_W(eur_W, eur_k, top_n = cfg$jaccard$top_n)
        afr_snps <- get_cluster_snp_set_from_W(afr_W, afr_k, top_n = cfg$jaccard$top_n)
        shared_snps <- intersect(eur_snps, afr_snps)
        if (length(shared_snps) > 0) {
          snp_df <- data.table(VAR_ID = shared_snps)
          snp_df[, c("chr", "pos") := tstrsplit(VAR_ID, "_", keep = 1:2)]
          snp_df[, chr := as.integer(chr)]
          snp_df[, pos := as.integer(pos)]
          mapped <- map_snps_to_genes(snp_df, gene_df)
          if (nrow(mapped) > 0 && "gene_name" %in% colnames(mapped)) {
            correspondence$shared_top_genes[row_i] <- paste(
              head(unique(mapped$gene_name), 10), collapse = ", "
            )
          }
        }
      }
    }, error = function(e) {
      cat(sprintf("  Gene mapping failed: %s\n", e$message))
      cat("  Skipping shared gene annotation.\n")
    })
  }

  corr_file <- file.path(figures_dir, "output3_correspondence_table.tsv")
  fwrite(correspondence, corr_file, sep = "\t")
  cat(sprintf("  Written: %s\n", corr_file))
  print(correspondence)

} else {
  cat("  Skipping: missing H matrices or no cross-ancestry results.\n")
}


# =====================================================================
# Output 4: Effect size scatter plots (trait-correlation-based matching)
# =====================================================================
#
# Matching strategy:
#   The Jaccard-based matching from Output 3 can be degenerate when NMF
#   clusters share nearly all SNPs (small W matrix + large top_n makes
#   K1 and K2 SNP sets identical → arbitrary Hungarian assignment).
#
#   Instead, we match EUR↔AFR clusters by maximising the Pearson
#   correlation of their full trait-weight profiles from the H matrices.
#   This uses continuous net weights (pos − neg) across ALL traits
#   (GWAS + reference), not binary set membership.
#
# Each scatter point = one trait.
#   x = EUR cluster's net weight for that trait
#   y = AFR cluster's net weight for that trait
#   Points near the y = x diagonal → clusters share the same trait profile.

cat("\n--- Output 4: Effect size scatter plots ---\n")

if (!is.null(eur_H) && !is.null(afr_H) && nrow(cross_results) > 0) {

  # ---- 1. Identify all common traits in H matrices ----
  common_h_cols <- intersect(setdiff(colnames(eur_H), "Cluster"),
                             setdiff(colnames(afr_H), "Cluster"))
  all_h_traits <- extract_trait_names(c("Cluster", common_h_cols))

  # Helper: net weight vector (pos − neg) for one cluster across traits
  get_net_vec <- function(H_dt, cluster_label, trait_vec) {
    row <- H_dt[Cluster == cluster_label]
    h_cols <- colnames(row)
    sapply(trait_vec, function(tr) {
      pc <- paste0(tr, "_pos"); nc <- paste0(tr, "_neg")
      pv <- if (pc %in% h_cols) row[[pc]] else 0
      nv <- if (nc %in% h_cols) row[[nc]] else 0
      pv - nv
    })
  }

  eur_clusters <- eur_H$Cluster
  afr_clusters <- afr_H$Cluster

  # ---- 2. Precompute all net-weight vectors (used by correlation + scatter) ----
  eur_net_vecs <- setNames(
    lapply(eur_clusters, function(k) get_net_vec(eur_H, k, all_h_traits)),
    eur_clusters
  )
  afr_net_vecs <- setNames(
    lapply(afr_clusters, function(k) get_net_vec(afr_H, k, all_h_traits)),
    afr_clusters
  )

  # ---- 3. Trait-correlation matrix: cor(EUR_i, AFR_j) ----
  cor_mat <- matrix(NA_real_,
                    nrow = length(eur_clusters), ncol = length(afr_clusters),
                    dimnames = list(eur_clusters, afr_clusters))
  for (i in seq_along(eur_clusters)) {
    for (j in seq_along(afr_clusters)) {
      cor_mat[i, j] <- cor(eur_net_vecs[[i]], afr_net_vecs[[j]])
    }
  }

  cat("  Trait-weight correlation matrix (EUR rows \u00d7 AFR cols):\n")
  print(round(cor_mat, 3))

  # ---- 4. Hungarian matching to maximise total correlation ----
  # solve_LSAP requires a square, non-negative matrix; pad if K differs
  # and shift cor ∈ [-1,1] to [0,2]
  shifted <- cor_mat + 1
  n_eur <- length(eur_clusters)
  n_afr <- length(afr_clusters)
  if (n_eur != n_afr) {
    max_k <- max(n_eur, n_afr)
    padded <- matrix(0, nrow = max_k, ncol = max_k)
    padded[seq_len(n_eur), seq_len(n_afr)] <- shifted
    assignment <- clue::solve_LSAP(padded, maximum = TRUE)[seq_len(n_eur)]
  } else {
    assignment <- clue::solve_LSAP(shifted, maximum = TRUE)
  }

  trait_matches <- data.table(
    EUR_cluster = eur_clusters,
    AFR_cluster = afr_clusters[assignment],
    trait_cor   = sapply(seq_along(eur_clusters),
                         function(i) cor_mat[i, assignment[i]])
  )

  cat("\n  Trait-correlation matching (used for scatter):\n")
  print(trait_matches)

  # Compare with Jaccard-based matching from Output 3
  if (!is.null(correspondence) && nrow(correspondence) > 0) {
    jaccard_str <- paste(correspondence$EUR_cluster, "\u2194",
                         correspondence$AFR_cluster, collapse = ", ")
    trait_str   <- paste(trait_matches$EUR_cluster, "\u2194",
                         trait_matches$AFR_cluster, collapse = ", ")
    cat(sprintf("\n  Jaccard matching  (Output 3): %s\n", jaccard_str))
    cat(sprintf("  Trait-cor matching (Output 4): %s\n", trait_str))
    if (jaccard_str != trait_str) {
      cat("  \u26a0 Matchings differ \u2014 Jaccard was likely degenerate.\n")
    }
  }

  # ---- 5. Generate scatter plots ----
  plots <- list()
  gwas_traits <- cfg$gwas$traits
  use_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

  for (row_i in seq_len(nrow(trait_matches))) {
    eur_k <- trait_matches$EUR_cluster[row_i]
    afr_k <- trait_matches$AFR_cluster[row_i]
    tc    <- trait_matches$trait_cor[row_i]

    eur_net <- eur_net_vecs[[eur_k]]
    afr_net <- afr_net_vecs[[afr_k]]

    scatter_df <- data.frame(
      trait      = all_h_traits,
      EUR        = eur_net,
      AFR        = afr_net,
      point_type = factor(ifelse(all_h_traits %in% gwas_traits,
                                 "GWAS trait", "Reference"),
                          levels = c("GWAS trait", "Reference")),
      stringsAsFactors = FALSE
    )

    max_abs <- max(abs(c(eur_net, afr_net)), na.rm = TRUE) * 1.2
    if (is.infinite(max_abs) || max_abs == 0) max_abs <- 1

    p <- ggplot(scatter_df, aes(x = EUR, y = AFR)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey50", linewidth = 0.5) +
      geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
      geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
      geom_point(aes(color = point_type, size = point_type)) +
      scale_color_manual(values = c("GWAS trait" = "firebrick3",
                                    "Reference"  = "steelblue"),
                         name = NULL) +
      scale_size_manual(values = c("GWAS trait" = 5, "Reference" = 3),
                        guide = "none") +
      coord_fixed(xlim = c(-max_abs, max_abs), ylim = c(-max_abs, max_abs)) +
      labs(
        title = sprintf("%s (EUR) \u2194 %s (AFR)  [r = %.2f]", eur_k, afr_k, tc),
        x = "EUR Net Trait Weight (pos \u2212 neg)",
        y = "AFR Net Trait Weight (pos \u2212 neg)"
      ) +
      theme_big_text(base_size = 12) +
      theme(legend.position = c(0.02, 0.98),
            legend.justification = c(0, 1),
            legend.background = element_rect(fill = alpha("white", 0.8), color = NA))

    if (use_ggrepel) {
      p <- p + ggrepel::geom_text_repel(
        aes(label = trait), size = 3.5,
        max.overlaps = Inf, seed = 42, segment.color = "grey60"
      )
    } else {
      p <- p + geom_text(aes(label = trait), vjust = -1.2, size = 3.5)
    }

    plots[[row_i]] <- p
  }

  if (length(plots) > 0) {
    combined <- wrap_plots(plots, ncol = length(plots))
    scatter_file <- file.path(figures_dir, "output4_effect_size_scatter")

    plot_width  <- 9 * length(plots)
    plot_height <- 9

    ggsave(paste0(scatter_file, ".png"), combined,
           width = plot_width, height = plot_height, dpi = 300)
    cat(sprintf("  Written: %s.png (%d panels)\n",
                scatter_file, length(plots)))
  } else {
    cat("  No valid cluster pairs for scatter plots.\n")
  }

} else {
  cat("  Skipping: missing H matrices or no cross-ancestry results.\n")
}


# =====================================================================
# Output 5: Cluster-trait heatmap (EUR K1/K2 + AFR K1/K2)
# =====================================================================

cat("\n--- Output 5: Cluster-trait heatmap ---\n")

if (!is.null(eur_H) && !is.null(afr_H) && nrow(eur_H) >= 2 && nrow(afr_H) >= 2) {

  # Extract trait names from H matrix columns (TRAIT_pos / TRAIT_neg pairs)
  trait_names <- extract_trait_names(colnames(eur_H))

  # Compute net weights (pos - neg) for each ancestry
  compute_net_weights <- function(H_dt, ancestry_label) {
    rows <- list()
    for (ri in seq_len(nrow(H_dt))) {
      cluster <- H_dt$Cluster[ri]
      net <- sapply(trait_names, function(tr) {
        pos_val <- if (paste0(tr, "_pos") %in% colnames(H_dt)) H_dt[[paste0(tr, "_pos")]][ri] else 0
        neg_val <- if (paste0(tr, "_neg") %in% colnames(H_dt)) H_dt[[paste0(tr, "_neg")]][ri] else 0
        pos_val - neg_val
      })
      rows[[ri]] <- data.table(
        ancestry = ancestry_label,
        cluster  = cluster,
        label    = sprintf("%s %s", ancestry_label, cluster),
        trait    = trait_names,
        weight   = net
      )
    }
    rbindlist(rows)
  }

  eur_label <- sprintf("EUR (%s)", eur_rep_label)
  heatmap_dt <- rbind(
    compute_net_weights(eur_H, eur_label),
    compute_net_weights(afr_H, "AFR")
  )

  # Order rows: EUR K1, EUR K2, AFR K1, AFR K2
  heatmap_dt[, label := factor(label, levels = rev(unique(label)))]

  p_heatmap <- ggplot(heatmap_dt, aes(x = trait, y = label, fill = weight)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.1f", weight)), size = 2.2, color = "black") +
    scale_fill_gradient2(
      low = "steelblue4", mid = "white", high = "firebrick3",
      midpoint = 0, name = "Net Weight\n(pos \u2212 neg)"
    ) +
    theme_big_text(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y  = element_text(size = 12),
      panel.grid   = element_blank()
    ) +
    labs(
      x     = "Trait",
      y     = NULL,
      title = "Cluster-Trait Weight Heatmap (EUR vs AFR)"
    ) +
    coord_fixed(ratio = 0.8)

  heatmap_file <- file.path(figures_dir, "output5_cluster_trait_heatmap")
  ggsave(paste0(heatmap_file, ".png"), p_heatmap, width = 14, height = 5, dpi = 300)
  cat(sprintf("  Written: %s.png\n", heatmap_file))

} else {
  cat("  Skipping: need K>=2 for both EUR and AFR H matrices.\n")
  if (!is.null(eur_H)) cat(sprintf("    EUR K=%d\n", nrow(eur_H)))
  if (!is.null(afr_H)) cat(sprintf("    AFR K=%d\n", nrow(afr_H)))
}


# =====================================================================
# Summary
# =====================================================================

cat("\n=== Visualization complete ===\n")
cat(sprintf("Figures in: %s\n", figures_dir))
cat("Files:\n")
cat("  output1_significance_table.tsv   - Per-comparison p-values (SNP + gene)\n")
cat("  output1_aggregate_pvalues.tsv    - Per-cluster aggregate p-values (SNP + gene)\n")
cat("  output2_null_distribution.png     - SNP + gene + pathway Jaccard null distribution histograms\n")
cat("  output3_correspondence_table.tsv - Matched cluster table with SNP/gene/pathway Jaccard + genes\n")
cat("  output4_effect_size_scatter.png   - EUR vs AFR trait weight scatter\n")
cat("  output5_cluster_trait_heatmap.png - EUR + AFR cluster-trait weight heatmap\n")
