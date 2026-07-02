#!/usr/bin/env Rscript
# 04_compute_jaccard.R
# Core C1.1 analysis: compute Jaccard distances between bNMF clusters,
# use Hungarian algorithm for optimal cluster matching, build within-EUR
# null distribution, and test cross-ancestry (EUR vs AFR) significance.
#
# Usage:
#   Rscript scripts/c1_ancestry_test/04_compute_jaccard.R
#   Rscript scripts/c1_ancestry_test/04_compute_jaccard.R --config config/c1_config.yaml

library(data.table)
library(yaml)

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
dir.create(jaccard_dir, recursive = TRUE, showWarnings = FALSE)

top_n           <- cfg$jaccard$top_n              # 100
n_subsamples    <- cfg$subsampling$n_subsamples   # 10
membership_method <- cfg$jaccard$membership_method # "top_n"
gene_flank_kb   <- if (is.null(cfg$jaccard$gene_flank_kb)) 50 else cfg$jaccard$gene_flank_kb
gene_flank_bp   <- gene_flank_kb * 1000
gtf_path        <- cfg$gtf_path
pathway_gmt_path <- cfg$jaccard$pathway_gmt
pathway_fdr      <- if (is.null(cfg$jaccard$pathway_fdr)) 0.05 else cfg$jaccard$pathway_fdr

# Source existing utilities
source(file.path(project_root, "scripts", "a1_analysis", "figure_utils.R"))
source(file.path(project_root, "scripts", "c1_ancestry_test", "jaccard_utils.R"))

cat("=== C1.1 Jaccard Analysis ===\n\n")


# =====================================================================
# SECTION 2: Load W matrices from all bNMF runs
# =====================================================================

cat("--- Loading W matrices ---\n")

# EUR subsamples
eur_W <- list()
eur_K <- integer(0)

for (i in seq_len(n_subsamples)) {
  sub_label <- sprintf("sub_%02d", i)
  w_file <- file.path(results_dir, "subsamples", sub_label,
                       sprintf("W_matrix_%s.tsv", sub_label))

  if (!file.exists(w_file)) {
    cat(sprintf("  WARNING: W matrix not found for %s: %s\n", sub_label, w_file))
    next
  }

  W <- fread(w_file)
  eur_W[[sub_label]] <- W
  k <- ncol(W) - 1  # subtract VAR_ID column
  eur_K <- c(eur_K, k)
  cat(sprintf("  %s: K=%d, %d variants\n", sub_label, k, nrow(W)))
}

cat(sprintf("\n  EUR K distribution: %s\n", paste(eur_K, collapse = ", ")))
cat(sprintf("  EUR modal K: %d\n",
            as.integer(names(sort(table(eur_K), decreasing = TRUE))[1])))

# AFR
afr_w_file <- file.path(results_dir, "afr", "W_matrix_AFR.tsv")
if (!file.exists(afr_w_file)) {
  stop(sprintf("AFR W matrix not found: %s", afr_w_file))
}
afr_W <- fread(afr_w_file)
afr_K <- ncol(afr_W) - 1
cat(sprintf("\n  AFR: K=%d, %d variants\n", afr_K, nrow(afr_W)))

if (afr_K == 0) {
  cat("\nWARNING: AFR has K=0. Cross-ancestry comparison not possible.\n")
  cat("Writing empty results.\n")
  fwrite(data.table(note = "AFR K=0, no clusters to compare"),
         file.path(jaccard_dir, "cross_ancestry_results.tsv"), sep = "\t")
  quit(save = "no", status = 0)
}


# =====================================================================
# SECTION 2b: Build precomputed SNP-to-gene mapping
# =====================================================================

cat("\n--- Building SNP-to-gene mapping ---\n")
cat(sprintf("  GTF: %s\n", gtf_path))
cat(sprintf("  Flanking window: %d kb (%d bp)\n", gene_flank_kb, gene_flank_bp))

gene_dt <- parse_gtf_genes_dt(gtf_path)
cat(sprintf("  Protein-coding genes: %d\n", nrow(gene_dt)))

# Collect all VAR_IDs across all W matrices
all_var_ids <- unique(c(
  unlist(lapply(eur_W, function(W) W$VAR_ID)),
  afr_W$VAR_ID
))
cat(sprintf("  Unique VAR_IDs across all runs: %d\n", length(all_var_ids)))

snp_to_gene_map <- build_snp_to_gene_map(all_var_ids, gene_dt, flank_bp = gene_flank_bp)
cat(sprintf("  SNPs mapped to at least one gene: %d / %d (%.1f%%)\n",
            length(unique(snp_to_gene_map$VAR_ID)),
            length(all_var_ids),
            100 * length(unique(snp_to_gene_map$VAR_ID)) / length(all_var_ids)))
cat(sprintf("  Unique genes in mapping: %d\n",
            length(unique(snp_to_gene_map$gene_name))))

fwrite(snp_to_gene_map, file.path(jaccard_dir, "snp_to_gene_map.tsv"), sep = "\t")


# =====================================================================
# SECTION 2c: Load pathway data and build background
# =====================================================================

cat("\n--- Loading pathway data ---\n")

pathway_map <- NULL
background_gene_ids <- NULL

if (!is.null(pathway_gmt_path) && file.exists(pathway_gmt_path)) {
  pathway_map <- parse_gmt(pathway_gmt_path)
  cat(sprintf("  GMT file: %s\n", pathway_gmt_path))
  cat(sprintf("  Pathways loaded: %d\n", length(unique(pathway_map$pathway))))
  cat(sprintf("  Unique genes in GMT: %d\n", length(unique(pathway_map$gene_id))))

  # Background = all Ensembl gene IDs in the SNP-to-gene mapping
  background_gene_ids <- unique(snp_to_gene_map$gene_id)
  gmt_genes <- unique(pathway_map$gene_id)
  overlap <- length(intersect(background_gene_ids, gmt_genes))
  cat(sprintf("  Background genes (from SNP mapping): %d\n", length(background_gene_ids)))
  cat(sprintf("  Overlap (background ∩ GMT): %d (%.1f%% of background)\n",
              overlap, 100 * overlap / length(background_gene_ids)))
  cat(sprintf("  FDR threshold: %g\n", pathway_fdr))
} else {
  cat("  WARNING: Pathway GMT file not found or not configured. Skipping pathway-level analysis.\n")
  if (!is.null(pathway_gmt_path)) cat(sprintf("  Path: %s\n", pathway_gmt_path))
}

has_pathway <- !is.null(pathway_map) && !is.null(background_gene_ids)


# =====================================================================
# SECTION 3: Build within-EUR null distribution
# =====================================================================

cat("\n--- Building within-EUR null distribution ---\n")

# All C(n,2) pairwise EUR-EUR comparisons
sub_names <- names(eur_W)
pairs <- combn(sub_names, 2, simplify = FALSE)
cat(sprintf("  Number of EUR-EUR pairs: %d\n", length(pairs)))

null_jaccards         <- numeric()
null_gene_jaccards    <- numeric()
null_snp_spearman     <- numeric()
null_gene_spearman    <- numeric()
null_pathway_jaccards <- numeric()
null_pathway_spearman <- numeric()
null_pathway_top10    <- numeric()
null_pathway_top25    <- numeric()
null_pathway_top50    <- numeric()
pairwise_results      <- list()

for (pair in pairs) {
  pair_label <- paste(pair, collapse = "_vs_")
  result <- hungarian_match_clusters(eur_W[[pair[1]]], eur_W[[pair[2]]],
                                     top_n = top_n)
  pairwise_results[[pair_label]] <- result

  if (nrow(result$matching) > 0) {
    null_jaccards <- c(null_jaccards, result$matching$jaccard_sim)

    # Gene-level Jaccard for each matched pair
    gene_jaccards_this_pair <- sapply(seq_len(nrow(result$matching)), function(i) {
      c1 <- result$matching$cluster_1[i]
      c2 <- result$matching$cluster_2[i]
      gene_level_jaccard(result$snp_sets_1[[c1]], result$snp_sets_2[[c2]],
                         snp_to_gene_map)
    })
    result$matching$gene_jaccard <- gene_jaccards_this_pair
    null_gene_jaccards <- c(null_gene_jaccards, gene_jaccards_this_pair)

    # SNP-level and gene-level Spearman for each matched pair
    snp_spearman_this_pair <- sapply(seq_len(nrow(result$matching)), function(i) {
      c1 <- result$matching$cluster_1[i]
      c2 <- result$matching$cluster_2[i]
      snp_level_spearman(eur_W[[pair[1]]], eur_W[[pair[2]]], c1, c2)
    })
    gene_spearman_this_pair <- sapply(seq_len(nrow(result$matching)), function(i) {
      c1 <- result$matching$cluster_1[i]
      c2 <- result$matching$cluster_2[i]
      gene_level_spearman(eur_W[[pair[1]]], eur_W[[pair[2]]], c1, c2, snp_to_gene_map)
    })
    result$matching$snp_spearman  <- snp_spearman_this_pair
    result$matching$gene_spearman <- gene_spearman_this_pair
    null_snp_spearman  <- c(null_snp_spearman, snp_spearman_this_pair)
    null_gene_spearman <- c(null_gene_spearman, gene_spearman_this_pair)

    # Pathway-level Jaccard for each matched pair
    if (has_pathway) {
      pathway_jaccards_this_pair <- sapply(seq_len(nrow(result$matching)), function(i) {
        c1 <- result$matching$cluster_1[i]
        c2 <- result$matching$cluster_2[i]
        pathway_level_jaccard(result$snp_sets_1[[c1]], result$snp_sets_2[[c2]],
                              snp_to_gene_map, pathway_map, background_gene_ids,
                              pathway_fdr)
      })
      result$matching$pathway_jaccard <- pathway_jaccards_this_pair
      null_pathway_jaccards <- c(null_pathway_jaccards, pathway_jaccards_this_pair)

      # Ranked pathway metrics for each matched pair
      ranked_this_pair <- lapply(seq_len(nrow(result$matching)), function(i) {
        c1 <- result$matching$cluster_1[i]
        c2 <- result$matching$cluster_2[i]
        ranked_pathway_similarity(result$snp_sets_1[[c1]], result$snp_sets_2[[c2]],
                                   snp_to_gene_map, pathway_map, background_gene_ids)
      })
      result$matching$pathway_spearman      <- sapply(ranked_this_pair, `[[`, "pathway_spearman")
      result$matching$pathway_top10_overlap  <- sapply(ranked_this_pair, `[[`, "pathway_top10_overlap")
      result$matching$pathway_top25_overlap  <- sapply(ranked_this_pair, `[[`, "pathway_top25_overlap")
      result$matching$pathway_top50_overlap  <- sapply(ranked_this_pair, `[[`, "pathway_top50_overlap")

      null_pathway_spearman <- c(null_pathway_spearman, result$matching$pathway_spearman)
      null_pathway_top10    <- c(null_pathway_top10, result$matching$pathway_top10_overlap)
      null_pathway_top25    <- c(null_pathway_top25, result$matching$pathway_top25_overlap)
      null_pathway_top50    <- c(null_pathway_top50, result$matching$pathway_top50_overlap)
    }

    pairwise_results[[pair_label]] <- result

    cat(sprintf("  %s: %d matched, SNP Jaccard=%.3f, Gene Jaccard=%.3f%s\n",
                pair_label, nrow(result$matching),
                mean(result$matching$jaccard_sim),
                mean(gene_jaccards_this_pair),
                if (has_pathway) sprintf(", Pathway Jaccard=%.3f, Spearman=%.3f",
                  mean(result$matching$pathway_jaccard),
                  mean(result$matching$pathway_spearman, na.rm = TRUE)) else ""))
    cat(sprintf("    SNP Spearman=%.3f, Gene Spearman=%.3f\n",
                mean(snp_spearman_this_pair, na.rm = TRUE),
                mean(gene_spearman_this_pair, na.rm = TRUE)))
  } else {
    cat(sprintf("  %s: no clusters matched\n", pair_label))
  }
}

cat(sprintf("\n  SNP-level null: %d values, mean=%.4f, median=%.4f, sd=%.4f\n",
            length(null_jaccards),
            mean(null_jaccards), median(null_jaccards), sd(null_jaccards)))
cat(sprintf("  SNP-level null range: [%.4f, %.4f]\n",
            min(null_jaccards), max(null_jaccards)))
cat(sprintf("  Gene-level null: %d values, mean=%.4f, median=%.4f, sd=%.4f\n",
            length(null_gene_jaccards),
            mean(null_gene_jaccards), median(null_gene_jaccards), sd(null_gene_jaccards)))
if (length(null_gene_jaccards) < 2 || isTRUE(sd(null_gene_jaccards) == 0)) {
  cat("  WARNING: Gene-level null has zero or undefined variance. Gene p-values will be uninformative.\n")
}
if (has_pathway) {
  cat(sprintf("  Pathway-level null: %d values, mean=%.4f, median=%.4f, sd=%.4f\n",
              length(null_pathway_jaccards),
              mean(null_pathway_jaccards), median(null_pathway_jaccards),
              sd(null_pathway_jaccards)))
  if (length(null_pathway_jaccards) < 2 || isTRUE(sd(null_pathway_jaccards) == 0)) {
    cat("  WARNING: Pathway-level null has zero or undefined variance. Pathway p-values will be uninformative.\n")
  }
}

# Save null distribution
null_dt <- data.table(snp_jaccard = null_jaccards, gene_jaccard = null_gene_jaccards,
                      snp_spearman = null_snp_spearman, gene_spearman = null_gene_spearman,
                      type = "EUR_EUR_null")
if (has_pathway) {
  null_dt[, pathway_jaccard := null_pathway_jaccards]
  null_dt[, pathway_spearman := null_pathway_spearman]
  null_dt[, pathway_top10_overlap := null_pathway_top10]
  null_dt[, pathway_top25_overlap := null_pathway_top25]
  null_dt[, pathway_top50_overlap := null_pathway_top50]
}
fwrite(null_dt, file.path(jaccard_dir, "null_distribution.tsv"), sep = "\t")


# =====================================================================
# SECTION 4: Cross-ancestry test (EUR vs AFR)
# =====================================================================

cat("\n--- Cross-ancestry comparison (EUR vs AFR) ---\n")

cross_results_list     <- list()
cross_jaccards         <- numeric()
cross_gene_jaccards    <- numeric()
cross_snp_spearman     <- numeric()
cross_gene_spearman    <- numeric()
cross_pathway_jaccards <- numeric()
cross_pathway_spearman <- numeric()
cross_pathway_top10    <- numeric()
cross_pathway_top25    <- numeric()
cross_pathway_top50    <- numeric()

for (sub_name in sub_names) {
  result <- hungarian_match_clusters(eur_W[[sub_name]], afr_W, top_n = top_n)

  if (nrow(result$matching) > 0) {
    # SNP-level p-values
    result$matching$p_value <- sapply(result$matching$jaccard_sim, function(obs) {
      # p = fraction of null values <= observed
      # (low Jaccard = more different = low p-value)
      mean(null_jaccards <= obs)
    })

    # Gene-level Jaccard for each matched pair
    result$matching$gene_jaccard <- sapply(seq_len(nrow(result$matching)), function(i) {
      c1 <- result$matching$cluster_1[i]
      c2 <- result$matching$cluster_2[i]
      gene_level_jaccard(result$snp_sets_1[[c1]], result$snp_sets_2[[c2]],
                         snp_to_gene_map)
    })

    # Gene-level p-values
    result$matching$gene_p_value <- sapply(result$matching$gene_jaccard, function(obs) {
      mean(null_gene_jaccards <= obs)
    })

    # LD ratio: gene_jaccard / snp_jaccard
    result$matching$ld_ratio <- ifelse(
      result$matching$jaccard_sim > 0,
      result$matching$gene_jaccard / result$matching$jaccard_sim,
      NA_real_
    )

    # SNP-level and gene-level Spearman for each matched pair
    result$matching$snp_spearman <- sapply(seq_len(nrow(result$matching)), function(i) {
      c1 <- result$matching$cluster_1[i]
      c2 <- result$matching$cluster_2[i]
      snp_level_spearman(eur_W[[sub_name]], afr_W, c1, c2)
    })
    result$matching$gene_spearman <- sapply(seq_len(nrow(result$matching)), function(i) {
      c1 <- result$matching$cluster_1[i]
      c2 <- result$matching$cluster_2[i]
      gene_level_spearman(eur_W[[sub_name]], afr_W, c1, c2, snp_to_gene_map)
    })
    # Spearman p-values
    result$matching$snp_spearman_p <- sapply(result$matching$snp_spearman, function(obs) {
      if (is.na(obs)) return(NA_real_)
      mean(null_snp_spearman <= obs, na.rm = TRUE)
    })
    result$matching$gene_spearman_p <- sapply(result$matching$gene_spearman, function(obs) {
      if (is.na(obs)) return(NA_real_)
      mean(null_gene_spearman <= obs, na.rm = TRUE)
    })

    # Pathway-level Jaccard for each matched pair
    if (has_pathway) {
      result$matching$pathway_jaccard <- sapply(seq_len(nrow(result$matching)), function(i) {
        c1 <- result$matching$cluster_1[i]
        c2 <- result$matching$cluster_2[i]
        pathway_level_jaccard(result$snp_sets_1[[c1]], result$snp_sets_2[[c2]],
                              snp_to_gene_map, pathway_map, background_gene_ids,
                              pathway_fdr)
      })

      result$matching$pathway_p_value <- sapply(result$matching$pathway_jaccard, function(obs) {
        mean(null_pathway_jaccards <= obs)
      })

      result$matching$pathway_ratio <- ifelse(
        result$matching$jaccard_sim > 0,
        result$matching$pathway_jaccard / result$matching$jaccard_sim,
        NA_real_
      )

      # Ranked pathway metrics for each matched pair
      ranked_cross <- lapply(seq_len(nrow(result$matching)), function(i) {
        c1 <- result$matching$cluster_1[i]
        c2 <- result$matching$cluster_2[i]
        ranked_pathway_similarity(result$snp_sets_1[[c1]], result$snp_sets_2[[c2]],
                                   snp_to_gene_map, pathway_map, background_gene_ids)
      })
      result$matching$pathway_spearman      <- sapply(ranked_cross, `[[`, "pathway_spearman")
      result$matching$pathway_top10_overlap  <- sapply(ranked_cross, `[[`, "pathway_top10_overlap")
      result$matching$pathway_top25_overlap  <- sapply(ranked_cross, `[[`, "pathway_top25_overlap")
      result$matching$pathway_top50_overlap  <- sapply(ranked_cross, `[[`, "pathway_top50_overlap")

      # Empirical p-values for ranked metrics
      result$matching$pathway_spearman_p <- sapply(result$matching$pathway_spearman, function(obs) {
        if (is.na(obs)) return(NA_real_)
        mean(null_pathway_spearman <= obs, na.rm = TRUE)
      })
      result$matching$pathway_top10_p <- sapply(result$matching$pathway_top10_overlap, function(obs) {
        if (is.na(obs)) return(NA_real_)
        mean(null_pathway_top10 <= obs, na.rm = TRUE)
      })
      result$matching$pathway_top25_p <- sapply(result$matching$pathway_top25_overlap, function(obs) {
        if (is.na(obs)) return(NA_real_)
        mean(null_pathway_top25 <= obs, na.rm = TRUE)
      })
      result$matching$pathway_top50_p <- sapply(result$matching$pathway_top50_overlap, function(obs) {
        if (is.na(obs)) return(NA_real_)
        mean(null_pathway_top50 <= obs, na.rm = TRUE)
      })
    }

    result$matching$eur_subsample <- sub_name
    cross_results_list[[sub_name]] <- result

    cross_jaccards      <- c(cross_jaccards, result$matching$jaccard_sim)
    cross_gene_jaccards <- c(cross_gene_jaccards, result$matching$gene_jaccard)
    cross_snp_spearman  <- c(cross_snp_spearman, result$matching$snp_spearman)
    cross_gene_spearman <- c(cross_gene_spearman, result$matching$gene_spearman)
    if (has_pathway) {
      cross_pathway_jaccards <- c(cross_pathway_jaccards, result$matching$pathway_jaccard)
      cross_pathway_spearman <- c(cross_pathway_spearman, result$matching$pathway_spearman)
      cross_pathway_top10    <- c(cross_pathway_top10, result$matching$pathway_top10_overlap)
      cross_pathway_top25    <- c(cross_pathway_top25, result$matching$pathway_top25_overlap)
      cross_pathway_top50    <- c(cross_pathway_top50, result$matching$pathway_top50_overlap)
    }

    cat(sprintf("  %s vs AFR: %d matched, SNP J=%.3f (p=%.4f), Gene J=%.3f (p=%.4f)%s\n",
                sub_name, nrow(result$matching),
                mean(result$matching$jaccard_sim),
                mean(result$matching$p_value),
                mean(result$matching$gene_jaccard),
                mean(result$matching$gene_p_value),
                if (has_pathway) sprintf(", Pathway J=%.3f (p=%.4f), Spearman=%.3f",
                  mean(result$matching$pathway_jaccard),
                  mean(result$matching$pathway_p_value),
                  mean(result$matching$pathway_spearman, na.rm = TRUE)) else ""))
  } else {
    cat(sprintf("  %s vs AFR: no clusters matched\n", sub_name))
  }
}

# Combine all cross-ancestry results
cross_results <- rbindlist(
  lapply(cross_results_list, function(r) as.data.table(r$matching)),
  fill = TRUE
)

cat(sprintf("\n  Cross-ancestry SNP Jaccard: %d values, mean=%.4f, median=%.4f\n",
            length(cross_jaccards),
            mean(cross_jaccards), median(cross_jaccards)))
cat(sprintf("  Cross-ancestry Gene Jaccard: %d values, mean=%.4f, median=%.4f\n",
            length(cross_gene_jaccards),
            mean(cross_gene_jaccards), median(cross_gene_jaccards)))
if (has_pathway) {
  cat(sprintf("  Cross-ancestry Pathway Jaccard: %d values, mean=%.4f, median=%.4f\n",
              length(cross_pathway_jaccards),
              mean(cross_pathway_jaccards), median(cross_pathway_jaccards)))
}

# Save cross-ancestry results
fwrite(cross_results, file.path(jaccard_dir, "cross_ancestry_results.tsv"),
       sep = "\t")

# Save cross-ancestry Jaccard values
cross_jacc_dt <- data.table(snp_jaccard = cross_jaccards,
                            gene_jaccard = cross_gene_jaccards,
                            snp_spearman = cross_snp_spearman,
                            gene_spearman = cross_gene_spearman,
                            type = "EUR_AFR_cross")
if (has_pathway) {
  cross_jacc_dt[, pathway_jaccard := cross_pathway_jaccards]
  cross_jacc_dt[, pathway_spearman := cross_pathway_spearman]
  cross_jacc_dt[, pathway_top10_overlap := cross_pathway_top10]
  cross_jacc_dt[, pathway_top25_overlap := cross_pathway_top25]
  cross_jacc_dt[, pathway_top50_overlap := cross_pathway_top50]
}
fwrite(cross_jacc_dt, file.path(jaccard_dir, "cross_ancestry_jaccards.tsv"),
       sep = "\t")


# =====================================================================
# SECTION 5: Aggregate per-cluster summary
# =====================================================================

cat("\n--- Per-cluster aggregate summary ---\n")

if (nrow(cross_results) > 0) {
  # For each AFR cluster, aggregate across EUR subsamples
  agg <- cross_results[, .(
    n_comparisons       = .N,
    # SNP-level
    median_snp_jaccard  = median(jaccard_sim, na.rm = TRUE),
    mean_snp_jaccard    = mean(jaccard_sim, na.rm = TRUE),
    sd_snp_jaccard      = sd(jaccard_sim, na.rm = TRUE),
    median_snp_p        = median(p_value, na.rm = TRUE),
    combined_snp_p      = mean(null_jaccards <= median(jaccard_sim, na.rm = TRUE)),
    # Gene-level
    median_gene_jaccard = median(gene_jaccard, na.rm = TRUE),
    mean_gene_jaccard   = mean(gene_jaccard, na.rm = TRUE),
    sd_gene_jaccard     = sd(gene_jaccard, na.rm = TRUE),
    median_gene_p       = median(gene_p_value, na.rm = TRUE),
    combined_gene_p     = mean(null_gene_jaccards <= median(gene_jaccard, na.rm = TRUE)),
    # LD ratio
    median_ld_ratio     = median(ld_ratio, na.rm = TRUE)
  ), by = .(cluster_2)]  # group by AFR cluster

  # Add pathway-level aggregates if available
  if (has_pathway && "pathway_jaccard" %in% colnames(cross_results)) {
    pw_agg <- cross_results[, .(
      median_pathway_jaccard = median(pathway_jaccard, na.rm = TRUE),
      mean_pathway_jaccard   = mean(pathway_jaccard, na.rm = TRUE),
      sd_pathway_jaccard     = sd(pathway_jaccard, na.rm = TRUE),
      median_pathway_p       = median(pathway_p_value, na.rm = TRUE),
      combined_pathway_p     = mean(null_pathway_jaccards <= median(pathway_jaccard, na.rm = TRUE)),
      median_pathway_ratio   = median(pathway_ratio, na.rm = TRUE),
      # Ranked pathway metrics
      median_pathway_spearman      = median(pathway_spearman, na.rm = TRUE),
      mean_pathway_spearman        = mean(pathway_spearman, na.rm = TRUE),
      sd_pathway_spearman          = sd(pathway_spearman, na.rm = TRUE),
      median_pathway_top10_overlap = median(pathway_top10_overlap, na.rm = TRUE),
      median_pathway_top25_overlap = median(pathway_top25_overlap, na.rm = TRUE),
      median_pathway_top50_overlap = median(pathway_top50_overlap, na.rm = TRUE)
    ), by = .(cluster_2)]
    agg <- merge(agg, pw_agg, by = "cluster_2", all.x = TRUE)
  }

  setnames(agg, "cluster_2", "afr_cluster")
  print(agg)

  fwrite(agg, file.path(jaccard_dir, "per_cluster_summary.tsv"), sep = "\t")
}


# =====================================================================
# SECTION 6: Save pairwise EUR-EUR matching details
# =====================================================================

cat("\n--- Saving pairwise EUR-EUR results ---\n")

eur_eur_all <- rbindlist(
  lapply(names(pairwise_results), function(pair_label) {
    r <- pairwise_results[[pair_label]]
    if (nrow(r$matching) > 0) {
      r$matching$pair <- pair_label
      as.data.table(r$matching)
    }
  }),
  fill = TRUE
)

fwrite(eur_eur_all, file.path(jaccard_dir, "eur_eur_pairwise.tsv"), sep = "\t")


# =====================================================================
# SECTION 7: Save K-value summary
# =====================================================================

k_summary <- data.table(
  run = c(sub_names, "AFR"),
  K   = c(eur_K, afr_K),
  type = c(rep("EUR_subsample", length(eur_K)), "AFR")
)
fwrite(k_summary, file.path(jaccard_dir, "k_values.tsv"), sep = "\t")


cat("\n=== Jaccard analysis complete ===\n")
cat(sprintf("Results in: %s\n", jaccard_dir))
cat("Files:\n")
cat("  snp_to_gene_map.tsv         - Precomputed SNP-to-gene mapping\n")
cat("  null_distribution.tsv       - EUR-EUR null (SNP + gene + pathway Jaccard)\n")
cat("  cross_ancestry_results.tsv  - All EUR-AFR matched cluster Jaccards + p-values\n")
cat("  cross_ancestry_jaccards.tsv - Cross-ancestry Jaccard values (SNP + gene + pathway)\n")
cat("  per_cluster_summary.tsv     - Aggregate per AFR cluster (SNP + gene + pathway)\n")
cat("  eur_eur_pairwise.tsv        - All EUR-EUR pairwise matching details\n")
cat("  k_values.tsv                - K values per run\n")
