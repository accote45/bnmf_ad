# figure_utils.R
# Shared utilities for multi-panel figure generation:
#   theme, data loading, cosine similarity, cluster matching,
#   GTF gene parsing, SNP-to-gene mapping.

# data.table is loaded before tidyverse so that tidyverse masks the overlapping
# verbs (first/last/between/transpose) for the tidyverse-style consumers here;
# the data.table-only helpers below call fread/`:=`/.SD/fifelse, which are not masked.
library(data.table)
library(tidyverse)


# --- Custom ggplot theme ---

theme_big_text <- function(base_size = 18, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title    = element_text(size = base_size * 1.2, face = "bold"),
      plot.subtitle = element_text(size = base_size, face = "italic"),
      axis.title    = element_text(size = base_size),
      axis.text     = element_text(size = base_size * 0.9),
      legend.title  = element_text(size = base_size),
      legend.text   = element_text(size = base_size * 0.9),
      strip.text    = element_text(size = base_size),
      panel.grid.major = element_line(linewidth = 0.5),
      panel.grid.minor = element_line(linewidth = 0.25)
    )
}


# --- Data loading ---

#' Load all H matrices from results directory
#'
#' @param results_dir Base results directory (e.g., "results/a1_analysis")
#' @param ancestries Character vector of ancestry codes (e.g., c("EUR", "AFR"))
#' @return Named list of tibbles, keyed by ancestry
load_h_matrices <- function(results_dir, ancestries) {
  h_list <- set_names(
    map(ancestries, function(anc) {
      path <- file.path(results_dir, anc, sprintf("H_matrix_%s.tsv", anc))
      if (!file.exists(path)) {
        warning(sprintf("H_matrix_%s.tsv not found for %s at %s", anc, anc, path))
        return(NULL)
      }
      read_tsv(path, show_col_types = FALSE)
    }),
    ancestries
  )
  h_list[!map_lgl(h_list, is.null)]
}


#' Extract base trait names from H-matrix column names
#'
#' @param h_colnames Column names of H-matrix (includes "Cluster" + trait_pos/trait_neg)
#' @return Character vector of unique base trait names (e.g., "T2D", "LDL")
extract_trait_names <- function(h_colnames) {
  trait_cols <- setdiff(h_colnames, "Cluster")
  unique(str_replace(trait_cols, "_(pos|neg)$", ""))
}


#' Build signed trait loadings (pos - neg) from an H matrix file
#'
#' For each cluster and trait, the signed loading is the positive-arm weight
#' minus the negative-arm weight. data.table-based for speed.
#'
#' @param h_file Path to an H_matrix_<anc>.tsv (Cluster column + trait_pos/trait_neg)
#' @param strip_source When TRUE, drop the trailing GWAS-source tag from trait names
#'   (e.g. "LDL_Willer2013" -> "LDL") and keep the max |loading| per (cluster, trait)
#' @return data.table with columns: cluster, trait, loading, abs_loading
build_trait_loadings <- function(h_file, strip_source = FALSE) {
  h <- fread(h_file)
  clusters <- h$Cluster
  h[, Cluster := NULL]

  cols <- names(h)
  pos_cols <- grep("_pos$", cols, value = TRUE)
  neg_cols <- grep("_neg$", cols, value = TRUE)

  pos_bases <- sub("_pos$", "", pos_cols)
  neg_bases <- sub("_neg$", "", neg_cols)
  shared <- intersect(pos_bases, neg_bases)

  result <- rbindlist(lapply(seq_along(clusters), function(i) {
    rbindlist(lapply(shared, function(trait) {
      pos_val <- as.numeric(h[i, get(paste0(trait, "_pos"))])
      neg_val <- as.numeric(h[i, get(paste0(trait, "_neg"))])
      signed <- pos_val - neg_val

      clean_name <- trait
      if (strip_source) {
        clean_name <- sub("_[A-Za-z]+\\d{4}$", "", trait)
      }

      data.table(cluster = clusters[i], trait = clean_name,
                 loading = signed, abs_loading = abs(signed))
    }))
  }))

  # Same trait from multiple GWAS → keep the one with max abs loading per cluster
  if (strip_source) {
    result <- result[, .SD[which.max(abs_loading)], by = .(cluster, trait)]
  }

  result
}


# --- Gene mapping ---

#' Parse a GTF file to extract protein-coding gene coordinates
#'
#' @param gtf_path Path to an Ensembl GTF file (GRCh37)
#' @return tibble with columns: chr, start, end, gene_name
parse_gtf_genes <- function(gtf_path) {
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
    mutate(gene_name = str_replace(attributes, '.*gene_name "([^"]+)".*', "\\1")) %>%
    filter(chr %in% c(as.character(1:22), "X", "Y")) %>%
    select(chr, start, end, gene_name)
}


#' Map SNPs to genes using gene-body overlap
#'
#' @param snp_df tibble with columns: chr, pos, VAR_ID, plus cluster weight columns
#' @param gene_df tibble from parse_gtf_genes()
#' @return tibble with gene_name joined to SNPs (SNPs not in any gene are dropped)
map_snps_to_genes <- function(snp_df, gene_df) {
  snp_df %>%
    inner_join(gene_df, by = "chr", relationship = "many-to-many") %>%
    filter(pos >= start, pos <= end)
}


#' Build a per-cluster gene scatter data frame from a W matrix
#'
#' For each cluster (W-matrix column) and each gene (via SNP-to-gene overlap):
#'   max_loading  = max W[snp, k] over the gene's SNPs
#'   specificity  = mean(W[snp, k] / sum_k W[snp, k]) over SNPs with W[snp, k] > 0
#' The top-N genes by max_loading per cluster are flagged for highlighting/labelling.
#' Facets are ordered by cluster-column order (K1, K2, ...).
#'
#' @param results_dir Base results directory (contains <anc>/W_matrix_<anc>.tsv)
#' @param anc Ancestry code (e.g. "META")
#' @param gtf_path Path to GTF used for gene coordinates
#' @param cluster_labels Named list (per ancestry) mapping "<anc> K#" -> display label
#' @param top_n Number of genes to flag per cluster (default 5)
#' @return tibble with columns: gene_name, max_loading, specificity, is_top, cluster_label
#'   (cluster_label is an ordered factor)
build_gene_scatter_df <- function(results_dir, anc, gtf_path, cluster_labels,
                                  top_n = 5) {
  gene_df <- parse_gtf_genes(gtf_path)
  w_path  <- file.path(results_dir, anc, sprintf("W_matrix_%s.tsv", anc))
  if (!file.exists(w_path)) stop(sprintf("W matrix not found: %s", w_path))

  w_df <- read_tsv(w_path, show_col_types = FALSE) %>%
    mutate(chr = str_replace(VAR_ID, "^([^_]+)_.*", "\\1"),
           pos = as.numeric(str_replace(VAR_ID, "^[^_]+_([^_]+)_.*", "\\1")))
  cluster_cols <- setdiff(colnames(w_df), c("VAR_ID", "chr", "pos"))

  w_df <- w_df %>%
    mutate(.row_total = rowSums(across(all_of(cluster_cols)), na.rm = TRUE))

  mapped <- map_snps_to_genes(w_df, gene_df)
  cat(sprintf("  %s: %d SNPs mapped to genes (of %d total)\n",
              anc, n_distinct(mapped$VAR_ID), nrow(w_df)))

  # One long data frame of (gene, max_loading, specificity) per cluster
  scatter_df <- map_dfr(cluster_cols, function(kcol) {
    gene_max <- mapped %>%
      group_by(gene_name) %>%
      summarise(max_loading = max(.data[[kcol]], na.rm = TRUE), .groups = "drop") %>%
      filter(max_loading > 0)

    gene_spec <- mapped %>%
      filter(.data[[kcol]] > 0, .row_total > 0) %>%
      mutate(.spec = .data[[kcol]] / .row_total) %>%
      group_by(gene_name) %>%
      summarise(specificity = mean(.spec, na.rm = TRUE), .groups = "drop")

    gene_combined <- gene_max %>%
      left_join(gene_spec, by = "gene_name") %>%
      mutate(specificity = replace_na(specificity, 0))

    # Flag the top-N genes by max loading for highlighting + labelling
    top_genes <- gene_combined %>%
      slice_max(max_loading, n = top_n, with_ties = FALSE) %>%
      pull(gene_name)

    facet_key <- paste(anc, kcol)
    display_label <- if (!is.null(cluster_labels[[anc]])) {
      cluster_labels[[anc]][[facet_key]] %||% facet_key
    } else {
      facet_key
    }

    gene_combined %>%
      mutate(is_top = gene_name %in% top_genes,
             cluster_label = display_label)
  })

  # Order facets by the cluster-column order (K1, K2, ...)
  facet_levels <- map_chr(cluster_cols, function(kcol) {
    facet_key <- paste(anc, kcol)
    if (!is.null(cluster_labels[[anc]])) cluster_labels[[anc]][[facet_key]] %||% facet_key
    else facet_key
  })
  scatter_df %>%
    mutate(cluster_label = factor(cluster_label, levels = facet_levels))
}


# --- Cosine similarity and cluster matching ---

#' Compute cosine similarity between two numeric vectors
#'
#' @param a Numeric vector
#' @param b Numeric vector (same length as a)
#' @return Scalar cosine similarity
cosine_similarity <- function(a, b) {
  denom <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  if (denom < 1e-15) return(0)
  sum(a * b) / denom
}


#' Match clusters between two ancestries using cosine similarity
#'
#' Greedy best-match: iterates over the ancestry with fewer clusters,
#' assigning each to its best unmatched partner in the other ancestry.
#'
#' @param h1 tibble H-matrix for ancestry 1 (with Cluster column)
#' @param h2 tibble H-matrix for ancestry 2 (with Cluster column)
#' @param sim_threshold Minimum cosine similarity for a valid match (default 0.3)
#' @return tibble with columns: cluster_1, cluster_2, cosine_sim
#'   Unmatched clusters have NA in the partner column.
match_clusters <- function(h1, h2, sim_threshold = 0.3) {
  clusters_1 <- h1$Cluster
  clusters_2 <- h2$Cluster
  mat1 <- h1 %>% select(-Cluster) %>% as.matrix()
  mat2 <- h2 %>% select(-Cluster) %>% as.matrix()

  common_cols <- intersect(colnames(mat1), colnames(mat2))
  mat1 <- mat1[, common_cols, drop = FALSE]
  mat2 <- mat2[, common_cols, drop = FALSE]

  sim_mat <- matrix(0, nrow = length(clusters_1), ncol = length(clusters_2))
  for (i in seq_along(clusters_1)) {
    for (j in seq_along(clusters_2)) {
      sim_mat[i, j] <- cosine_similarity(mat1[i, ], mat2[j, ])
    }
  }
  rownames(sim_mat) <- clusters_1
  colnames(sim_mat) <- clusters_2

  if (length(clusters_1) <= length(clusters_2)) {
    src_clusters <- clusters_1
    tgt_clusters <- clusters_2
    sm <- sim_mat
    flipped <- FALSE
  } else {
    src_clusters <- clusters_2
    tgt_clusters <- clusters_1
    sm <- t(sim_mat)
    flipped <- TRUE
  }

  matches <- list()
  used_tgt <- character()

  for (i in seq_along(src_clusters)) {
    available <- setdiff(tgt_clusters, used_tgt)
    if (length(available) == 0) {
      matches[[length(matches) + 1]] <- tibble(
        cluster_1 = if (!flipped) src_clusters[i] else NA_character_,
        cluster_2 = if (flipped) src_clusters[i] else NA_character_,
        cosine_sim = NA_real_
      )
      next
    }
    avail_idx <- match(available, tgt_clusters)
    sims <- sm[i, avail_idx]
    best <- which.max(sims)
    best_sim <- sims[best]
    best_tgt <- available[best]

    if (best_sim >= sim_threshold) {
      used_tgt <- c(used_tgt, best_tgt)
      matches[[length(matches) + 1]] <- tibble(
        cluster_1 = if (!flipped) src_clusters[i] else best_tgt,
        cluster_2 = if (flipped) src_clusters[i] else best_tgt,
        cosine_sim = best_sim
      )
    } else {
      matches[[length(matches) + 1]] <- tibble(
        cluster_1 = if (!flipped) src_clusters[i] else NA_character_,
        cluster_2 = if (flipped) src_clusters[i] else NA_character_,
        cosine_sim = NA_real_
      )
    }
  }

  unmatched_tgt <- setdiff(tgt_clusters, used_tgt)
  if (length(unmatched_tgt) > 0) {
    matches[[length(matches) + 1]] <- tibble(
      cluster_1 = if (!flipped) NA_character_ else unmatched_tgt,
      cluster_2 = if (flipped) NA_character_ else unmatched_tgt,
      cosine_sim = NA_real_
    )
  }

  bind_rows(matches)
}


#' Match clusters across all ancestry pairs
#'
#' @param h_list Named list of H-matrix tibbles
#' @param sim_threshold Minimum cosine similarity for valid match
#' @param reference Optional: name of reference ancestry. When set, only pairs
#'   that include the reference ancestry are generated. NULL (default) generates
#'   all pairwise combinations.
#' @return tibble with columns: anc1, anc2, cluster_1, cluster_2, cosine_sim
match_all_clusters <- function(h_list, sim_threshold = 0.3, reference = NULL) {
  anc_names <- names(h_list)
  if (length(anc_names) < 2) {
    warning("Need at least 2 ancestries for cluster matching")
    return(tibble(anc1 = character(), anc2 = character(),
                  cluster_1 = character(), cluster_2 = character(),
                  cosine_sim = numeric()))
  }

  if (!is.null(reference)) {
    if (!reference %in% anc_names) {
      stop(sprintf("Reference ancestry '%s' not found in h_list. Available: %s",
                   reference, paste(anc_names, collapse = ", ")))
    }
    others <- setdiff(anc_names, reference)
    pairs <- map(others, ~ c(reference, .x))
  } else {
    pairs <- combn(anc_names, 2, simplify = FALSE)
  }

  map_dfr(pairs, function(pair) {
    anc1_name <- pair[1]
    anc2_name <- pair[2]
    m <- match_clusters(h_list[[anc1_name]], h_list[[anc2_name]],
                        sim_threshold = sim_threshold)
    m$anc1 <- anc1_name
    m$anc2 <- anc2_name
    m
  }) %>%
    select(anc1, anc2, cluster_1, cluster_2, cosine_sim)
}


#' Compute net weight (pos - neg) for one cluster across all traits
#'
#' @param H_dt data.table/tibble H-matrix with Cluster column
#' @param cluster_label Cluster name (e.g., "K1")
#' @param trait_vec Character vector of base trait names
#' @return Named numeric vector of net weights
get_net_vec <- function(H_dt, cluster_label, trait_vec) {
  row <- H_dt[H_dt$Cluster == cluster_label, ]
  h_cols <- colnames(row)
  sapply(trait_vec, function(tr) {
    pc <- paste0(tr, "_pos")
    nc <- paste0(tr, "_neg")
    pv <- if (pc %in% h_cols) as.numeric(row[[pc]]) else 0
    nv <- if (nc %in% h_cols) as.numeric(row[[nc]]) else 0
    pv - nv
  })
}


#' Match clusters between two ancestries using trait-correlation + Hungarian algorithm
#'
#' Uses Pearson correlation of net trait-weight profiles and Hungarian assignment
#' to find optimal matching. Falls back to greedy matching if clue is unavailable.
#'
#' @param h1 H-matrix tibble/data.table for ancestry 1
#' @param h2 H-matrix tibble/data.table for ancestry 2
#' @return List with $assignment (tibble), $cor_mat (matrix), $method (character)
match_clusters_by_trait_correlation <- function(h1, h2) {
  traits1 <- extract_trait_names(colnames(h1))
  traits2 <- extract_trait_names(colnames(h2))
  common_traits <- intersect(traits1, traits2)

  clusters_1 <- h1$Cluster
  clusters_2 <- h2$Cluster

  # Precompute net vectors
  vecs1 <- setNames(
    lapply(clusters_1, function(k) get_net_vec(h1, k, common_traits)),
    clusters_1
  )
  vecs2 <- setNames(
    lapply(clusters_2, function(k) get_net_vec(h2, k, common_traits)),
    clusters_2
  )

  # Build correlation matrix
  cor_mat <- matrix(0, nrow = length(clusters_1), ncol = length(clusters_2),
                    dimnames = list(clusters_1, clusters_2))
  for (i in seq_along(clusters_1)) {
    for (j in seq_along(clusters_2)) {
      cor_mat[i, j] <- cor(vecs1[[i]], vecs2[[j]], method = "pearson")
    }
  }

  # Hungarian assignment (requires non-negative costs, shift by +1)
  n1 <- length(clusters_1)
  n2 <- length(clusters_2)
  shifted <- cor_mat + 1  # map [-1,1] -> [0,2]

  use_hungarian <- requireNamespace("clue", quietly = TRUE)
  if (use_hungarian) {
    if (n1 != n2) {
      max_k <- max(n1, n2)
      padded <- matrix(0, nrow = max_k, ncol = max_k)
      padded[seq_len(n1), seq_len(n2)] <- shifted
      assignment_idx <- clue::solve_LSAP(padded, maximum = TRUE)[seq_len(n1)]
    } else {
      assignment_idx <- clue::solve_LSAP(shifted, maximum = TRUE)
    }
    valid_pairs <- assignment_idx <= n2
    assignment <- tibble(
      cluster_1 = clusters_1[valid_pairs],
      cluster_2 = clusters_2[assignment_idx[valid_pairs]],
      pearson_r = mapply(function(i, j) cor_mat[i, j],
                         which(valid_pairs), assignment_idx[valid_pairs])
    )
    method <- "hungarian"
  } else {
    # Greedy fallback
    assignment <- list()
    used <- character()
    for (i in seq_along(clusters_1)) {
      avail <- setdiff(clusters_2, used)
      if (length(avail) == 0) break
      sims <- cor_mat[i, avail, drop = TRUE]
      best <- avail[which.max(sims)]
      used <- c(used, best)
      assignment[[length(assignment) + 1]] <- tibble(
        cluster_1 = clusters_1[i],
        cluster_2 = best,
        pearson_r = cor_mat[i, best]
      )
    }
    assignment <- bind_rows(assignment)
    method <- "greedy"
  }

  list(assignment = assignment, cor_mat = cor_mat, method = method,
       net_vecs_1 = vecs1, net_vecs_2 = vecs2, common_traits = common_traits)
}
