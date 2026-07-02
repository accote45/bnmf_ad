# jaccard_utils.R
# Shared utility functions for Jaccard analysis across C1.1 pipeline scripts.
# Sourced by: 04_compute_jaccard.R, 05_visualize_results.R, 06_simulate_ld_null.R
#
# Functions:
#   - jaccard_similarity()          SNP-level Jaccard
#   - get_cluster_snp_set()         Extract top-N SNPs from W matrix
#   - parse_gtf_genes_dt()          Parse GTF for protein-coding genes
#   - build_snp_to_gene_map()       Map VAR_IDs to genes via foverlaps
#   - gene_level_jaccard()          Gene-level Jaccard
#   - parse_gmt()                   Parse GMT pathway file
#   - get_enriched_pathways()       Binary enrichment (FDR threshold)
#   - get_pathway_scores()          Continuous enrichment scores (-log10 FDR)
#   - pathway_level_jaccard()       Pathway-level Jaccard (binary)
#   - ranked_pathway_similarity()   Ranked pathway metrics (Spearman + top-k)
#   - hungarian_match_clusters()    Optimal cluster matching via Hungarian algorithm

library(data.table)


# =====================================================================
# SNP-level functions
# =====================================================================

#' Compute Jaccard similarity between two sets
#'
#' @param set_a Character vector
#' @param set_b Character vector
#' @return Scalar Jaccard similarity in [0, 1]
jaccard_similarity <- function(set_a, set_b) {
  if (length(set_a) == 0 && length(set_b) == 0) return(0)
  n_intersect <- length(intersect(set_a, set_b))
  n_union     <- length(union(set_a, set_b))
  if (n_union == 0) return(0)
  n_intersect / n_union
}


#' Extract top-N SNP set for a cluster from W matrix
#'
#' Filters out NMF residual noise using a relative threshold: only SNPs with
#' weight >= rel_threshold * max(weight) for this cluster are considered.
#' This prevents near-zero NMF artefacts (e.g., 1e-65) from making every
#' cluster return the same full SNP set when the W matrix is small.
#'
#' @param W_matrix data.table or matrix with VAR_ID column/rownames and cluster columns
#' @param cluster_col Name of the cluster column (e.g., "K1")
#' @param top_n Number of top-weighted variants to return
#' @param rel_threshold Minimum weight as fraction of cluster max (default 0.01)
#' @return Character vector of VAR_IDs
get_cluster_snp_set <- function(W_matrix, cluster_col, top_n = 100,
                                rel_threshold = 0.01) {
  if (is.data.table(W_matrix) || is.data.frame(W_matrix)) {
    var_ids <- W_matrix$VAR_ID
    weights <- W_matrix[[cluster_col]]
  } else {
    # Matrix with rownames
    var_ids <- rownames(W_matrix)
    weights <- W_matrix[, cluster_col]
  }

  names(weights) <- var_ids
  weights <- weights[weights > 0]

  if (length(weights) == 0) return(character(0))

  # Relative threshold to filter NMF residual noise (disabled for now —
  # the Spearman-based metrics already handle this by using continuous weights
  # rather than binary set membership; re-enable if Jaccard matching degenerates)
  # max_w <- max(weights)
  # if (max_w > 0) {
  #   weights <- weights[weights >= rel_threshold * max_w]
  # }
  # if (length(weights) == 0) return(character(0))

  sorted <- sort(weights, decreasing = TRUE)
  names(head(sorted, min(top_n, length(sorted))))
}


# =====================================================================
# Gene-level functions
# =====================================================================

#' Parse GTF for protein-coding gene coordinates (data.table version)
#'
#' @param gtf_path Path to Ensembl GTF (GRCh37)
#' @return data.table with chr (integer), start, end, gene_name, gene_id
parse_gtf_genes_dt <- function(gtf_path) {
  gtf <- fread(gtf_path, sep = "\t", header = FALSE,
               select = c(1, 3, 4, 5, 9),
               col.names = c("chr", "feature", "start", "end", "attributes"))
  gtf <- gtf[feature == "gene" & grepl('gene_biotype "protein_coding"', attributes)]
  gtf[, gene_name := sub('.*gene_name "([^"]+)".*', "\\1", attributes)]
  gtf[, gene_id := sub("\\..*", "", sub('.*gene_id "([^"]+)".*', "\\1", attributes))]
  gtf <- gtf[chr %in% as.character(1:22)]
  gtf[, chr := as.integer(chr)]
  gtf[, .(chr, start, end, gene_name, gene_id)]
}


#' Build SNP-to-gene mapping table from a set of VAR_IDs
#'
#' Uses foverlaps for efficient interval overlap with flanking window.
#'
#' @param var_ids Character vector of VAR_IDs (format: CHR_POS_A1_A2)
#' @param gene_dt data.table from parse_gtf_genes_dt()
#' @param flank_bp Flanking window in bp (added to both sides of gene body)
#' @return data.table with columns: VAR_ID, gene_name, gene_id
build_snp_to_gene_map <- function(var_ids, gene_dt, flank_bp = 50000) {
  snp_dt <- data.table(VAR_ID = unique(var_ids))
  snp_dt[, c("chr", "pos") := tstrsplit(VAR_ID, "_", keep = 1:2)]
  snp_dt[, chr := as.integer(chr)]
  snp_dt[, pos := as.integer(pos)]
  snp_dt[, pos_end := pos]  # point interval for foverlaps

  gene_expanded <- copy(gene_dt)
  gene_expanded[, start_f := pmax(0L, as.integer(start - flank_bp))]
  gene_expanded[, end_f   := as.integer(end + flank_bp)]

  # foverlaps requires keyed intervals
  setkey(gene_expanded, chr, start_f, end_f)
  setkey(snp_dt, chr, pos, pos_end)

  mapped <- foverlaps(snp_dt, gene_expanded, type = "within", nomatch = NULL)
  unique(mapped[, .(VAR_ID, gene_name, gene_id)])
}


#' Compute gene-level Jaccard between two SNP sets
#'
#' @param snp_set_1 Character vector of VAR_IDs
#' @param snp_set_2 Character vector of VAR_IDs
#' @param snp_to_gene_map data.table with VAR_ID, gene_name columns
#' @return Scalar Jaccard similarity on gene sets [0, 1]
gene_level_jaccard <- function(snp_set_1, snp_set_2, snp_to_gene_map) {
  genes_1 <- unique(snp_to_gene_map[VAR_ID %in% snp_set_1, gene_name])
  genes_2 <- unique(snp_to_gene_map[VAR_ID %in% snp_set_2, gene_name])
  jaccard_similarity(genes_1, genes_2)
}


#' Compute SNP-level Spearman correlation between matched clusters
#'
#' Extracts full weight vectors from W matrices for two matched clusters,
#' aligns on the union of all VAR_IDs (0 for absent), and computes Spearman rho.
#'
#' @param W1 W matrix data.table with VAR_ID + cluster columns
#' @param W2 W matrix data.table with VAR_ID + cluster columns
#' @param cluster_1 Name of cluster column in W1
#' @param cluster_2 Name of cluster column in W2
#' @return Scalar Spearman rho, or NA if fewer than 3 non-zero SNPs
snp_level_spearman <- function(W1, W2, cluster_1, cluster_2) {
  v1 <- W1[[cluster_1]]; names(v1) <- W1$VAR_ID
  v2 <- W2[[cluster_2]]; names(v2) <- W2$VAR_ID
  all_snps <- union(names(v1), names(v2))
  w1 <- setNames(numeric(length(all_snps)), all_snps)
  w2 <- w1
  w1[names(v1)] <- v1
  w2[names(v2)] <- v2
  # Keep only SNPs with non-zero weight in at least one cluster
  keep <- w1 > 0 | w2 > 0
  if (sum(keep) < 3) return(NA_real_)
  suppressWarnings(cor(w1[keep], w2[keep], method = "spearman"))
}


#' Compute gene-level Spearman correlation between matched clusters
#'
#' Aggregates SNP weights to genes (max weight per gene), aligns on the
#' union of genes, and computes Spearman rho.
#'
#' @param W1 W matrix data.table with VAR_ID + cluster columns
#' @param W2 W matrix data.table with VAR_ID + cluster columns
#' @param cluster_1 Name of cluster column in W1
#' @param cluster_2 Name of cluster column in W2
#' @param snp_to_gene_map data.table with VAR_ID and gene_name columns
#' @return Scalar Spearman rho, or NA if fewer than 3 genes
gene_level_spearman <- function(W1, W2, cluster_1, cluster_2, snp_to_gene_map) {
  # Build SNP-to-weight for each cluster
  dt1 <- data.table(VAR_ID = W1$VAR_ID, w = W1[[cluster_1]])
  dt2 <- data.table(VAR_ID = W2$VAR_ID, w = W2[[cluster_2]])
  # Map to genes and take max weight per gene
  g1 <- merge(dt1, snp_to_gene_map[, .(VAR_ID, gene_name)], by = "VAR_ID", allow.cartesian = TRUE)
  g1 <- g1[, .(w = max(w)), by = gene_name]
  g2 <- merge(dt2, snp_to_gene_map[, .(VAR_ID, gene_name)], by = "VAR_ID", allow.cartesian = TRUE)
  g2 <- g2[, .(w = max(w)), by = gene_name]
  # Align on union of genes
  all_genes <- union(g1$gene_name, g2$gene_name)
  w1 <- setNames(numeric(length(all_genes)), all_genes)
  w2 <- w1
  w1[g1$gene_name] <- g1$w
  w2[g2$gene_name] <- g2$w
  keep <- w1 > 0 | w2 > 0
  if (sum(keep) < 3) return(NA_real_)
  suppressWarnings(cor(w1[keep], w2[keep], method = "spearman"))
}


# =====================================================================
# Pathway-level functions
# =====================================================================

#' Parse GMT pathway file into long-format data.table
#'
#' @param gmt_path Path to GMT file (tab-separated: pathway, description, gene1, gene2, ...)
#' @return data.table with columns: pathway (character), gene_id (character, version-stripped ENSG)
parse_gmt <- function(gmt_path) {
  lines <- readLines(gmt_path)
  result <- rbindlist(lapply(lines, function(line) {
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) < 3) return(NULL)
    pathway_name <- fields[1]
    # fields[2] is description (PLACEHOLDER), skip it
    gene_ids <- fields[3:length(fields)]
    gene_ids <- gene_ids[nzchar(gene_ids)]
    # Strip version suffix (ENSG00000186092.4 -> ENSG00000186092)
    gene_ids <- sub("\\..*", "", gene_ids)
    data.table(pathway = pathway_name, gene_id = gene_ids)
  }))
  result
}


#' Get enriched pathways for a gene set using Fisher's exact test (binary)
#'
#' @param cluster_gene_ids Character vector of Ensembl gene IDs in the cluster
#' @param background_gene_ids Character vector of all Ensembl gene IDs in the background
#' @param pathway_map data.table with columns: pathway, gene_id (from parse_gmt)
#' @param fdr_threshold FDR cutoff for enrichment (default 0.05)
#' @return Character vector of significantly enriched pathway names
get_enriched_pathways <- function(cluster_gene_ids, background_gene_ids,
                                  pathway_map, fdr_threshold = 0.05) {
  if (length(cluster_gene_ids) == 0) return(character(0))

  # Restrict to genes that exist in the background
  cluster_in_bg <- intersect(cluster_gene_ids, background_gene_ids)
  if (length(cluster_in_bg) == 0) return(character(0))

  n_bg <- length(background_gene_ids)
  n_cluster <- length(cluster_in_bg)

  # Get unique pathways
  pathways <- unique(pathway_map$pathway)

  # Pre-compute pathway gene sets (only genes in background)
  pathway_genes <- pathway_map[gene_id %in% background_gene_ids,
                                .(genes = list(unique(gene_id))), by = pathway]

  # Fisher's exact test for each pathway
  pvals <- sapply(seq_len(nrow(pathway_genes)), function(i) {
    pw_genes <- pathway_genes$genes[[i]]
    n_pw <- length(pw_genes)
    if (n_pw == 0) return(1)

    # 2x2 contingency table
    a <- length(intersect(cluster_in_bg, pw_genes))  # in cluster AND pathway
    if (a == 0) return(1)  # no overlap, skip expensive test

    b <- n_cluster - a                                 # in cluster, NOT pathway
    c <- n_pw - a                                      # in pathway, NOT cluster
    d <- n_bg - a - b - c                              # neither

    # One-sided Fisher's exact (over-representation)
    fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value
  })

  names(pvals) <- pathway_genes$pathway

  # BH correction
  fdr <- p.adjust(pvals, method = "BH")

  names(fdr)[fdr < fdr_threshold]
}


#' Get continuous pathway enrichment scores for a gene set
#'
#' Returns -log10(FDR) for every pathway tested, enabling ranked comparison
#' without dependence on an arbitrary FDR threshold.
#'
#' @param cluster_gene_ids Character vector of Ensembl gene IDs in the cluster
#' @param background_gene_ids Character vector of all Ensembl gene IDs in background
#' @param pathway_map data.table with columns: pathway, gene_id (from parse_gmt)
#' @return Named numeric vector: -log10(FDR) per pathway, names = pathway names
get_pathway_scores <- function(cluster_gene_ids, background_gene_ids, pathway_map) {
  if (length(cluster_gene_ids) == 0) return(numeric(0))

  cluster_in_bg <- intersect(cluster_gene_ids, background_gene_ids)
  if (length(cluster_in_bg) == 0) return(numeric(0))

  n_bg <- length(background_gene_ids)
  n_cluster <- length(cluster_in_bg)

  pathway_genes <- pathway_map[gene_id %in% background_gene_ids,
                                .(genes = list(unique(gene_id))), by = pathway]

  pvals <- sapply(seq_len(nrow(pathway_genes)), function(i) {
    pw_genes <- pathway_genes$genes[[i]]
    n_pw <- length(pw_genes)
    if (n_pw == 0) return(1)

    a <- length(intersect(cluster_in_bg, pw_genes))
    if (a == 0) return(1)

    b <- n_cluster - a
    c <- n_pw - a
    d <- n_bg - a - b - c

    fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value
  })

  names(pvals) <- pathway_genes$pathway

  fdr <- p.adjust(pvals, method = "BH")
  scores <- -log10(pmax(fdr, .Machine$double.xmin))  # avoid -log10(0)
  names(scores) <- pathway_genes$pathway
  scores
}


#' Compute pathway-level Jaccard between two SNP sets using enrichment (binary)
#'
#' @param snp_set_1 Character vector of VAR_IDs
#' @param snp_set_2 Character vector of VAR_IDs
#' @param snp_to_gene_map data.table with VAR_ID, gene_name, gene_id columns
#' @param pathway_map data.table with pathway, gene_id columns (from parse_gmt)
#' @param background_gene_ids Character vector of all Ensembl gene IDs in background
#' @param fdr_threshold FDR cutoff for pathway enrichment
#' @return Scalar Jaccard similarity on enriched pathway sets [0, 1]
pathway_level_jaccard <- function(snp_set_1, snp_set_2, snp_to_gene_map,
                                   pathway_map, background_gene_ids,
                                   fdr_threshold = 0.05) {
  # Map SNPs to Ensembl gene IDs
  genes_1 <- unique(snp_to_gene_map[VAR_ID %in% snp_set_1, gene_id])
  genes_2 <- unique(snp_to_gene_map[VAR_ID %in% snp_set_2, gene_id])

  # Get enriched pathways for each set
  pw_1 <- get_enriched_pathways(genes_1, background_gene_ids, pathway_map, fdr_threshold)
  pw_2 <- get_enriched_pathways(genes_2, background_gene_ids, pathway_map, fdr_threshold)

  jaccard_similarity(pw_1, pw_2)
}


#' Compute ranked pathway similarity between two SNP sets
#'
#' Uses continuous -log10(FDR) scores instead of binary thresholds.
#' Returns Spearman correlation across all shared pathways and
#' top-k overlap for k = 10, 25, 50.
#'
#' @param snp_set_1 Character vector of VAR_IDs
#' @param snp_set_2 Character vector of VAR_IDs
#' @param snp_to_gene_map data.table with VAR_ID, gene_name, gene_id columns
#' @param pathway_map data.table with pathway, gene_id columns
#' @param background_gene_ids Character vector of all Ensembl gene IDs in background
#' @return Named list: pathway_spearman, pathway_top10_overlap,
#'   pathway_top25_overlap, pathway_top50_overlap
ranked_pathway_similarity <- function(snp_set_1, snp_set_2, snp_to_gene_map,
                                       pathway_map, background_gene_ids) {
  genes_1 <- unique(snp_to_gene_map[VAR_ID %in% snp_set_1, gene_id])
  genes_2 <- unique(snp_to_gene_map[VAR_ID %in% snp_set_2, gene_id])

  scores_1 <- get_pathway_scores(genes_1, background_gene_ids, pathway_map)
  scores_2 <- get_pathway_scores(genes_2, background_gene_ids, pathway_map)

  # Need at least some pathways scored in both sets
  if (length(scores_1) == 0 || length(scores_2) == 0) {
    return(list(
      pathway_spearman       = NA_real_,
      pathway_top10_overlap  = NA_real_,
      pathway_top25_overlap  = NA_real_,
      pathway_top50_overlap  = NA_real_
    ))
  }

  # Shared pathways (both scored)
  shared_pw <- intersect(names(scores_1), names(scores_2))

  # Spearman correlation on shared pathways
  spearman <- if (length(shared_pw) >= 3) {
    suppressWarnings(cor(scores_1[shared_pw], scores_2[shared_pw], method = "spearman"))
  } else {
    NA_real_
  }

  # Top-k overlap: fraction of top-k pathways in set 1 that appear in top-k of set 2
  top_k_overlap <- function(s1, s2, k) {
    if (length(s1) < k || length(s2) < k) {
      k <- min(length(s1), length(s2))
    }
    if (k == 0) return(NA_real_)
    top1 <- names(sort(s1, decreasing = TRUE))[1:k]
    top2 <- names(sort(s2, decreasing = TRUE))[1:k]
    length(intersect(top1, top2)) / k
  }

  list(
    pathway_spearman       = spearman,
    pathway_top10_overlap  = top_k_overlap(scores_1, scores_2, 10),
    pathway_top25_overlap  = top_k_overlap(scores_1, scores_2, 25),
    pathway_top50_overlap  = top_k_overlap(scores_1, scores_2, 50)
  )
}


# =====================================================================
# Cluster matching
# =====================================================================

#' Match clusters between two bNMF runs using Hungarian algorithm
#'
#' @param W1 W matrix from run 1 (data.table with VAR_ID + cluster columns)
#' @param W2 W matrix from run 2
#' @param top_n Top N variants per cluster for Jaccard
#' @return List: matching (data.frame), jaccard_matrix, cluster_snp_sets
hungarian_match_clusters <- function(W1, W2, top_n = 100) {
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("Package 'clue' required. Install with: install.packages('clue')")
  }

  # Identify cluster columns (everything except VAR_ID)
  if (is.data.table(W1) || is.data.frame(W1)) {
    clusters_1 <- setdiff(colnames(W1), "VAR_ID")
    clusters_2 <- setdiff(colnames(W2), "VAR_ID")
  } else {
    clusters_1 <- colnames(W1)
    clusters_2 <- colnames(W2)
  }

  k1 <- length(clusters_1)
  k2 <- length(clusters_2)

  if (k1 == 0 || k2 == 0) {
    cat("  WARNING: One or both runs have K=0. Cannot match.\n")
    return(list(
      matching = data.frame(
        cluster_1 = character(0), cluster_2 = character(0),
        jaccard_sim = numeric(0), stringsAsFactors = FALSE
      ),
      jaccard_matrix = matrix(nrow = 0, ncol = 0),
      snp_sets_1 = list(), snp_sets_2 = list()
    ))
  }

  # Extract SNP sets for each cluster
  sets1 <- setNames(
    lapply(clusters_1, function(k) get_cluster_snp_set(W1, k, top_n)),
    clusters_1
  )
  sets2 <- setNames(
    lapply(clusters_2, function(k) get_cluster_snp_set(W2, k, top_n)),
    clusters_2
  )

  # Build Jaccard similarity matrix
  jaccard_mat <- matrix(0, nrow = k1, ncol = k2,
                        dimnames = list(clusters_1, clusters_2))
  for (i in seq_len(k1)) {
    for (j in seq_len(k2)) {
      jaccard_mat[i, j] <- jaccard_similarity(sets1[[i]], sets2[[j]])
    }
  }

  # Cost matrix = 1 - Jaccard (minimize distance)
  # Pad to square if K differs
  max_k <- max(k1, k2)
  cost_mat <- matrix(1, nrow = max_k, ncol = max_k)
  cost_mat[1:k1, 1:k2] <- 1 - jaccard_mat

  # Solve assignment problem
  assignment <- clue::solve_LSAP(cost_mat)

  # Extract matching
  matches <- list()
  for (i in seq_len(k1)) {
    j <- assignment[i]
    if (j <= k2) {
      matches[[length(matches) + 1]] <- data.frame(
        cluster_1   = clusters_1[i],
        cluster_2   = clusters_2[j],
        jaccard_sim = jaccard_mat[i, j],
        stringsAsFactors = FALSE
      )
    }
  }

  matching <- do.call(rbind, matches)
  if (is.null(matching)) {
    matching <- data.frame(
      cluster_1 = character(0), cluster_2 = character(0),
      jaccard_sim = numeric(0), stringsAsFactors = FALSE
    )
  }

  list(
    matching       = matching,
    jaccard_matrix = jaccard_mat,
    snp_sets_1     = sets1,
    snp_sets_2     = sets2
  )
}
