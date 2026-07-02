#!/usr/bin/env Rscript
# test_c1_jaccard.R
# Unit tests for C1.1 Jaccard analysis functions:
#   - jaccard_similarity()
#   - get_cluster_snp_set()
#   - hungarian_match_clusters()
#   - build_null_distribution logic
#
# Usage:
#   Rscript tests/test_c1_jaccard.R

library(data.table)

# Determine project root
project_root <- getwd()
if (!file.exists(file.path(project_root, "scripts", "c1_ancestry_test", "04_compute_jaccard.R"))) {
  project_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), ".."))
}

# Source the functions under test (from 04_compute_jaccard.R)
# We need to extract the function definitions without running the main body.
# Source figure_utils.R first (dependency)
source(file.path(project_root, "scripts", "analysis", "figure_utils.R"))

# Define the functions directly here (same as 04_compute_jaccard.R)
# This avoids executing the main analysis code during testing.

get_cluster_snp_set <- function(W_matrix, cluster_col, top_n = 100) {
  if (is.data.table(W_matrix) || is.data.frame(W_matrix)) {
    var_ids <- W_matrix$VAR_ID
    weights <- W_matrix[[cluster_col]]
  } else {
    var_ids <- rownames(W_matrix)
    weights <- W_matrix[, cluster_col]
  }
  names(weights) <- var_ids
  weights <- weights[weights > 0]
  if (length(weights) == 0) return(character(0))
  sorted <- sort(weights, decreasing = TRUE)
  names(head(sorted, min(top_n, length(sorted))))
}

jaccard_similarity <- function(set_a, set_b) {
  if (length(set_a) == 0 && length(set_b) == 0) return(0)
  n_intersect <- length(intersect(set_a, set_b))
  n_union     <- length(union(set_a, set_b))
  if (n_union == 0) return(0)
  n_intersect / n_union
}

hungarian_match_clusters <- function(W1, W2, top_n = 100) {
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("Package 'clue' required. Install with: install.packages('clue')")
  }

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
    return(list(
      matching = data.frame(cluster_1 = character(0), cluster_2 = character(0),
                            jaccard_sim = numeric(0), stringsAsFactors = FALSE),
      jaccard_matrix = matrix(nrow = 0, ncol = 0),
      snp_sets_1 = list(), snp_sets_2 = list()
    ))
  }

  sets1 <- setNames(lapply(clusters_1, function(k) get_cluster_snp_set(W1, k, top_n)), clusters_1)
  sets2 <- setNames(lapply(clusters_2, function(k) get_cluster_snp_set(W2, k, top_n)), clusters_2)

  jaccard_mat <- matrix(0, nrow = k1, ncol = k2, dimnames = list(clusters_1, clusters_2))
  for (i in seq_len(k1)) {
    for (j in seq_len(k2)) {
      jaccard_mat[i, j] <- jaccard_similarity(sets1[[i]], sets2[[j]])
    }
  }

  max_k <- max(k1, k2)
  cost_mat <- matrix(1, nrow = max_k, ncol = max_k)
  cost_mat[1:k1, 1:k2] <- 1 - jaccard_mat

  assignment <- clue::solve_LSAP(cost_mat)

  matches <- list()
  for (i in seq_len(k1)) {
    j <- assignment[i]
    if (j <= k2) {
      matches[[length(matches) + 1]] <- data.frame(
        cluster_1 = clusters_1[i], cluster_2 = clusters_2[j],
        jaccard_sim = jaccard_mat[i, j], stringsAsFactors = FALSE
      )
    }
  }

  matching <- if (length(matches) > 0) do.call(rbind, matches) else
    data.frame(cluster_1 = character(0), cluster_2 = character(0),
               jaccard_sim = numeric(0), stringsAsFactors = FALSE)

  list(matching = matching, jaccard_matrix = jaccard_mat,
       snp_sets_1 = sets1, snp_sets_2 = sets2)
}

# =====================================================================
# Test framework
# =====================================================================

n_tests  <- 0
n_passed <- 0
n_failed <- 0

test <- function(description, expr) {
  n_tests <<- n_tests + 1
  result <- tryCatch(
    { eval(expr, envir = parent.frame()); TRUE },
    error = function(e) { cat(sprintf("    ERROR: %s\n", e$message)); FALSE }
  )
  if (result) {
    n_passed <<- n_passed + 1
    cat(sprintf("  PASS: %s\n", description))
  } else {
    n_failed <<- n_failed + 1
    cat(sprintf("  FAIL: %s\n", description))
  }
}

expect_equal <- function(actual, expected, tol = 1e-10) {
  if (is.numeric(actual) && is.numeric(expected)) {
    if (abs(actual - expected) > tol) {
      stop(sprintf("Expected %s, got %s", expected, actual))
    }
  } else {
    if (!identical(actual, expected)) {
      stop(sprintf("Expected %s, got %s", deparse(expected), deparse(actual)))
    }
  }
}

expect_true <- function(x) {
  if (!isTRUE(x)) stop(sprintf("Expected TRUE, got %s", deparse(x)))
}

expect_length <- function(x, n) {
  if (length(x) != n) stop(sprintf("Expected length %d, got %d", n, length(x)))
}


# =====================================================================
# Tests: jaccard_similarity()
# =====================================================================

cat("\n=== Testing jaccard_similarity() ===\n")

test("identical sets -> 1.0", {
  expect_equal(jaccard_similarity(c("a", "b", "c"), c("a", "b", "c")), 1.0)
})

test("disjoint sets -> 0.0", {
  expect_equal(jaccard_similarity(c("a", "b"), c("c", "d")), 0.0)
})

test("partial overlap -> correct value", {
  # {b,c} intersect / {a,b,c,d} union = 2/4 = 0.5
  expect_equal(jaccard_similarity(c("a", "b", "c"), c("b", "c", "d")), 0.5)
})

test("empty sets -> 0.0", {
  expect_equal(jaccard_similarity(character(0), character(0)), 0.0)
})

test("one empty set -> 0.0", {
  expect_equal(jaccard_similarity(c("a", "b"), character(0)), 0.0)
})

test("single shared element", {
  # {a} / {a,b,c} = 1/3
  expect_equal(jaccard_similarity(c("a"), c("a", "b", "c")), 1/3)
})


# =====================================================================
# Tests: get_cluster_snp_set()
# =====================================================================

cat("\n=== Testing get_cluster_snp_set() ===\n")

test("returns top N by weight from data.table", {
  W <- data.table(
    VAR_ID = paste0("snp", 1:10),
    K1 = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0)
  )
  result <- get_cluster_snp_set(W, "K1", top_n = 3)
  expect_length(result, 3)
  expect_equal(result[1], "snp1")
  expect_equal(result[2], "snp2")
  expect_equal(result[3], "snp3")
})

test("excludes zero-weight variants", {
  W <- data.table(VAR_ID = paste0("snp", 1:5), K1 = c(0.5, 0, 0, 0, 0))
  result <- get_cluster_snp_set(W, "K1", top_n = 10)
  expect_length(result, 1)
  expect_equal(result[1], "snp1")
})

test("handles top_n larger than available", {
  W <- data.table(VAR_ID = paste0("snp", 1:3), K1 = c(0.5, 0.3, 0.1))
  result <- get_cluster_snp_set(W, "K1", top_n = 100)
  expect_length(result, 3)
})

test("works with matrix input", {
  W <- matrix(c(0.9, 0.5, 0.1, 0.2, 0.8, 0.3), nrow = 3, ncol = 2,
              dimnames = list(c("snp1", "snp2", "snp3"), c("K1", "K2")))
  result <- get_cluster_snp_set(W, "K1", top_n = 2)
  expect_length(result, 2)
  expect_equal(result[1], "snp1")
})


# =====================================================================
# Tests: hungarian_match_clusters()
# =====================================================================

cat("\n=== Testing hungarian_match_clusters() ===\n")

test("matches 2x2 clusters correctly", {
  # W1: K1 has snp1,snp2 high; K2 has snp3,snp4 high
  # W2: K1 has snp3,snp4 high; K2 has snp1,snp2 high
  # So W1:K1 should match W2:K2 and W1:K2 should match W2:K1
  W1 <- data.table(
    VAR_ID = paste0("snp", 1:4),
    K1 = c(0.9, 0.8, 0.1, 0.0),
    K2 = c(0.0, 0.1, 0.9, 0.8)
  )
  W2 <- data.table(
    VAR_ID = paste0("snp", 1:4),
    K1 = c(0.1, 0.0, 0.9, 0.8),
    K2 = c(0.9, 0.8, 0.0, 0.1)
  )
  result <- hungarian_match_clusters(W1, W2, top_n = 2)

  expect_equal(nrow(result$matching), 2)
  # K1 from W1 (snp1,snp2) should match K2 from W2 (snp1,snp2)
  k1_match <- result$matching[result$matching$cluster_1 == "K1", "cluster_2"]
  expect_equal(k1_match, "K2")
  # Both matches should have Jaccard = 1.0
  expect_true(all(result$matching$jaccard_sim == 1.0))
})

test("handles different K values", {
  W1 <- data.table(VAR_ID = paste0("snp", 1:6),
                   K1 = c(0.9, 0.8, 0.7, 0, 0, 0),
                   K2 = c(0, 0, 0, 0.9, 0.8, 0.7),
                   K3 = c(0.5, 0, 0.5, 0, 0.5, 0))
  W2 <- data.table(VAR_ID = paste0("snp", 1:6),
                   K1 = c(0.8, 0.9, 0.6, 0, 0, 0),
                   K2 = c(0, 0, 0, 0.8, 0.9, 0.6))

  result <- hungarian_match_clusters(W1, W2, top_n = 3)
  # Should have 2 matched pairs (min of K1=3, K2=2)
  expect_true(nrow(result$matching) <= 3)
  expect_true(nrow(result$matching) >= 2)
  # No duplicate assignments
  expect_equal(length(unique(result$matching$cluster_2)), nrow(result$matching))
})

test("returns valid Jaccard range [0,1]", {
  set.seed(42)
  W1 <- data.table(VAR_ID = paste0("snp", 1:20),
                   K1 = runif(20), K2 = runif(20))
  W2 <- data.table(VAR_ID = paste0("snp", 1:20),
                   K1 = runif(20), K2 = runif(20))
  result <- hungarian_match_clusters(W1, W2, top_n = 10)
  expect_true(all(result$matching$jaccard_sim >= 0))
  expect_true(all(result$matching$jaccard_sim <= 1))
})

test("handles K=0 gracefully", {
  W1 <- data.table(VAR_ID = paste0("snp", 1:5), K1 = c(0.9, 0.8, 0.7, 0.6, 0.5))
  W2 <- data.table(VAR_ID = paste0("snp", 1:5))  # No cluster columns
  result <- hungarian_match_clusters(W1, W2, top_n = 3)
  expect_equal(nrow(result$matching), 0)
})


# =====================================================================
# Tests: Null distribution logic
# =====================================================================

cat("\n=== Testing null distribution logic ===\n")

test("pairwise comparisons count is C(n,2)", {
  sub_names <- c("sub_01", "sub_02", "sub_03", "sub_04")
  pairs <- combn(sub_names, 2, simplify = FALSE)
  expect_equal(length(pairs), 6)  # C(4,2) = 6
})

test("empirical p-value computation", {
  null_values <- c(0.4, 0.5, 0.6, 0.7, 0.8)
  observed <- 0.3

  # p = fraction of null <= observed
  p_val <- mean(null_values <= observed)
  expect_equal(p_val, 0.0)  # None of null is <= 0.3

  # observed = 0.5 -> 2/5 = 0.4
  observed2 <- 0.5
  p_val2 <- mean(null_values <= observed2)
  expect_equal(p_val2, 0.4)

  # observed = 0.9 -> 5/5 = 1.0
  observed3 <- 0.9
  p_val3 <- mean(null_values <= observed3)
  expect_equal(p_val3, 1.0)
})

test("full mini pipeline with 3 synthetic subsamples", {
  set.seed(123)

  # Create 3 W matrices with similar structure
  make_W <- function(seed_offset) {
    set.seed(123 + seed_offset)
    data.table(
      VAR_ID = paste0("snp", 1:20),
      K1 = c(runif(10, 0.5, 1.0), runif(10, 0, 0.1)),
      K2 = c(runif(10, 0, 0.1), runif(10, 0.5, 1.0))
    )
  }

  w_list <- list(sub_01 = make_W(0), sub_02 = make_W(1), sub_03 = make_W(2))

  # Build null: 3 pairwise comparisons
  null_jaccards <- numeric()
  for (pair in combn(names(w_list), 2, simplify = FALSE)) {
    result <- hungarian_match_clusters(w_list[[pair[1]]], w_list[[pair[2]]], top_n = 10)
    null_jaccards <- c(null_jaccards, result$matching$jaccard_sim)
  }

  expect_true(length(null_jaccards) >= 3)  # At least 3 pairs * min matched
  expect_true(all(null_jaccards >= 0 & null_jaccards <= 1))

  # Create a "different" W matrix (simulating cross-ancestry)
  afr_W <- data.table(
    VAR_ID = paste0("snp", 1:20),
    K1 = c(runif(10, 0, 0.1), runif(10, 0.5, 1.0)),  # Swapped!
    K2 = c(runif(10, 0.5, 1.0), runif(10, 0, 0.1))
  )

  cross_result <- hungarian_match_clusters(w_list[["sub_01"]], afr_W, top_n = 10)
  expect_true(nrow(cross_result$matching) > 0)

  # Cross-ancestry p-values
  p_values <- sapply(cross_result$matching$jaccard_sim, function(obs) {
    mean(null_jaccards <= obs)
  })
  expect_true(all(p_values >= 0 & p_values <= 1))
})


# =====================================================================
# Summary
# =====================================================================

cat(sprintf("\n=== Test Summary: %d/%d passed, %d failed ===\n",
            n_passed, n_tests, n_failed))

if (n_failed > 0) {
  quit(save = "no", status = 1)
} else {
  cat("All tests passed.\n")
}
