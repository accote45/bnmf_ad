# bnmf_algorithm.R
# Core Bayesian Non-negative Matrix Factorization with automatic relevance determination.
# Adapted from gwas-partitioning/bnmf-clustering (Broad Institute, BSD license).
# Reference: Tan, V.Y. & Fevotte, C. - Automatic Relevance Determination in NMF.

#' BayesNMF with L2 norm and Euclidean distance
#'
#' @param V0 Non-negative input matrix (variants x 2*traits)
#' @param n.iter Maximum number of iterations
#' @param a0 Prior parameter for automatic relevance determination
#' @param tol Convergence tolerance on lambda changes
#' @param K Number of clusters (columns in W, rows in H)
#' @param K0 Initial number of clusters for prior calculation
#' @param phi Regularization scaling factor (multiplied by variance of V)
#' @return List with W, H, n.like, n.evid, n.lambda, n.error
BayesNMF.L2EU <- function(V0, n.iter = 10000, a0 = 10, tol = 1e-7,
                           K = 15, K0 = 15, phi = 1.0,
                           window_size = 25, min_iter = 100) {
  eps <- 1.e-50
  del <- 1.0

  # Remove zero-sum columns
  active_nodes <- colSums(V0) != 0
  V0 <- V0[, active_nodes]
  V <- V0 - min(V0)
  Vmax <- max(V)
  N <- dim(V)[1]
  M <- dim(V)[2]

  # Random initialization
  W <- matrix(runif(N * K) * Vmax, ncol = K)
  H <- matrix(runif(M * K) * Vmax, ncol = M)

  V.ap <- W %*% H + eps

  # Prior parameters
  phi <- sd(V)^2 * phi
  C <- (N + M) / 2 + a0 + 1
  b0 <- 3.14 * (a0 - 1) * mean(V) / (2 * K0)
  lambda.bound <- b0 / C
  lambda <- (0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / C
  lambda.cut <- lambda.bound * 1.5

  n.like <- list()
  n.evid <- list()
  n.error <- list()
  n.lambda <- list()
  n.lambda[[1]] <- lambda
  iter <- 2

  # Sliding-window convergence check on reconstruction error
  check_convergence <- function(errors, window_size) {
    n <- length(errors)
    if (n < window_size) return(FALSE)
    recent <- tail(errors, window_size)
    if (any(!is.finite(recent))) return(FALSE)
    mean_recent <- mean(recent)
    if (mean_recent == 0) return(TRUE)
    rel_change <- abs(mean(diff(recent)) / mean_recent)
    return(rel_change < tol)
  }

  converged <- FALSE
  lambda_converged <- FALSE
  error_converged <- FALSE

  while (!converged & iter < n.iter) {
    # Update H
    H <- H * (t(W) %*% V) /
      (t(W) %*% V.ap + phi * H * matrix(rep(1 / lambda, M), ncol = M) + eps)
    V.ap <- W %*% H + eps

    # Update W
    W <- W * (V %*% t(H)) /
      (V.ap %*% t(H) + phi * W * t(matrix(rep(1 / lambda, N), ncol = N)) + eps)
    V.ap <- W %*% H + eps

    # Update lambda (automatic relevance determination)
    lambda <- (0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / C
    del <- max(abs(lambda - n.lambda[[iter - 1]]) / n.lambda[[iter - 1]])

    like <- sum((V - V.ap)^2) / 2
    n.like[[iter]] <- like
    n.evid[[iter]] <- like + phi * sum((0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) /
                                         lambda + C * log(lambda))
    n.lambda[[iter]] <- lambda
    n.error[[iter]] <- sum((V - V.ap)^2)

    # Check convergence: require min_iter AND (lambda-based OR error-window-based)
    lambda_converged <- (del < tol)
    error_converged <- check_convergence(unlist(n.error), window_size)
    if (iter >= min_iter && (lambda_converged || error_converged)) {
      converged <- TRUE
    }

    if (iter %% 100 == 0) {
      cat(sprintf("iter %d | evid %.2f | like %.2f | error %.2f | del %.2e | active %d | K_eff %d\n",
                  iter, n.evid[[iter]], n.like[[iter]], n.error[[iter]], del,
                  sum(colSums(W) != 0), sum(lambda >= lambda.cut)))
    }
    iter <- iter + 1
  }

  if (converged) {
    conv_type <- if (lambda_converged && error_converged) "lambda+error"
                 else if (lambda_converged) "lambda"
                 else "error-window"
    cat(sprintf("Converged at iteration %d (del=%.2e, method=%s). Effective K = %d\n",
                iter - 1, del, conv_type, sum(lambda >= lambda.cut)))
  } else {
    cat(sprintf("Did NOT converge after %d iterations (del=%.2e). Effective K = %d\n",
                iter - 1, del, sum(lambda >= lambda.cut)))
  }

  return(list(
    W = W,
    H = H,
    n.like = n.like,
    n.evid = n.evid,
    n.lambda = n.lambda,
    n.error = n.error,
    K_converged = sum(lambda >= lambda.cut),
    lambda = lambda,
    lambda.cut = lambda.cut,
    active_nodes = active_nodes,
    n.active = sum(colSums(W) > 0),
    iterations = iter - 1,
    converged = converged
  ))
}


#' Run bNMF multiple times from different random initializations
#'
#' @param z_matrix Non-negative z-score matrix (variants x 2*traits)
#' @param n_reps Number of independent runs
#' @param K Initial number of clusters
#' @param K0 Clusters for prior calculation
#' @param seed Random seed for reproducibility
#' @param tolerance Convergence tolerance
#' @return List of results from each replicate
run_bnmf <- function(z_matrix, n_reps = 5, K = 10, K0 = 10,
                     seed = 1, tolerance = 1e-7, phi = 1.0) {
  cat(sprintf("Running bNMF clustering (%d replicates, K=%d, phi=%.1f)...\n", n_reps, K, phi))

  set.seed(seed)

  results <- lapply(seq_len(n_reps), function(r) {
    cat(sprintf("\n=== Replicate %d/%d ===\n", r, n_reps))
    res <- BayesNMF.L2EU(V0 = z_matrix, K = K, K0 = K0, tol = tolerance, phi = phi)
    res
  })

  # Summary of K values across replicates
  k_values <- sapply(results, function(r) r$K_converged)
  cat(sprintf("\nK values across replicates: %s\n", paste(k_values, collapse = ", ")))
  cat(sprintf("Most frequent K: %d\n", as.integer(names(which.max(table(k_values))))))

  results
}


#' Test multiple phi values to find optimal regularization strength
#'
#' For each phi value, runs bNMF n_reps times and keeps the best run
#' (lowest final reconstruction error). Returns all best results for comparison.
#'
#' @param V0 Non-negative input matrix (variants x 2*traits)
#' @param phi_values Numeric vector of phi values to test
#' @param n_reps Number of replicates per phi value
#' @param seed Base random seed (each replicate uses seed + i)
#' @param ... Additional arguments passed to BayesNMF.L2EU (K, K0, tol, etc.)
#' @return Named list (keyed by phi value) with best_result, all_K_values, all_final_errors
test_phi_values <- function(V0, phi_values = c(1.0, 2.0, 5.0, 10.0),
                            n_reps = 10, seed = 1, ...) {
  cat(sprintf("Testing %d phi values with %d reps each...\n",
              length(phi_values), n_reps))

  results <- list()
  for (phi in phi_values) {
    cat(sprintf("\n--- phi = %.1f ---\n", phi))
    best_solution <- NULL
    best_error <- Inf
    k_values <- integer(n_reps)
    final_errors <- numeric(n_reps)

    for (i in seq_len(n_reps)) {
      set.seed(seed + i)
      cat(sprintf("  Rep %d/%d: ", i, n_reps))
      result <- BayesNMF.L2EU(V0, phi = phi, ...)

      err_list <- result$n.error
      final_error <- err_list[[length(err_list)]]

      if (!is.finite(final_error)) {
        cat(sprintf("WARNING: non-finite error (%.2e), skipping.\n", final_error))
        k_values[i] <- NA_integer_
        final_errors[i] <- NA_real_
        next
      }

      k_values[i] <- result$K_converged
      final_errors[i] <- final_error

      if (final_error < best_error) {
        best_solution <- result
        best_error <- final_error
      }
      cat(sprintf("K=%d, error=%.2f, converged=%s\n",
                  result$K_converged, final_error, result$converged))
    }

    results[[as.character(phi)]] <- list(
      best_result = best_solution,
      all_K_values = k_values,
      all_final_errors = final_errors,
      phi = phi
    )

    valid <- !is.na(k_values)
    if (any(valid)) {
      cat(sprintf("  phi=%.1f summary: K values=%s, median error=%.2f\n",
                  phi, paste(k_values[valid], collapse = ","),
                  median(final_errors[valid])))
    }
  }

  # Cross-phi summary
  cat("\n=== Phi comparison ===\n")
  for (phi_str in names(results)) {
    r <- results[[phi_str]]
    valid <- !is.na(r$all_final_errors)
    if (any(valid)) {
      cat(sprintf("phi=%s: best_error=%.2f, best_K=%d, modal_K=%d\n",
                  phi_str,
                  min(r$all_final_errors[valid]),
                  r$best_result$K_converged,
                  as.integer(names(which.max(table(r$all_K_values[valid]))))))
    }
  }

  results
}
