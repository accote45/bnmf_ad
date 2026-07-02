#!/usr/bin/env Rscript
# run_tests.R
# Test runner for the bNMF toy pipeline.
# Usage: Rscript tests/run_tests.R

library(testthat)

cat("=== bNMF Toy Pipeline Test Suite ===\n\n")

# Ensure we're in the project root
if (!file.exists("scripts/analysis/bnmf_algorithm.R")) {
  if (file.exists("../scripts/analysis/bnmf_algorithm.R")) {
    setwd("..")
  } else {
    stop("Please run this script from the project root directory.")
  }
}

cat(sprintf("Working directory: %s\n\n", getwd()))

# Run individual test files
test_files <- list.files("tests", pattern = "^test_.*\\.R$", full.names = TRUE)
cat(sprintf("Found %d test files:\n", length(test_files)))
cat(paste(" ", test_files, collapse = "\n"), "\n\n")

results <- lapply(test_files, function(f) {
  cat(sprintf("--- Running %s ---\n", basename(f)))
  test_file(f, reporter = "summary")
})

# Combine results
all_results <- do.call(rbind, lapply(results, as.data.frame))

cat("\n=== Test Summary ===\n")
n_pass <- sum(all_results$passed)
n_fail <- sum(all_results$failed)
n_warn <- sum(all_results$warning)
n_skip <- sum(all_results$skipped)
cat(sprintf("Passed: %d | Failed: %d | Warnings: %d | Skipped: %d\n",
            n_pass, n_fail, n_warn, n_skip))

if (n_fail == 0) {
  cat("All tests passed!\n")
} else {
  cat("Some tests failed. See details above.\n")
  quit(status = 1)
}
