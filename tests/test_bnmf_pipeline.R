# test_bnmf_pipeline.R
# Integration tests: synthetic data -> preprocessing -> bNMF -> output files.

library(testthat)
library(data.table)

# Resolve project root (works whether run from project root or tests/)
project_root <- if (file.exists("scripts/analysis/bnmf_algorithm.R")) {
  getwd()
} else if (file.exists("../scripts/analysis/bnmf_algorithm.R")) {
  normalizePath("..")
} else {
  stop("Cannot find project root. Run from project root or tests/ directory.")
}

source(file.path(project_root, "scripts", "analysis", "bnmf_algorithm.R"))
source(file.path(project_root, "scripts", "analysis", "prep_bnmf_toy.R"))
source(file.path(project_root, "scripts", "analysis", "format_results.R"))
source(file.path(project_root, "tests", "generate_synthetic_data.R"))

# Ensure synthetic data exists
test_data_dir <- file.path(project_root, "tests", "test_data")
generate_synthetic_gwas(output_dir = test_data_dir)

# Test output directory
test_output_dir <- file.path(project_root, "tests", "test_output")
dir.create(test_output_dir, recursive = TRUE, showWarnings = FALSE)

context("Variant Filtering")

test_that("filter_variants produces non-empty output", {
  cad_file <- file.path(test_data_dir, "synthetic_CAD_EUR.txt.gz")
  filtered <- filter_variants(cad_file, p_threshold = 5e-8, maf_threshold = 0.001)

  expect_s3_class(filtered, "data.table")
  expect_gt(nrow(filtered), 0)
  expect_true("VAR_ID" %in% names(filtered))
  expect_true("BETA" %in% names(filtered))
  expect_true("SE" %in% names(filtered))
})

test_that("filter_variants removes low-MAF variants", {
  cad_file <- file.path(test_data_dir, "synthetic_CAD_EUR.txt.gz")
  filtered <- filter_variants(cad_file, p_threshold = 1, maf_threshold = 0.001)

  if ("MAF" %in% names(filtered)) {
    expect_true(all(filtered$MAF >= 0.001))
  }
})

test_that("filter_variants removes strand-ambiguous variants", {
  cad_file <- file.path(test_data_dir, "synthetic_CAD_EUR.txt.gz")
  filtered <- filter_variants(cad_file, p_threshold = 1, maf_threshold = 0)

  # Parse alleles from VAR_ID
  allele_pairs <- sapply(strsplit(filtered$VAR_ID, "_"), function(p) paste0(p[3], p[4]))
  ambiguous <- allele_pairs %in% c("AT", "TA", "CG", "GC")
  expect_equal(sum(ambiguous), 0)
})

context("Z-score Matrix Construction")

test_that("build_z_matrix produces correct dimensions", {
  cad_file <- file.path(test_data_dir, "synthetic_CAD_EUR.txt.gz")
  filtered <- filter_variants(cad_file, p_threshold = 5e-8, maf_threshold = 0.001)

  trait_files <- list(
    T2D = file.path(test_data_dir, "synthetic_T2D_EUR.txt.gz"),
    LDL = file.path(test_data_dir, "synthetic_LDL_EUR.txt.gz"),
    BMI = file.path(test_data_dir, "synthetic_BMI_EUR.txt.gz")
  )

  matrices <- build_z_matrix(filtered, trait_files)

  expect_equal(nrow(matrices$z_matrix), nrow(filtered))
  expect_equal(ncol(matrices$z_matrix), length(trait_files))
  expect_equal(colnames(matrices$z_matrix), names(trait_files))
  expect_equal(rownames(matrices$z_matrix), filtered$VAR_ID)
})

test_that("build_z_matrix handles missing variants gracefully", {
  cad_file <- file.path(test_data_dir, "synthetic_CAD_EUR.txt.gz")
  filtered <- filter_variants(cad_file, p_threshold = 5e-8, maf_threshold = 0.001)

  trait_files <- list(
    T2D = file.path(test_data_dir, "synthetic_T2D_EUR.txt.gz")
  )

  matrices <- build_z_matrix(filtered, trait_files)

  # Should still produce a matrix even if some variants don't match
  expect_equal(nrow(matrices$z_matrix), nrow(filtered))
  expect_false(all(is.na(matrices$z_matrix)))
})

context("Non-negative Expansion")

test_that("expand_to_nonneg doubles columns and produces non-negative values", {
  z <- matrix(c(1.5, -0.5, 0, 2.0, -1.0, 0.3), nrow = 3, ncol = 2,
              dimnames = list(c("v1", "v2", "v3"), c("T2D", "LDL")))

  nonneg <- expand_to_nonneg(z)

  expect_equal(nrow(nonneg), 3)
  expect_equal(ncol(nonneg), 4)
  expect_true(all(nonneg >= 0))
  expect_equal(colnames(nonneg), c("T2D_pos", "T2D_neg", "LDL_pos", "LDL_neg"))

  # Check specific values
  expect_equal(nonneg["v1", "T2D_pos"], 1.5)
  expect_equal(nonneg["v1", "T2D_neg"], 0)
  expect_equal(nonneg["v2", "T2D_pos"], 0)
  expect_equal(nonneg["v2", "T2D_neg"], 0.5)
})

context("bNMF Algorithm")

test_that("BayesNMF.L2EU produces W and H matrices with correct dimensions", {
  set.seed(123)
  # Small synthetic non-negative matrix
  V <- matrix(abs(rnorm(50 * 6)), nrow = 50, ncol = 6)
  colnames(V) <- paste0("trait_", 1:6)

  result <- BayesNMF.L2EU(V, n.iter = 500, K = 5, K0 = 5, tol = 1e-5)

  expect_true(!is.null(result$W))
  expect_true(!is.null(result$H))
  expect_equal(nrow(result$W), 50)
  expect_equal(ncol(result$W), 5)
  expect_equal(ncol(result$H), 6)  # after removing zero-sum columns
  expect_true(all(result$W >= 0))
  expect_true(all(result$H >= 0))
  expect_true(result$K_converged >= 1)
  expect_true(result$K_converged <= 5)
})

test_that("run_bnmf produces results for each replicate", {
  set.seed(123)
  V <- matrix(abs(rnorm(30 * 4)), nrow = 30, ncol = 4)

  results <- run_bnmf(V, n_reps = 2, K = 4, K0 = 4, seed = 1, tolerance = 1e-5)

  expect_length(results, 2)
  expect_true(!is.null(results[[1]]$W))
  expect_true(!is.null(results[[2]]$W))
})

test_that("bNMF is reproducible with same seed", {
  V <- matrix(abs(rnorm(30 * 4)), nrow = 30, ncol = 4)

  res1 <- run_bnmf(V, n_reps = 1, K = 4, seed = 42, tolerance = 1e-5)
  res2 <- run_bnmf(V, n_reps = 1, K = 4, seed = 42, tolerance = 1e-5)

  expect_equal(res1[[1]]$K_converged, res2[[1]]$K_converged)
})

context("Results Formatting and Output Files")

test_that("full pipeline produces all expected output files", {
  # Run mini pipeline on synthetic data
  cad_file <- file.path(test_data_dir, "synthetic_CAD_EUR.txt.gz")
  filtered <- filter_variants(cad_file, p_threshold = 5e-8, maf_threshold = 0.001)

  trait_files <- list(
    T2D = file.path(test_data_dir, "synthetic_T2D_EUR.txt.gz"),
    LDL = file.path(test_data_dir, "synthetic_LDL_EUR.txt.gz"),
    BMI = file.path(test_data_dir, "synthetic_BMI_EUR.txt.gz")
  )

  matrices <- build_z_matrix(filtered, trait_files)
  nonneg <- expand_to_nonneg(matrices$z_matrix)

  # Run with minimal settings for speed
  results <- run_bnmf(nonneg, n_reps = 2, K = 5, K0 = 5, seed = 1, tolerance = 1e-4)

  # Summarize
  summary_result <- summarize_bnmf(
    results_list = results,
    trait_names = names(trait_files),
    variant_ids = rownames(nonneg),
    output_dir = test_output_dir
  )

  # Check output files
  expect_true(file.exists(file.path(test_output_dir, "W_matrix.tsv")))
  expect_true(file.exists(file.path(test_output_dir, "H_matrix.tsv")))
  expect_true(file.exists(file.path(test_output_dir, "run_summary.csv")))

  # Check W_matrix structure
  w <- fread(file.path(test_output_dir, "W_matrix.tsv"))
  expect_true("VAR_ID" %in% names(w))
  expect_gt(ncol(w), 1)

  # Check H_matrix structure
  h <- fread(file.path(test_output_dir, "H_matrix.tsv"))
  expect_true("Cluster" %in% names(h))

  # Generate summary report
  generate_summary_report(summary_result, test_output_dir)
  expect_true(file.exists(file.path(test_output_dir, "summary.txt")))

  # Generate heatmaps (may skip if pheatmap not installed)
  plot_heatmaps(summary_result$W, summary_result$H, test_output_dir)
})

# Cleanup
teardown({
  unlink(test_output_dir, recursive = TRUE)
})
