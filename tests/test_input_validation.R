# test_input_validation.R
# Tests for validate_gwas_format() and input handling.

library(testthat)
library(data.table)

# Resolve project root (works whether run from project root or tests/)
project_root <- if (file.exists("scripts/analysis/prep_bnmf_toy.R")) {
  getwd()
} else if (file.exists("../scripts/analysis/prep_bnmf_toy.R")) {
  normalizePath("..")
} else {
  stop("Cannot find project root. Run from project root or tests/ directory.")
}

source(file.path(project_root, "scripts", "analysis", "prep_bnmf_toy.R"))
source(file.path(project_root, "tests", "generate_synthetic_data.R"))

# Ensure synthetic data exists
test_data_dir <- file.path(project_root, "tests", "test_data")
generate_synthetic_gwas(output_dir = test_data_dir)

context("GWAS Input Validation")

test_that("valid GWAS file passes validation", {
  result <- validate_gwas_format(file.path(test_data_dir, "synthetic_CAD_EUR.txt.gz"))
  expect_true(result$valid)
  expect_equal(result$messages, "All checks passed")
})

test_that("missing file is caught", {
  result <- validate_gwas_format("nonexistent_file.txt")
  expect_false(result$valid)
  expect_true(grepl("not found", result$messages[1]))
})

test_that("missing columns are detected", {
  # Create a file missing P_VALUE
  bad_dt <- data.table(
    VAR_ID = "1_100000_A_C",
    BETA = 0.05,
    SE = 0.01,
    N = 10000,
    MAF = 0.3
  )
  bad_file <- file.path(test_data_dir, "bad_missing_col.txt")
  fwrite(bad_dt, bad_file, sep = "\t")

  result <- validate_gwas_format(bad_file)
  expect_false(result$valid)
  expect_true(any(grepl("Missing columns.*P_VALUE", result$messages)))

  unlink(bad_file)
})

test_that("bad VAR_ID format is detected", {
  bad_dt <- data.table(
    VAR_ID = c("rs12345", "rs67890"),
    P_VALUE = c(1e-10, 0.5),
    BETA = c(0.1, 0.01),
    SE = c(0.01, 0.02),
    N = c(10000, 10000),
    MAF = c(0.3, 0.2)
  )
  bad_file <- file.path(test_data_dir, "bad_var_ids.txt")
  fwrite(bad_dt, bad_file, sep = "\t")

  result <- validate_gwas_format(bad_file)
  expect_false(result$valid)
  expect_true(any(grepl("VAR_ID format invalid", result$messages)))

  unlink(bad_file)
})

test_that("empty file is handled", {
  empty_file <- file.path(test_data_dir, "empty.txt")
  file.create(empty_file)

  result <- validate_gwas_format(empty_file)
  expect_false(result$valid)
  expect_true(any(grepl("empty", result$messages, ignore.case = TRUE)))

  unlink(empty_file)
})

test_that("gzipped file is supported", {
  result <- validate_gwas_format(file.path(test_data_dir, "synthetic_CAD_EUR.txt.gz"))
  expect_true(result$valid)
})

test_that("file with all-NA column is detected", {
  bad_dt <- data.table(
    VAR_ID = c("1_100000_A_C", "2_200000_T_G"),
    P_VALUE = c(NA_real_, NA_real_),
    BETA = c(0.1, 0.2),
    SE = c(0.01, 0.02),
    N = c(10000, 10000),
    MAF = c(0.3, 0.2)
  )
  bad_file <- file.path(test_data_dir, "bad_all_na.txt")
  fwrite(bad_dt, bad_file, sep = "\t")

  result <- validate_gwas_format(bad_file)
  expect_false(result$valid)
  expect_true(any(grepl("entirely NA", result$messages)))

  unlink(bad_file)
})
