#!/usr/bin/env Rscript
# run_toy_pipeline.R
# End-to-end orchestrator for the toy bNMF analysis.
# Usage:
#   Rscript scripts/pipeline/run_toy_pipeline.R                         # Run with real data
#   Rscript scripts/pipeline/run_toy_pipeline.R --use-synthetic         # Run with synthetic data
#   Rscript scripts/pipeline/run_toy_pipeline.R --ancestry EUR          # Single ancestry
#   Rscript scripts/pipeline/run_toy_pipeline.R --use-synthetic --nreps 3 --K 8 --phi 2.0
#   Rscript scripts/pipeline/run_toy_pipeline.R --ancestry AFR --p-threshold 1e-5

library(data.table)

# Parse command-line args
args <- commandArgs(trailingOnly = TRUE)
use_synthetic <- "--use-synthetic" %in% args
ancestry_arg <- NULL
if ("--ancestry" %in% args) {
  ancestry_arg <- args[which(args == "--ancestry") + 1]
}
n_reps <- 10
if ("--nreps" %in% args) {
  n_reps <- as.integer(args[which(args == "--nreps") + 1])
}
K_init <- 20
if ("--K" %in% args) {
  K_init <- as.integer(args[which(args == "--K") + 1])
}
p_threshold <- 5e-8
if ("--p-threshold" %in% args) {
  p_threshold <- as.numeric(args[which(args == "--p-threshold") + 1])
}
phi_value <- 1.0
if ("--phi" %in% args) {
  phi_value <- as.numeric(args[which(args == "--phi") + 1])
}

# Determine project root (assumes script is run from project root or via Rscript)
project_root <- getwd()
if (!file.exists(file.path(project_root, "scripts", "a1_analysis", "bnmf_algorithm.R"))) {
  # Try finding it relative to this script's location
  script_dir <- dirname(sys.frame(1)$ofile %||% ".")
  project_root <- normalizePath(file.path(script_dir, "..", ".."))
}

# Source analysis scripts (now in a1_analysis/)
source(file.path(project_root, "scripts", "a1_analysis", "bnmf_algorithm.R"))
source(file.path(project_root, "scripts", "a1_analysis", "prep_bnmf.R"))
source(file.path(project_root, "scripts", "a1_analysis", "format_results.R"))

# ===== CONFIGURATION =====

if (use_synthetic) {
  cat("=== RUNNING WITH SYNTHETIC DATA ===\n\n")

  # Generate synthetic data if it doesn't exist
  test_data_dir <- file.path(project_root, "tests", "test_data")
  if (!file.exists(file.path(test_data_dir, "synthetic_CAD_EUR.txt.gz"))) {
    source(file.path(project_root, "tests", "generate_synthetic_data.R"))
    generate_synthetic_gwas(output_dir = test_data_dir)
  }

  sumstats_dir <- test_data_dir

  ref_gwas <- list(
    EUR = c(
      file.path(sumstats_dir, "synthetic_CAD_EUR.txt.gz"),
      file.path(sumstats_dir, "synthetic_T2D_EUR.txt.gz")
    ),
    AFR = c(
      file.path(sumstats_dir, "synthetic_CAD_AFR.txt.gz"),
      file.path(sumstats_dir, "synthetic_T2D_AFR.txt.gz")
    )
  )

  trait_gwas <- list(
    EUR = list(
      LDL    = file.path(sumstats_dir, "synthetic_LDL_EUR.txt.gz"),
      Stroke = file.path(sumstats_dir, "synthetic_Stroke_EUR.txt.gz"),
      BMI    = file.path(sumstats_dir, "synthetic_BMI_EUR.txt.gz"),
      HF     = file.path(sumstats_dir, "synthetic_HF_EUR.txt.gz")
    ),
    AFR = list(
      LDL    = file.path(sumstats_dir, "synthetic_LDL_AFR.txt.gz"),
      Stroke = file.path(sumstats_dir, "synthetic_Stroke_AFR.txt.gz"),
      BMI    = file.path(sumstats_dir, "synthetic_BMI_AFR.txt.gz"),
      HF     = file.path(sumstats_dir, "synthetic_HF_AFR.txt.gz")
    )
  )
} else {
  cat("=== RUNNING WITH REAL DATA ===\n\n")
  sumstats_dir <- file.path(project_root, "sumstats")
  harmonized_dir <- file.path(sumstats_dir, "harmonized")

  ref_gwas <- list(
    EUR = c(
      file.path(harmonized_dir, "Tcheandjieu_NatureMed_2023.CAD.EUR.GRCh37.txt.gz"),
      file.path(harmonized_dir, "Suzuki_Nature_2024.t2d.EUR.GRCh37.txt.gz")
    ),
    AFR = c(
      file.path(harmonized_dir, "Tcheandjieu_NatureMed_2023.CAD.AFR.GRCh37.txt.gz"),
      file.path(harmonized_dir, "Suzuki_Nature_2024.t2d.AFR.GRCh37.txt.gz")
    ),
    META = c(
      file.path(harmonized_dir, "Tcheandjieu_NatureMed_2023.CAD.META.GRCh37.processed.txt.gz"),
      file.path(harmonized_dir, "Suzuki_Nature_2024.t2d.META.GRCh37.processed.txt.gz")
    )
  )

  trait_gwas <- list(
    EUR = list(
      LDL    = file.path(harmonized_dir, "Graham_Nature_2021.LDL.EUR.GRCh37.txt.gz"),
      Stroke = file.path(harmonized_dir, "Mishra_Nature_2022.stroke.EUR.GRCh37.harmonized.tsv.gz"),
      BMI    = file.path(harmonized_dir, "Huang_NatureComms_2022.BMI.EUR.GRCh37.txt.gz"),
      HF     = file.path(harmonized_dir, "Henry_NatureGenetics_2025.HF.EUR.GRCh37.txt.gz"),
      SBP    = file.path(harmonized_dir, "Keaton_NatureGenetics_2024.SBP.EUR.GRCh37.harmonized.tsv.gz")
    ),
    AFR = list(
      LDL    = file.path(harmonized_dir, "Graham_Nature_2021.LDL.AFR.GRCh37.txt.gz"),
      Stroke = file.path(harmonized_dir, "Mishra_Nature_2022.stroke.AFR.GRCh37.harmonized.tsv.gz"),
      BMI    = file.path(harmonized_dir, "Huang_NatureComms_2022.BMI.AFR.GRCh37.txt.gz"),
      HF     = file.path(harmonized_dir, "Joseph_NatureComms_2022.HF.AFR.GRCh37.txt.gz"),
      SBP    = file.path(harmonized_dir, "Singh_NatureComms_2023.SBP.AFR.GRCh37.harmonized.tsv.gz")
    ),
    META = list(
      LDL    = file.path(harmonized_dir, "Graham_Nature_2021.LDL.META.GRCh37.processed.txt.gz"),
      Stroke = file.path(harmonized_dir, "Mishra_Nature_2022.stroke.META.GRCh37.processed.txt.gz"),
      BMI    = file.path(harmonized_dir, "Huang_NatureComms_2022.BMI.META.GRCh37.processed.txt.gz"),
      HF     = file.path(harmonized_dir, "Henry_NatureGenetics_2025.HF.META.GRCh37.processed.txt.gz"),
      SBP    = file.path(harmonized_dir, "Giri_NatureGenetics_2018.SBP.META.GRCh37.processed.txt.gz")
    )
  )
}

ancestries <- c("EUR", "AFR", "META")
if (!is.null(ancestry_arg)) {
  if (!ancestry_arg %in% ancestries) {
    stop(sprintf("Invalid ancestry: %s. Must be one of: %s",
                 ancestry_arg, paste(ancestries, collapse = ", ")))
  }
  ancestries <- ancestry_arg
}
traits <- c("LDL", "Stroke", "BMI", "HF", "SBP")

# LD clumping configuration: ancestry-matched reference panels
ref_panel_map <- list(
  EUR  = file.path(project_root, "reference", "1kg_eur", "1000G.EUR.QC"),
  AFR  = file.path(project_root, "reference", "1kg_afr", "1000G.AFR.QC"),
  META = file.path(project_root, "reference", "1kg_eur", "1000G.EUR.QC")
)
clump_r2 <- 0.05
clump_kb <- 500

# ===== PIPELINE =====

for (anc in ancestries) {
  cat(sprintf("\n%s\n=== Processing ancestry: %s ===\n%s\n",
              strrep("=", 50), anc, strrep("=", 50)))

  output_dir <- file.path(project_root, "results", anc)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # 1. Validate input files
  cat("\n--- Step 1: Validating input files ---\n")
  ref_files <- ref_gwas[[anc]]
  skip_ancestry <- FALSE
  for (rf in ref_files) {
    validation <- validate_gwas_format(rf)
    if (!validation$valid) {
      cat(sprintf("VALIDATION FAILED for %s:\n  %s\n",
                  basename(rf), paste(validation$messages, collapse = "\n  ")))
      cat("Skipping this ancestry.\n")
      skip_ancestry <- TRUE
      break
    }
    cat(sprintf("Reference GWAS validation passed: %s\n", basename(rf)))
  }
  if (skip_ancestry) next

  # 2. QC and filter variants (union from all reference GWAS)
  cat("\n--- Step 2: Variant QC and filtering ---\n")

  # Determine LD reference panel for this ancestry
  ref_prefix <- ref_panel_map[[anc]]
  if (!is.null(ref_prefix) && !file.exists(paste0(ref_prefix, ".1.bed"))) {
    cat(sprintf("WARNING: LD reference panel not found at %s\n", ref_prefix))
    cat("LD clumping will be SKIPPED. Run scripts/pipeline/download_1kg_reference.sh first.\n")
    ref_prefix <- NULL
  }

  qc_result <- qc_variants_multi(ref_files, p_threshold = p_threshold,
                                  maf_threshold = 0.001,
                                  ref_panel_prefix = ref_prefix,
                                  clump_r2 = clump_r2, clump_kb = clump_kb)
  filtered <- qc_result$data

  if (nrow(filtered) < 10) {
    cat(sprintf("WARNING: Only %d variants passed QC. Need at least 10.\n", nrow(filtered)))
    cat("Skipping this ancestry.\n")
    next
  }

  # Save QC report
  qc_report_file <- file.path(output_dir, sprintf("qc_report_%s.tsv", anc))
  fwrite(qc_result$qc_report, qc_report_file, sep = "\t")
  cat(sprintf("Saved QC report: %s\n", qc_report_file))

  # Save filtered variants
  filtered_file <- file.path(output_dir, sprintf("filtered_variants_%s.tsv", anc))
  fwrite(filtered, filtered_file, sep = "\t")
  cat(sprintf("Saved filtered variants: %s\n", filtered_file))

  # 3. Build Z-score matrix
  cat("\n--- Step 3: Building Z-score matrix ---\n")
  trait_files <- trait_gwas[[anc]]
  matrices <- build_z_matrix(filtered, trait_files)

  cat(sprintf("Z-matrix dimensions: %d variants x %d traits\n",
              nrow(matrices$z_matrix), ncol(matrices$z_matrix)))

  # 4. Expand to non-negative
  cat("\n--- Step 4: Expanding to non-negative matrix ---\n")
  nonneg_matrix <- expand_to_nonneg(matrices$z_matrix)
  cat(sprintf("Non-negative matrix: %d variants x %d columns\n",
              nrow(nonneg_matrix), ncol(nonneg_matrix)))

  # Save the prepared matrix
  prep_file <- file.path(output_dir, sprintf("prepared_matrix_%s.tsv", anc))
  fwrite(data.table(VAR_ID = rownames(nonneg_matrix), nonneg_matrix),
         prep_file, sep = "\t")
  cat(sprintf("Saved prepared matrix: %s\n", prep_file))

  # 5. Run bNMF
  cat(sprintf("\n--- Step 5: Running bNMF (n_reps=%d, K=%d, phi=%.1f) ---\n", n_reps, K_init, phi_value))
  results <- run_bnmf(nonneg_matrix, n_reps = n_reps, K = K_init, K0 = K_init, seed = 42, phi = phi_value)

  # 6. Summarize and format results
  cat("\n--- Step 6: Formatting results ---\n")
  summary_result <- summarize_bnmf(
    results_list = results,
    trait_names = traits,
    variant_ids = rownames(nonneg_matrix),
    output_dir = output_dir,
    ancestry = anc
  )

  # 7. Generate heatmaps
  cat("\n--- Step 7: Generating heatmaps ---\n")
  plot_heatmaps(summary_result$W, summary_result$H, output_dir, ancestry = anc)

  # 8. Generate summary report
  cat("\n--- Step 8: Writing summary report ---\n")
  generate_summary_report(summary_result, output_dir, ancestry = anc)

  cat(sprintf("\n=== Completed %s. Results in: %s ===\n", anc, output_dir))
}

cat("\n=== Pipeline complete! ===\n")
cat("Output directories:\n")
for (anc in ancestries) {
  out <- file.path(project_root, "results", anc)
  if (dir.exists(out)) {
    cat(sprintf("  %s: %s\n", anc, out))
    cat(sprintf("    Files: %s\n", paste(list.files(out), collapse = ", ")))
  }
}
