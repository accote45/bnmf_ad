#!/usr/bin/env Rscript
# 02_run_bnmf_subsample.R
# Run the full bNMF pipeline for a single EUR subsample.
# Variant selection from UKB subsample GWAS (genotyped variant space),
# Z-score matrix from published reference GWAS, bNMF factorization.
#
# Usage:
#   Rscript scripts/c1_ancestry_test/02_run_bnmf_subsample.R --subsample-id 1
#   Rscript scripts/c1_ancestry_test/02_run_bnmf_subsample.R --subsample-id 5 --config config/c1_config.yaml

library(data.table)
library(yaml)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)

subsample_id <- NULL
config_path  <- "config/c1_config.yaml"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--subsample-id") {
    subsample_id <- as.integer(args[i + 1]); i <- i + 2
  } else if (args[i] == "--config") {
    config_path <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.null(subsample_id)) {
  stop("Usage: Rscript 02_run_bnmf_subsample.R --subsample-id N [--config path]")
}

# --- Load config ---
cfg <- yaml::read_yaml(config_path)

# --- Determine project root ---
project_root <- getwd()
if (!file.exists(file.path(project_root, "scripts", "a1_analysis", "bnmf_algorithm.R"))) {
  script_dir <- dirname(sys.frame(1)$ofile %||% ".")
  project_root <- normalizePath(file.path(script_dir, "..", ".."))
}

# --- Source existing analysis modules ---
source(file.path(project_root, "scripts", "a1_analysis", "bnmf_algorithm.R"))
source(file.path(project_root, "scripts", "a1_analysis", "prep_bnmf.R"))
source(file.path(project_root, "scripts", "a1_analysis", "format_results.R"))

# --- Configuration ---
sub_label  <- sprintf("sub_%02d", subsample_id)
output_dir <- file.path(project_root, cfg$results_dir, "subsamples", sub_label)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_reps      <- cfg$bnmf$nreps       # 100
K_init      <- cfg$bnmf$K           # 20
phi_value   <- cfg$bnmf$phi         # 1.0
p_threshold <- as.numeric(cfg$bnmf$p_threshold_eur)  # 5e-8
clump_r2    <- cfg$ld_clump$r2      # 0.05
clump_kb    <- cfg$ld_clump$kb      # 500
hapmap3_file <- if (!is.null(cfg$hapmap3_snp_file)) {
  file.path(project_root, cfg$hapmap3_snp_file)
} else NULL

# LD reference panel (EUR)
ref_prefix <- file.path(project_root, cfg$ld_clump$ref_panel$EUR)
if (!file.exists(paste0(ref_prefix, ".1.bed"))) {
  cat(sprintf("WARNING: LD reference panel not found at %s\n", ref_prefix))
  cat("LD clumping will be SKIPPED.\n")
  ref_prefix <- NULL
}

# UKB subsample trait GWAS for variant selection (genotyped variant space)
traits <- cfg$gwas$traits
variant_gwas_files <- sapply(traits, function(tr) {
  file.path(output_dir, sprintf("%s_formatted.txt.gz", tr))
})

# Published reference GWAS for Z-score trait columns
z_trait_files <- setNames(
  lapply(cfg$ref_gwas$EUR, function(f) file.path(project_root, f)),
  names(cfg$ref_gwas$EUR)
)

cat(sprintf("\n%s\n=== C1.1 bNMF: EUR %s ===\n%s\n",
            strrep("=", 50), sub_label, strrep("=", 50)))
cat(sprintf("  nreps: %d, K: %d, phi: %.1f, p-threshold: %s\n",
            n_reps, K_init, phi_value, format(p_threshold, scientific = TRUE)))

# --- Step 1: Validate input files ---
cat("\n--- Step 1: Validating input files ---\n")
for (vf in variant_gwas_files) {
  if (!file.exists(vf)) {
    stop(sprintf("UKB trait GWAS not found: %s\nRun 01_run_gwas.sh first.", vf))
  }
  cat(sprintf("  Variant GWAS OK: %s\n", basename(vf)))
}

for (tr in names(z_trait_files)) {
  tf <- z_trait_files[[tr]]
  if (!file.exists(tf)) {
    stop(sprintf("Reference GWAS not found: %s", tf))
  }
  validation <- validate_gwas_format(tf)
  if (!validation$valid) {
    stop(sprintf("Validation failed for %s: %s",
                 basename(tf), paste(validation$messages, collapse = "; ")))
  }
  cat(sprintf("  Z-score GWAS OK: %s (%s)\n", tr, basename(tf)))
}

# --- Step 2: Variant QC (from UKB trait GWAS) ---
cat("\n--- Step 2: Variant QC (UKB trait GWAS) ---\n")
qc_result <- qc_variants_multi(
  variant_gwas_files,
  p_threshold     = p_threshold,
  maf_threshold   = cfg$bnmf$maf_threshold,
  ref_panel_prefix = ref_prefix,
  clump_r2        = clump_r2,
  clump_kb        = clump_kb,
  hapmap3_file    = hapmap3_file
)
filtered <- qc_result$data

cat(sprintf("  Variants after QC: %d\n", nrow(filtered)))

if (nrow(filtered) < 10) {
  stop(sprintf("Only %d variants passed QC. Need at least 10.", nrow(filtered)))
}

# Save QC report
qc_file <- file.path(output_dir, sprintf("qc_report_%s.tsv", sub_label))
fwrite(qc_result$qc_report, qc_file, sep = "\t")

filtered_file <- file.path(output_dir, sprintf("filtered_variants_%s.tsv", sub_label))
fwrite(filtered, filtered_file, sep = "\t")

# --- Step 3: Build Z-score matrix from published reference GWAS ---
cat("\n--- Step 3: Building Z-score matrix ---\n")
matrices <- build_z_matrix(filtered, z_trait_files)
cat(sprintf("  Z-matrix: %d variants x %d traits\n",
            nrow(matrices$z_matrix), ncol(matrices$z_matrix)))

# --- Step 4: Expand to non-negative ---
cat("\n--- Step 4: Non-negative expansion ---\n")
nonneg_matrix <- expand_to_nonneg(matrices$z_matrix)
cat(sprintf("  Non-negative matrix: %d x %d\n",
            nrow(nonneg_matrix), ncol(nonneg_matrix)))

prep_file <- file.path(output_dir, sprintf("prepared_matrix_%s.tsv", sub_label))
fwrite(data.table(VAR_ID = rownames(nonneg_matrix), nonneg_matrix),
       prep_file, sep = "\t")

# --- Step 5: Run bNMF ---
cat(sprintf("\n--- Step 5: Running bNMF (n_reps=%d, K=%d, phi=%.1f) ---\n",
            n_reps, K_init, phi_value))
results <- run_bnmf(
  nonneg_matrix,
  n_reps = n_reps,
  K      = K_init,
  K0     = K_init,
  seed   = 42,
  phi    = phi_value
)

# --- Step 6: Summarize results ---
cat("\n--- Step 6: Formatting results ---\n")
summary_result <- summarize_bnmf(
  results_list = results,
  trait_names  = names(z_trait_files),
  variant_ids  = rownames(nonneg_matrix),
  output_dir   = output_dir,
  ancestry     = sub_label
)

cat(sprintf("\n=== Completed EUR %s. Optimal K = %d ===\n",
            sub_label, summary_result$optimal_k))
cat(sprintf("  Results in: %s\n", output_dir))
