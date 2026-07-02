#!/usr/bin/env Rscript
# b2_2_phewas.R
# PheWAS: logistic regression of META cluster PRS vs ICD10 binary phenotypes
# from UK Biobank. Each 3-character ICD10 block in selected chapters is tested
# against each cluster PRS, adjusted for age, age^2, sex, Batch, PC1-10.
#
# Usage:
#   Rscript scripts/b2_analysis/b2_2_phewas.R --config config/b2_2_config.yaml
#   Rscript scripts/b2_analysis/b2_2_phewas.R --config config/b2_2_config.yaml --hard-assignment

library(data.table)
library(yaml)
library(RSQLite)
library(parallel)

# =====================================================================
# 1. Setup
# =====================================================================

args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL
i <- 1
while (i <= length(args)) {
  if (args[i] == "--config") {
    config_path <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}
if (is.null(config_path)) stop("--config is required.")

cfg <- read_yaml(config_path)

hard_assignment <- "--hard-assignment" %in% args

if (hard_assignment) {
  cfg$b1_results$prs_file <- "results/b1_analysis/prs_hard_assignment/prs/cluster_prs_all.tsv"
  cfg$results_dir <- "results/b2_2_analysis/prs_hard_assignment"
  cat("  ** HARD ASSIGNMENT MODE **\n")
  cat(sprintf("  PRS file: %s\n", cfg$b1_results$prs_file))
  cat(sprintf("  Output dir: %s\n\n", cfg$results_dir))
}

results_dir    <- cfg$results_dir
clusters       <- cfg$analysis$clusters
sample_group   <- cfg$analysis$sample_group
icd10_chapters <- cfg$analysis$icd10_chapters
min_cases      <- cfg$analysis$min_cases
prs_cols       <- paste0("PRS_", clusters)
cluster_labels <- unlist(cfg$cluster_labels)

cat("=== B2.2 PheWAS ===\n")
cat(sprintf("  Sample group:    %s\n", sample_group))
cat(sprintf("  Clusters:        %s\n", paste(clusters, collapse = ", ")))
cat(sprintf("  ICD10 chapters:  %s\n", paste(icd10_chapters, collapse = ", ")))
cat(sprintf("  Min cases:       %d\n", min_cases))
cat(sprintf("  Covariates:      age, age2, sex, Batch (factor), PC1-PC10\n"))
cat(sprintf("  PRS:             z-scored (ORs per SD)\n\n"))

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# =====================================================================
# 2. Load PRS and covariates
# =====================================================================

cat("--- Loading PRS and covariates ---\n")

prs   <- fread(cfg$b1_results$prs_file)
covar <- fread(cfg$phenotypes$covariate_file)

prs <- prs[group == sample_group]
cat(sprintf("  PRS (%s): %d individuals\n", sample_group, nrow(prs)))
cat(sprintf("  Covariates: %d individuals\n", nrow(covar)))

# Merge PRS + covariates
dt <- merge(prs[, c("FID", "IID", prs_cols), with = FALSE],
            covar[, .(IID, age, age2, sex, PC1, PC2, PC3, PC4, PC5,
                       PC6, PC7, PC8, PC9, PC10, Batch)],
            by = "IID")

cat(sprintf("  Merged: %d individuals\n", nrow(dt)))

# Convert Batch to factor
dt[, Batch := as.factor(Batch)]

# Z-score each PRS
for (col in prs_cols) {
  dt[, (col) := scale(get(col))[, 1]]
}

# =====================================================================
# 3. Extract ICD10 phenotypes from UKB SQLite
# =====================================================================

cat("\n--- Extracting ICD10 codes from UKB database ---\n")

con <- dbConnect(SQLite(), dbname = cfg$ukb$db_path)

# Read withdrawn IDs
withdrawn_ids <- fread(cfg$ukb$withdrawn_file, header = FALSE)$V1
cat(sprintf("  Withdrawn participants: %d\n", length(withdrawn_ids)))

# Query ICD10 tables
icd10_list <- lapply(cfg$ukb$icd10_tables, function(tbl) {
  cat(sprintf("  Reading %s...", tbl))
  out <- as.data.table(dbGetQuery(con, sprintf("SELECT sample_id, pheno FROM %s", tbl)))
  cat(sprintf(" %d records\n", nrow(out)))
  out
})
dbDisconnect(con)

icd10_all <- rbindlist(icd10_list)
cat(sprintf("  Total ICD10 records: %d\n", nrow(icd10_all)))

# Truncate to 3-character blocks
icd10_all[, code3 := substr(pheno, 1, 3)]

# Filter to target chapters
icd10_all <- icd10_all[substr(code3, 1, 1) %in% icd10_chapters]
cat(sprintf("  Records in target chapters: %d\n", nrow(icd10_all)))

# Remove withdrawn
icd10_all <- icd10_all[!sample_id %in% withdrawn_ids]

# Unique (sample_id, code3) pairs
icd10_unique <- unique(icd10_all[, .(sample_id, code3)])
cat(sprintf("  Unique (individual, code) pairs: %d\n", nrow(icd10_unique)))
cat(sprintf("  Unique ICD10 3-char codes: %d\n", uniqueN(icd10_unique$code3)))

# Filter to codes with enough cases among EUR validation IIDs
eur_ids <- dt$IID
icd10_eur <- icd10_unique[sample_id %in% eur_ids]
code_counts <- icd10_eur[, .N, by = code3]
codes_pass <- code_counts[N >= min_cases, code3]
cat(sprintf("  Codes with >= %d cases in %s: %d\n",
            min_cases, sample_group, length(codes_pass)))

# Build wide binary matrix (only EUR validation, only passing codes)
icd10_wide <- dcast(icd10_eur[code3 %in% codes_pass],
                    sample_id ~ code3, fun.aggregate = length, fill = 0L)

# Clamp to 0/1
pheno_cols <- setdiff(colnames(icd10_wide), "sample_id")
for (col in pheno_cols) {
  set(icd10_wide, j = col, value = as.integer(icd10_wide[[col]] > 0))
}

# Prefix columns to avoid R naming issues (e.g., I10 -> icd_I10)
new_names <- paste0("icd_", pheno_cols)
setnames(icd10_wide, pheno_cols, new_names)
pheno_cols_prefixed <- new_names

cat(sprintf("  Phenotype matrix: %d individuals x %d codes\n",
            nrow(icd10_wide), length(pheno_cols_prefixed)))

# =====================================================================
# 4. Merge all data
# =====================================================================

cat("\n--- Merging PRS + covariates + ICD10 phenotypes ---\n")

# Join on IID = sample_id
icd10_wide[, IID := sample_id]
dt <- merge(dt, icd10_wide[, c("IID", pheno_cols_prefixed), with = FALSE],
            by = "IID")

cat(sprintf("  Final sample: %d individuals\n", nrow(dt)))

# Recount cases after merge
case_counts <- sapply(pheno_cols_prefixed, function(col) sum(dt[[col]] == 1))
codes_final <- pheno_cols_prefixed[case_counts >= min_cases]
cat(sprintf("  Codes with >= %d cases after merge: %d\n",
            min_cases, length(codes_final)))

# Report chapter breakdown
code_chapters <- substr(gsub("^icd_", "", codes_final), 1, 1)
chapter_tab <- table(code_chapters)
cat("\n  Codes per chapter:\n")
for (ch in sort(names(chapter_tab))) {
  cat(sprintf("    %s: %d codes\n", ch, chapter_tab[ch]))
}

# =====================================================================
# 5. Run logistic regressions
# =====================================================================

cat(sprintf("\n--- Running PheWAS: %d codes x %d PRS = %d tests ---\n",
            length(codes_final), length(prs_cols),
            length(codes_final) * length(prs_cols)))

covar_terms <- c("age", "age2", "sex", "Batch", paste0("PC", 1:10))

run_phewas_one_code <- function(pheno_col) {
  code_raw <- gsub("^icd_", "", pheno_col)
  results <- vector("list", length(prs_cols))

  for (k in seq_along(prs_cols)) {
    prs_col <- prs_cols[k]
    cluster_id <- clusters[k]

    n_cases <- sum(dt[[pheno_col]] == 1)
    n_controls <- sum(dt[[pheno_col]] == 0)

    rhs <- paste(c(prs_col, covar_terms), collapse = " + ")
    form <- as.formula(paste(pheno_col, "~", rhs))

    fit <- tryCatch(
      suppressWarnings(glm(form, data = dt, family = binomial)),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      results[[k]] <- data.table(
        icd10_code = code_raw, chapter = substr(code_raw, 1, 1),
        cluster = cluster_id, cluster_label = cluster_labels[cluster_id],
        n_cases = n_cases, n_controls = n_controls,
        beta = NA_real_, se = NA_real_, OR = NA_real_,
        CI_lower = NA_real_, CI_upper = NA_real_,
        p_value = NA_real_, converged = FALSE
      )
      next
    }

    converged <- fit$converged
    coef_summ <- summary(fit)$coefficients

    if (!(prs_col %in% rownames(coef_summ))) {
      results[[k]] <- data.table(
        icd10_code = code_raw, chapter = substr(code_raw, 1, 1),
        cluster = cluster_id, cluster_label = cluster_labels[cluster_id],
        n_cases = n_cases, n_controls = n_controls,
        beta = NA_real_, se = NA_real_, OR = NA_real_,
        CI_lower = NA_real_, CI_upper = NA_real_,
        p_value = NA_real_, converged = converged
      )
      next
    }

    prs_row <- coef_summ[prs_col, ]
    beta_val <- prs_row["Estimate"]
    se_val   <- prs_row["Std. Error"]

    results[[k]] <- data.table(
      icd10_code = code_raw,
      chapter = substr(code_raw, 1, 1),
      cluster = cluster_id,
      cluster_label = cluster_labels[cluster_id],
      n_cases = n_cases,
      n_controls = n_controls,
      beta = round(beta_val, 6),
      se = round(se_val, 6),
      OR = round(exp(beta_val), 4),
      CI_lower = round(exp(beta_val - 1.96 * se_val), 4),
      CI_upper = round(exp(beta_val + 1.96 * se_val), 4),
      p_value = prs_row["Pr(>|z|)"],
      converged = converged
    )
  }
  rbindlist(results)
}

n_cores <- min(detectCores(), 4)
cat(sprintf("  Using %d cores\n", n_cores))

t0 <- proc.time()
all_results <- rbindlist(mclapply(codes_final, run_phewas_one_code,
                                   mc.cores = n_cores))
elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("  Completed in %.1f seconds\n", elapsed))

# =====================================================================
# 6. Multiple testing correction and output
# =====================================================================

cat("\n--- Results summary ---\n")

# Remove rows where model failed entirely
valid <- all_results[!is.na(p_value)]
cat(sprintf("  Total tests: %d (of %d attempted)\n", nrow(valid), nrow(all_results)))

# Multiple testing corrections
n_tests <- nrow(valid)
valid[, p_bonferroni := pmin(p_value * n_tests, 1)]
valid[, p_fdr := p.adjust(p_value, method = "BH")]

# Sort by p-value
valid <- valid[order(p_value)]

# Save full results
out_file <- file.path(results_dir, "phewas_results.csv")
fwrite(valid, out_file)
cat(sprintf("  Saved: %s (%d rows)\n", out_file, nrow(valid)))

# Summary statistics
n_nom    <- sum(valid$p_value < 0.05)
n_bonf   <- sum(valid$p_bonferroni < 0.05)
n_fdr    <- sum(valid$p_fdr < 0.05)
cat(sprintf("\n  Significant at p < 0.05 (nominal):     %d\n", n_nom))
cat(sprintf("  Significant at Bonferroni < 0.05:       %d\n", n_bonf))
cat(sprintf("  Significant at FDR < 0.05 (BH):         %d\n", n_fdr))

# Top 20 associations
cat("\n--- Top 20 associations ---\n")
top20 <- head(valid, 20)
cat(sprintf("%-6s %-20s %6s %6s %8s %15s %12s %10s\n",
            "Code", "Cluster", "cases", "ctrls", "OR", "95% CI", "p-value", "p_fdr"))
cat(paste(rep("-", 95), collapse = ""), "\n")
for (r in seq_len(nrow(top20))) {
  row <- top20[r]
  cat(sprintf("%-6s %-20s %6d %6d %8.3f (%6.3f-%6.3f) %12.2e %10.2e\n",
              row$icd10_code, row$cluster_label,
              row$n_cases, row$n_controls,
              row$OR, row$CI_lower, row$CI_upper,
              row$p_value, row$p_fdr))
}

# Per-chapter summary
cat("\n--- Per-chapter summary (FDR < 0.05) ---\n")
chapter_names <- c(
  E = "Endocrine/metabolic", F = "Mental/behavioural", G = "Nervous system",
  H = "Eyes/ears", I = "Circulatory", J = "Respiratory",
  K = "Digestive", L = "Skin", M = "Musculoskeletal", N = "Genitourinary"
)
for (ch in sort(unique(valid$chapter))) {
  ch_data <- valid[chapter == ch]
  n_sig <- sum(ch_data$p_fdr < 0.05)
  cat(sprintf("  %s (%s): %d codes tested, %d FDR-significant\n",
              ch, chapter_names[ch], uniqueN(ch_data$icd10_code), n_sig))
}

cat("\n=== B2.2 PheWAS complete ===\n")
