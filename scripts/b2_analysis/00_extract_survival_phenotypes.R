#!/usr/bin/env Rscript
# 00_extract_survival_phenotypes.R
# Extract date-of-diagnosis for T2D and CAD from UKB SQLite DB.
# Build survival dataframe for time-to-comorbidity analysis.
#
# Fields used:
#   f130708  - Date of first T2D (E11) diagnosis
#   f131298  - Date of first CAD (I25) diagnosis
#   f40000   - Date of death
#   f41262   - Hospital episode dates (for last-record censoring)
#   f53      - Assessment centre date (for age-at-diagnosis)
#
# Usage:
#   Rscript scripts/b2_analysis/00_extract_survival_phenotypes.R
#   Rscript scripts/b2_analysis/00_extract_survival_phenotypes.R --config config/b2_config.yaml

library(RSQLite)
library(data.table)
library(yaml)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b2_config.yaml"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--config") {
    config_path <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== B2 Step 0: Extract Survival Phenotypes ===\n")
cat(sprintf("  Config: %s\n", config_path))
cat(sprintf("  DB:     %s\n", cfg$ukb$db_path))
cat(sprintf("  Output: %s\n\n", results_dir))

# --- Connect to UKB SQLite ---
con <- dbConnect(SQLite(), dbname = cfg$ukb$db_path)
on.exit(dbDisconnect(con), add = TRUE)

# --- Load withdrawn IDs ---
cat("--- Loading withdrawn participants ---\n")
withdrawn_ids <- fread(cfg$ukb$withdrawn_file, header = FALSE)$V1
cat(sprintf("  Withdrawn: %d participants\n\n", length(withdrawn_ids)))

# --- Query date fields ---
query_field <- function(field_table, col_name, aggregate = FALSE) {
  cat(sprintf("  Querying %s (%s)...", col_name, field_table))
  if (aggregate) {
    # For f41262: get max date per person (last hospital record)
    sql <- sprintf(
      "SELECT sample_id, MAX(pheno) AS %s FROM %s GROUP BY sample_id",
      col_name, field_table
    )
  } else {
    sql <- sprintf(
      "SELECT sample_id, pheno AS %s FROM %s WHERE instance = 0",
      col_name, field_table
    )
  }
  dt <- as.data.table(dbGetQuery(con, sql))
  cat(sprintf(" %d records\n", nrow(dt)))
  dt
}

cat("--- Querying UKB date fields ---\n")
t2d_dates  <- query_field(cfg$fields$t2d_date, "T2D_date")
cad_dates  <- query_field(cfg$fields$cad_date, "CAD_date")
death_dates <- query_field(cfg$fields$death_date, "death_date")
last_hosp   <- query_field(cfg$fields$hosp_dates, "last_hosp_date", aggregate = TRUE)
assess_dates <- query_field(cfg$fields$assessment_date, "assessment_date")

# --- Convert to Date type ---
cat("\n--- Converting to Date format ---\n")
t2d_dates[, T2D_date := as.Date(T2D_date)]
cad_dates[, CAD_date := as.Date(CAD_date)]
death_dates[, death_date := as.Date(death_date)]
last_hosp[, last_hosp_date := as.Date(last_hosp_date)]
assess_dates[, assessment_date := as.Date(assessment_date)]

# --- Merge all date fields ---
cat("--- Merging date fields ---\n")
surv <- merge(t2d_dates, cad_dates, by = "sample_id", all = TRUE)
surv <- merge(surv, death_dates, by = "sample_id", all.x = TRUE)
surv <- merge(surv, last_hosp, by = "sample_id", all.x = TRUE)
surv <- merge(surv, assess_dates, by = "sample_id", all.x = TRUE)
cat(sprintf("  Total individuals with T2D or CAD date: %d\n", nrow(surv)))

# --- Exclude withdrawn ---
n_before <- nrow(surv)
surv <- surv[!sample_id %in% withdrawn_ids]
cat(sprintf("  After excluding withdrawn: %d (removed %d)\n",
            nrow(surv), n_before - nrow(surv)))

# --- Filter: must have at least one diagnosis ---
surv <- surv[!is.na(T2D_date) | !is.na(CAD_date)]
cat(sprintf("  With at least one of T2D/CAD: %d\n", nrow(surv)))

# --- Determine first diagnosis ---
cat("\n--- Determining first diagnosis ---\n")
surv[, first_dx := fifelse(
  !is.na(T2D_date) & !is.na(CAD_date) & T2D_date == CAD_date, "dual",
  fifelse(
    !is.na(T2D_date) & (is.na(CAD_date) | T2D_date < CAD_date), "T2D",
    fifelse(
      !is.na(CAD_date) & (is.na(T2D_date) | CAD_date < T2D_date), "CAD",
      NA_character_
    )
  )
)]

# Tabulate
tab <- surv[, .N, by = first_dx]
cat("  First diagnosis counts:\n")
for (r in seq_len(nrow(tab))) {
  cat(sprintf("    %s: %d\n", tab$first_dx[r], tab$N[r]))
}

# Write tabulation before removing dual
fwrite(tab, file.path(results_dir, "cohort_tabulation.tsv"), sep = "\t")
cat(sprintf("  Saved: %s\n", file.path(results_dir, "cohort_tabulation.tsv")))

# --- Remove dual-diagnosed ---
n_dual <- sum(surv$first_dx == "dual", na.rm = TRUE)
cat(sprintf("\n  Removing %d dual-diagnosed (same-day) individuals\n", n_dual))
surv <- surv[first_dx != "dual"]
cat(sprintf("  Remaining: %d\n", nrow(surv)))

# --- First diagnosis date ---
surv[, first_dx_date := fifelse(first_dx == "T2D", T2D_date, CAD_date)]

# --- Second diagnosis (comorbidity event) ---
surv[, second_dx_date := fifelse(first_dx == "T2D", CAD_date, T2D_date)]
surv[, event := fifelse(!is.na(second_dx_date), 1L, 0L)]

# --- Censoring date ---
# Use min(death_date, last_hosp_date) for censored individuals
surv[, censor_date := pmin(death_date, last_hosp_date, na.rm = TRUE)]

# For individuals with no death or hospital record after first dx,
# use administrative end-of-follow-up (max of f41262 across all individuals)
admin_censor <- max(last_hosp$last_hosp_date, na.rm = TRUE)
cat(sprintf("\n  Administrative censoring date: %s\n", admin_censor))
surv[is.na(censor_date), censor_date := admin_censor]

# --- End date ---
surv[, end_date := fifelse(event == 1L, second_dx_date, censor_date)]

# --- Time in years ---
surv[, time_years := as.numeric(difftime(end_date, first_dx_date, units = "days")) / 365.25]

# --- Age at first diagnosis ---
# Load baseline age from covariates
covar <- fread(cfg$phenotypes$covariate_file)
covar_age <- covar[, .(IID, age)]
surv <- merge(surv, covar_age, by.x = "sample_id", by.y = "IID", all.x = TRUE)

# Compute age at first diagnosis
surv[, age_at_first_dx := fifelse(
  !is.na(assessment_date) & !is.na(age),
  age + as.numeric(difftime(first_dx_date, assessment_date, units = "days")) / 365.25,
  NA_real_
)]

# --- Sanity checks ---
cat("\n--- Sanity checks ---\n")
n_neg_time <- sum(surv$time_years <= 0, na.rm = TRUE)
n_na_time <- sum(is.na(surv$time_years))
cat(sprintf("  Rows with time_years <= 0: %d\n", n_neg_time))
cat(sprintf("  Rows with time_years NA: %d\n", n_na_time))

# Remove problematic rows
if (n_neg_time > 0 || n_na_time > 0) {
  surv <- surv[!is.na(time_years) & time_years > 0]
  cat(sprintf("  After removing: %d individuals remain\n", nrow(surv)))
}

# --- Summary statistics ---
cat("\n--- Summary ---\n")
cat(sprintf("  Total eligible individuals: %d\n", nrow(surv)))
cat(sprintf("  CAD-first: %d (events: %d, %.1f%%)\n",
            sum(surv$first_dx == "CAD"),
            sum(surv$first_dx == "CAD" & surv$event == 1),
            100 * mean(surv$event[surv$first_dx == "CAD"])))
cat(sprintf("  T2D-first: %d (events: %d, %.1f%%)\n",
            sum(surv$first_dx == "T2D"),
            sum(surv$first_dx == "T2D" & surv$event == 1),
            100 * mean(surv$event[surv$first_dx == "T2D"])))
cat(sprintf("  Median follow-up: %.1f years\n", median(surv$time_years)))
cat(sprintf("  Median age at first dx: %.1f years\n",
            median(surv$age_at_first_dx, na.rm = TRUE)))

# --- Add PLINK-format IDs ---
surv[, FID := sample_id]
surv[, IID := sample_id]

# --- Write output ---
out_cols <- c("FID", "IID", "sample_id", "T2D_date", "CAD_date", "death_date",
              "last_hosp_date", "assessment_date", "first_dx", "first_dx_date",
              "second_dx_date", "censor_date", "end_date", "time_years", "event",
              "age", "age_at_first_dx")
out_file <- file.path(results_dir, "survival_phenotypes.tsv")
fwrite(surv[, ..out_cols], out_file, sep = "\t")
cat(sprintf("\n  Saved: %s (%d rows)\n", out_file, nrow(surv)))

cat("\n=== Step 0 complete ===\n")
