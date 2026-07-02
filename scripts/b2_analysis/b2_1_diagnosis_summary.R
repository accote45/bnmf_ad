#!/usr/bin/env Rscript
# b2_1_diagnosis_summary.R
# Tabulate first diagnosis frequency and inter-diagnosis interval statistics
# for EUR validation sample.
#
# Usage:
#   Rscript scripts/b2_analysis/b2_1_diagnosis_summary.R --config config/b2_1_config.yaml

library(data.table)
library(yaml)

args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL
i <- 1
while (i <= length(args)) {
  if (args[i] == "--config") { config_path <- args[i + 1]; i <- i + 2 }
  else { i <- i + 1 }
}
if (is.null(config_path)) stop("--config is required.")

cfg <- read_yaml(config_path)
results_dir  <- cfg$results_dir
sample_group <- cfg$analysis$sample_group

cat("=== B2.1 Diagnosis Summary (EUR Validation) ===\n\n")

# --- Load and merge ---
prs  <- fread(cfg$b1_results$prs_file)
surv <- fread(cfg$survival_phenotypes)
prs  <- prs[group == sample_group]

dt <- merge(prs[, .(FID, IID)],
            surv[, .(FID, IID, first_dx, time_years, event,
                     T2D_date, CAD_date, assessment_date, age, age_at_first_dx)],
            by = c("FID", "IID"))

cat(sprintf("  Total individuals with T2D or CAD: %d\n", nrow(dt)))

# --- Age at each disease's diagnosis (per disease, not just the first) ---
# age_at_disease = baseline age + years from assessment to that diagnosis.
# UK Biobank uses placeholder dates (e.g. 1902-02-02) for some records, which
# yield implausible ages; we drop ages outside [18, 100] before summarizing.
dt[, age_at_T2D := age + as.numeric(difftime(T2D_date, assessment_date, units = "days")) / 365.25]
dt[, age_at_CAD := age + as.numeric(difftime(CAD_date, assessment_date, units = "days")) / 365.25]
age_plausible <- function(x) !is.na(x) & x >= 18 & x <= 100

# --- First diagnosis frequency ---
cat("\n--- First Diagnosis Frequency ---\n")
first_dx_tab <- dt[, .N, by = first_dx][order(-N)]
first_dx_tab[, pct := round(100 * N / sum(N), 1)]
print(first_dx_tab)

# --- Among those who developed comorbidity ---
comorbid <- dt[event == 1]
cat(sprintf("\n--- Comorbidity Events ---\n"))
cat(sprintf("  Developed comorbidity: %d / %d (%.1f%%)\n",
            nrow(comorbid), nrow(dt), 100 * nrow(comorbid) / nrow(dt)))

comorbid_dx_tab <- comorbid[, .N, by = first_dx][order(-N)]
comorbid_dx_tab[, pct := round(100 * N / sum(N), 1)]
cat("\n  First diagnosis among comorbid individuals:\n")
print(comorbid_dx_tab)

# --- Inter-diagnosis interval statistics ---
cat("\n--- Inter-Diagnosis Interval (years) ---\n")

# Overall
cat("\n  Overall (all comorbid):\n")
cat(sprintf("    n:      %d\n", nrow(comorbid)))
cat(sprintf("    Mean:   %.2f\n", mean(comorbid$time_years)))
cat(sprintf("    Median: %.2f\n", median(comorbid$time_years)))
cat(sprintf("    SD:     %.2f\n", sd(comorbid$time_years)))
cat(sprintf("    IQR:    %.2f - %.2f\n",
            quantile(comorbid$time_years, 0.25),
            quantile(comorbid$time_years, 0.75)))
cat(sprintf("    Range:  %.2f - %.2f\n",
            min(comorbid$time_years), max(comorbid$time_years)))

# By first diagnosis
cat("\n  By first diagnosis:\n")
interval_stats <- list()
for (dx in unique(comorbid$first_dx)) {
  sub <- comorbid[first_dx == dx]
  second_dx <- if (dx == "T2D") "CAD" else "T2D"
  # Age at the second (comorbidity) diagnosis = age at first dx + interval.
  age_second <- sub$age_at_first_dx + sub$time_years
  age_second <- age_second[age_plausible(age_second)]
  cat(sprintf("\n    %s first -> %s (n=%d):\n", dx, second_dx, nrow(sub)))
  cat(sprintf("      Interval mean:   %.2f years\n", mean(sub$time_years)))
  cat(sprintf("      Interval median: %.2f years\n", median(sub$time_years)))
  cat(sprintf("      Interval SD:     %.2f\n", sd(sub$time_years)))
  cat(sprintf("      Interval IQR:    %.2f - %.2f\n",
              quantile(sub$time_years, 0.25),
              quantile(sub$time_years, 0.75)))
  cat(sprintf("      Progressed within  1 yr: %.1f%%\n", 100 * mean(sub$time_years <= 1)))
  cat(sprintf("      Progressed within  5 yr: %.1f%%\n", 100 * mean(sub$time_years <= 5)))
  cat(sprintf("      Progressed within 10 yr: %.1f%%\n", 100 * mean(sub$time_years <= 10)))
  cat(sprintf("      Age at second dx (%s): mean %.1f, median %.1f\n",
              second_dx, mean(age_second), median(age_second)))

  interval_stats[[dx]] <- data.table(
    first_dx = dx,
    second_dx = second_dx,
    n_total = sum(dt$first_dx == dx),
    n_comorbid = nrow(sub),
    comorbidity_rate_pct = round(100 * nrow(sub) / sum(dt$first_dx == dx), 2),
    interval_mean = round(mean(sub$time_years), 2),
    interval_median = round(median(sub$time_years), 2),
    interval_sd = round(sd(sub$time_years), 2),
    interval_q25 = round(quantile(sub$time_years, 0.25), 2),
    interval_q75 = round(quantile(sub$time_years, 0.75), 2),
    interval_min = round(min(sub$time_years), 2),
    interval_max = round(max(sub$time_years), 2),
    pct_within_1yr = round(100 * mean(sub$time_years <= 1), 1),
    pct_within_5yr = round(100 * mean(sub$time_years <= 5), 1),
    pct_within_10yr = round(100 * mean(sub$time_years <= 10), 1),
    age_at_second_dx_mean = round(mean(age_second), 1),
    age_at_second_dx_median = round(median(age_second), 1)
  )
}

# --- Censored individuals ---
censored <- dt[event == 0]
cat(sprintf("\n--- Censored (no comorbidity) ---\n"))
cat(sprintf("  n: %d (%.1f%%)\n", nrow(censored), 100 * nrow(censored) / nrow(dt)))
cat(sprintf("  Mean follow-up: %.2f years\n", mean(censored$time_years)))
cat(sprintf("  Median follow-up: %.2f years\n", median(censored$time_years)))

# --- Age at diagnosis (per disease) ---
# Across all individuals carrying each diagnosis, regardless of order.
cat("\n--- Age at Diagnosis (per disease) ---\n")
age_t2d <- dt$age_at_T2D[age_plausible(dt$age_at_T2D)]
age_cad <- dt$age_at_CAD[age_plausible(dt$age_at_CAD)]
n_t2d_raw <- sum(!is.na(dt$age_at_T2D))
n_cad_raw <- sum(!is.na(dt$age_at_CAD))
cat(sprintf("  Dropped implausible (outside [18,100]): T2D %d, CAD %d\n",
            n_t2d_raw - length(age_t2d), n_cad_raw - length(age_cad)))

age_stats <- rbindlist(list(
  data.table(disease = "T2D", n = length(age_t2d),
             age_mean = round(mean(age_t2d), 2), age_median = round(median(age_t2d), 2),
             age_sd = round(sd(age_t2d), 2),
             age_q25 = round(quantile(age_t2d, 0.25), 2),
             age_q75 = round(quantile(age_t2d, 0.75), 2)),
  data.table(disease = "CAD", n = length(age_cad),
             age_mean = round(mean(age_cad), 2), age_median = round(median(age_cad), 2),
             age_sd = round(sd(age_cad), 2),
             age_q25 = round(quantile(age_cad, 0.25), 2),
             age_q75 = round(quantile(age_cad, 0.75), 2))
))
for (r in seq_len(nrow(age_stats))) {
  a <- age_stats[r]
  cat(sprintf("  %s (n=%d): mean %.1f, median %.1f, SD %.1f, IQR %.1f-%.1f\n",
              a$disease, a$n, a$age_mean, a$age_median, a$age_sd, a$age_q25, a$age_q75))
}

# --- Directional balance among comorbid individuals ---
cat("\n--- Directional Balance (comorbid only) ---\n")
n_t2d_first <- nrow(comorbid[first_dx == "T2D"])
n_cad_first <- nrow(comorbid[first_dx == "CAD"])
cat(sprintf("  T2D-first -> CAD: %d (%.1f%%)\n",
            n_t2d_first, 100 * n_t2d_first / nrow(comorbid)))
cat(sprintf("  CAD-first -> T2D: %d (%.1f%%)\n",
            n_cad_first, 100 * n_cad_first / nrow(comorbid)))

# --- Requested headline stats ---
cat("\n--- Headline Temporal Stats ---\n")
med_t2d_cad <- median(comorbid[first_dx == "T2D", time_years])
med_cad_t2d <- median(comorbid[first_dx == "CAD", time_years])
cat(sprintf("  Median time T2D -> CAD: %.2f years\n", med_t2d_cad))
cat(sprintf("  Median time CAD -> T2D: %.2f years\n", med_cad_t2d))
cat(sprintf("  Mean age at T2D diagnosis: %.1f years\n", mean(age_t2d)))
cat(sprintf("  Mean age at CAD diagnosis: %.1f years\n", mean(age_cad)))

# --- Save CSVs ---
out <- rbindlist(interval_stats)
out_file <- file.path(results_dir, "diagnosis_interval_summary.csv")
fwrite(out, out_file)
cat(sprintf("\nSaved: %s\n", out_file))

age_file <- file.path(results_dir, "age_at_diagnosis_summary.csv")
fwrite(age_stats, age_file)
cat(sprintf("Saved: %s\n", age_file))
