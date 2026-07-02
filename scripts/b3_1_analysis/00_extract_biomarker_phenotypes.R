#!/usr/bin/env Rscript
# 00_extract_biomarker_phenotypes.R
# Extract continuous biomarker and anthropometric traits from UKB SQLite,
# identify medication users, apply corrections, compute derived traits
# (eGFR, non-HDL, WHtR, PREVENT-ASCVD 10yr).
#
# Usage:
#   Rscript scripts/b3_1_analysis/00_extract_biomarker_phenotypes.R
#   Rscript scripts/b3_1_analysis/00_extract_biomarker_phenotypes.R --config config/b3_1_config.yaml

library(RSQLite)
library(data.table)
library(yaml)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b3_1_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== B3 Step 0: Extract Biomarker Phenotypes ===\n")
cat(sprintf("  Config:  %s\n", config_path))
cat(sprintf("  DB:      %s\n", cfg$ukb$db_path))
cat(sprintf("  Output:  %s\n\n", results_dir))

# --- Connect to UKB database ---
con <- dbConnect(SQLite(), dbname = cfg$ukb$db_path)
on.exit(dbDisconnect(con), add = TRUE)

# --- Helpers ---
get_baseline <- function(field_id, col_name) {
  tbl_name <- paste0("f", field_id)
  if (!dbExistsTable(con, tbl_name)) {
    warning(sprintf("Table '%s' not found — skipping", tbl_name))
    return(NULL)
  }
  query <- sprintf(
    "SELECT sample_id, CAST(pheno AS REAL) AS %s FROM %s WHERE instance = 0",
    col_name, tbl_name
  )
  as.data.table(dbGetQuery(con, query))
}

get_baseline_averaged <- function(field_id, col_name) {
  tbl_name <- paste0("f", field_id)
  if (!dbExistsTable(con, tbl_name)) {
    warning(sprintf("Table '%s' not found — skipping", tbl_name))
    return(NULL)
  }
  query <- sprintf(
    "SELECT sample_id, AVG(CAST(pheno AS REAL)) AS %s FROM %s WHERE instance = 0 GROUP BY sample_id",
    col_name, tbl_name
  )
  as.data.table(dbGetQuery(con, query))
}

# ============================================================
# 1. Get non-withdrawn participant IDs
# ============================================================
cat("--- Extracting participant IDs ---\n")
if (dbExistsTable(con, "Participant")) {
  ids <- as.data.table(dbGetQuery(con,
    "SELECT sample_id FROM Participant WHERE withdrawn = 0"))
} else {
  ids <- as.data.table(dbGetQuery(con,
    "SELECT DISTINCT sample_id FROM f31"))
}
cat(sprintf("  Non-withdrawn participants: %d\n", nrow(ids)))

# ============================================================
# 2. Extract direct trait fields
# ============================================================
cat("\n--- Extracting direct traits ---\n")

trait_tables <- list()
for (trait_name in names(cfg$traits$direct)) {
  field_id <- cfg$traits$direct[[trait_name]]$field
  dt <- get_baseline(field_id, trait_name)
  if (!is.null(dt)) {
    n_valid <- sum(!is.na(dt[[trait_name]]))
    cat(sprintf("  %s (f%s): %d non-NA out of %d\n",
                trait_name, field_id, n_valid, nrow(dt)))
    trait_tables[[trait_name]] <- dt
  }
}

# ============================================================
# 3. Extract averaged fields (SBP)
# ============================================================
cat("\n--- Extracting averaged traits ---\n")
for (trait_name in names(cfg$traits$averaged)) {
  field_id <- cfg$traits$averaged[[trait_name]]$field
  dt <- get_baseline_averaged(field_id, trait_name)
  if (!is.null(dt)) {
    n_valid <- sum(!is.na(dt[[trait_name]]))
    cat(sprintf("  %s (f%s, averaged): %d non-NA out of %d\n",
                trait_name, field_id, n_valid, nrow(dt)))
    trait_tables[[trait_name]] <- dt
  }
}

# ============================================================
# 4. Merge all traits
# ============================================================
cat("\n--- Merging trait tables ---\n")
pheno <- ids
for (nm in names(trait_tables)) {
  pheno <- merge(pheno, trait_tables[[nm]], by = "sample_id", all.x = TRUE)
}
cat(sprintf("  Merged phenotype table: %d individuals\n", nrow(pheno)))

# ============================================================
# 5. Identify medication users
# ============================================================
cat("\n--- Identifying medication users ---\n")

# -- Statins --
statin_codes <- as.character(cfg$medication$statin_codes)
statin_f20003 <- as.data.table(dbGetQuery(con, sprintf(
  "SELECT DISTINCT sample_id FROM f%s WHERE instance = 0 AND pheno IN (%s)",
  cfg$medication$treatment_codes,
  paste(sprintf("'%s'", statin_codes), collapse = ", ")
)))

# Self-reported cholesterol-lowering: f6177 (male) code 1, f6153 (female) code 1
statin_sr_m <- as.data.table(dbGetQuery(con, sprintf(
  "SELECT DISTINCT sample_id FROM f%s WHERE instance = 0 AND CAST(pheno AS INTEGER) = 1",
  cfg$medication$self_report_male
)))
statin_sr_f <- as.data.table(dbGetQuery(con, sprintf(
  "SELECT DISTINCT sample_id FROM f%s WHERE instance = 0 AND CAST(pheno AS INTEGER) = 1",
  cfg$medication$self_report_female
)))

statin_ids <- unique(rbind(statin_f20003, statin_sr_m, statin_sr_f))
cat(sprintf("  Statin users: %d\n", nrow(statin_ids)))

# -- Antihypertensives --
antihyp_sr_m <- as.data.table(dbGetQuery(con, sprintf(
  "SELECT DISTINCT sample_id FROM f%s WHERE instance = 0 AND CAST(pheno AS INTEGER) = 2",
  cfg$medication$self_report_male
)))
antihyp_sr_f <- as.data.table(dbGetQuery(con, sprintf(
  "SELECT DISTINCT sample_id FROM f%s WHERE instance = 0 AND CAST(pheno AS INTEGER) = 2",
  cfg$medication$self_report_female
)))
antihyp_ids <- unique(rbind(antihyp_sr_m, antihyp_sr_f))
cat(sprintf("  Antihypertensive users: %d\n", nrow(antihyp_ids)))

# -- Diabetes medications --
diabetes_codes <- as.character(cfg$medication$diabetes_med_codes)
dm_f20003 <- as.data.table(dbGetQuery(con, sprintf(
  "SELECT DISTINCT sample_id FROM f%s WHERE instance = 0 AND pheno IN (%s)",
  cfg$medication$treatment_codes,
  paste(sprintf("'%s'", diabetes_codes), collapse = ", ")
)))

# Self-reported insulin: f6177/f6153 code 3
dm_sr_m <- as.data.table(dbGetQuery(con, sprintf(
  "SELECT DISTINCT sample_id FROM f%s WHERE instance = 0 AND CAST(pheno AS INTEGER) = 3",
  cfg$medication$self_report_male
)))
dm_sr_f <- as.data.table(dbGetQuery(con, sprintf(
  "SELECT DISTINCT sample_id FROM f%s WHERE instance = 0 AND CAST(pheno AS INTEGER) = 3",
  cfg$medication$self_report_female
)))

# Diabetes diagnosed with medication: f2443 = 1
dm_diagnosed <- if (dbExistsTable(con, paste0("f", cfg$medication$diabetes_diagnosed))) {
  as.data.table(dbGetQuery(con, sprintf(
    "SELECT DISTINCT sample_id FROM f%s WHERE instance = 0 AND CAST(pheno AS INTEGER) = 1",
    cfg$medication$diabetes_diagnosed
  )))
} else {
  data.table(sample_id = integer(0))
}

dm_ids <- unique(rbind(dm_f20003, dm_sr_m, dm_sr_f, dm_diagnosed))
cat(sprintf("  Diabetes medication users: %d\n", nrow(dm_ids)))

# Add medication flags
pheno[, on_statin := as.integer(sample_id %in% statin_ids$sample_id)]
pheno[, on_antihypertensive := as.integer(sample_id %in% antihyp_ids$sample_id)]
pheno[, on_diabetes_med := as.integer(sample_id %in% dm_ids$sample_id)]

cat(sprintf("\n  Medication prevalence in cohort:\n"))
cat(sprintf("    Statins: %d (%.1f%%)\n",
            sum(pheno$on_statin), 100 * mean(pheno$on_statin)))
cat(sprintf("    Antihypertensives: %d (%.1f%%)\n",
            sum(pheno$on_antihypertensive), 100 * mean(pheno$on_antihypertensive)))
cat(sprintf("    Diabetes meds: %d (%.1f%%)\n",
            sum(pheno$on_diabetes_med), 100 * mean(pheno$on_diabetes_med)))

# ============================================================
# 6. Apply medication corrections
# ============================================================
cat("\n--- Applying medication corrections ---\n")

# LDL: divide by 0.7 for statin users
pheno[, LDL_corrected := fifelse(on_statin == 1L & !is.na(LDL), LDL / 0.7, LDL)]
cat(sprintf("  LDL corrected: %d statin users adjusted\n",
            sum(pheno$on_statin == 1L & !is.na(pheno$LDL))))

# Total cholesterol: divide by 0.7 for statin users (for non-HDL computation)
pheno[, TC_corrected := fifelse(on_statin == 1L & !is.na(total_cholesterol),
                                total_cholesterol / 0.7, total_cholesterol)]
cat(sprintf("  TC corrected: %d statin users adjusted\n",
            sum(pheno$on_statin == 1L & !is.na(pheno$total_cholesterol))))

# SBP: add 15 for antihypertensive users
pheno[, SBP_corrected := fifelse(on_antihypertensive == 1L & !is.na(SBP),
                                 SBP + 15, SBP)]
cat(sprintf("  SBP corrected: %d antihypertensive users adjusted\n",
            sum(pheno$on_antihypertensive == 1L & !is.na(pheno$SBP))))

# HbA1c: add 10.93 mmol/mol for diabetes medication users
pheno[, HbA1c_corrected := fifelse(on_diabetes_med == 1L & !is.na(HbA1c),
                                   HbA1c + 10.93, HbA1c)]
cat(sprintf("  HbA1c corrected: %d diabetes med users adjusted\n",
            sum(pheno$on_diabetes_med == 1L & !is.na(pheno$HbA1c))))

# Fasting glucose: add 1.0 mmol/L for diabetes medication users
pheno[, glucose_corrected := fifelse(on_diabetes_med == 1L & !is.na(fasting_glucose),
                                     fasting_glucose + 1.0, fasting_glucose)]
cat(sprintf("  Glucose corrected: %d diabetes med users adjusted\n",
            sum(pheno$on_diabetes_med == 1L & !is.na(pheno$fasting_glucose))))

# Triglycerides: divide by 0.85 for statin users
pheno[, triglycerides_corrected := fifelse(on_statin == 1L & !is.na(triglycerides),
                                           triglycerides / 0.85, triglycerides)]
cat(sprintf("  Triglycerides corrected: %d statin users adjusted\n",
            sum(pheno$on_statin == 1L & !is.na(pheno$triglycerides))))

# ============================================================
# 7. Compute derived traits
# ============================================================
cat("\n--- Computing derived traits ---\n")

# -- eGFR: CKD-EPI 2021 creatinine-cystatin C equation --
# Needs age and sex from covariate file or from DB
age_dt <- get_baseline(21022, "age")
sex_dt <- get_baseline(31, "sex")
pheno <- merge(pheno, age_dt, by = "sample_id", all.x = TRUE)
pheno <- merge(pheno, sex_dt, by = "sample_id", all.x = TRUE)
pheno[, age := as.numeric(age)]
pheno[, sex := as.integer(sex)]

compute_egfr_ckdepi2021 <- function(creat_umol, cys_c, age, sex) {
  # CKD-EPI 2021 combined creatinine-cystatin C equation (race-free)
  # sex: 0=female, 1=male (UKB field 31)
  scr <- creat_umol / 88.4  # convert umol/L to mg/dL

  is_female <- (sex == 0L)
  kappa <- fifelse(is_female, 0.7, 0.9)
  alpha <- fifelse(is_female, -0.219, -0.144)

  scr_kappa <- scr / kappa
  cys_0.8 <- cys_c / 0.8

  egfr <- 135 *
    pmin(scr_kappa, 1)^alpha *
    pmax(scr_kappa, 1)^(-0.544) *
    pmin(cys_0.8, 1)^(-0.323) *
    pmax(cys_0.8, 1)^(-0.778) *
    0.9961^age *
    fifelse(is_female, 0.963, 1)

  egfr
}

pheno[, eGFR := compute_egfr_ckdepi2021(creatinine, cystatin_c, age, sex)]
cat(sprintf("  eGFR: %d computed\n", sum(!is.na(pheno$eGFR))))

# -- non-HDL cholesterol (using corrected TC) --
pheno[, non_HDL_chol := TC_corrected - HDL]
cat(sprintf("  non-HDL cholesterol: %d computed\n", sum(!is.na(pheno$non_HDL_chol))))

# -- WHtR (waist-to-height ratio) --
pheno[, WHtR := waist_circumference / standing_height]
cat(sprintf("  WHtR: %d computed\n", sum(!is.na(pheno$WHtR))))

# -- PREVENT-ASCVD 10-year risk score (Khan et al. 2024, Circulation) --
# Extract smoking status
smoking_dt <- get_baseline(cfg$medication$smoking_status, "smoking_status")
if (!is.null(smoking_dt)) {
  pheno <- merge(pheno, smoking_dt, by = "sample_id", all.x = TRUE)
  pheno[, current_smoker := as.integer(smoking_status == 2)]
} else {
  pheno[, current_smoker := NA_integer_]
}
cat(sprintf("  Current smokers: %d\n", sum(pheno$current_smoker == 1, na.rm = TRUE)))

compute_prevent_ascvd_10yr <- function(age, tc, hdl, sbp, egfr,
                                        has_diabetes, current_smoker, bmi,
                                        on_statin, on_antihyp, sex) {
  # AHA PREVENT 2023 — 10-year total ASCVD risk (logistic model)
  # Khan et al. 2024, Circulation 149:e1144-e1156
  # sex: 0=female, 1=male (UKB field 31)
  #
  # Uses raw (uncorrected) SBP; treatment status handled via antihtn flag.
  # For statin users, statin=1 is assumed (they are on statins at baseline).
  # Returns risk as percentage (0-100).

  is_female <- (sex == 0L)
  statin <- on_statin
  antihtn <- on_antihyp

  # SBP: treated vs untreated. For UKB, if on_antihyp=1 then SBP is treated.
  treated_sbp   <- fifelse(antihtn == 1L, sbp, 0)
  untreated_sbp <- fifelse(antihtn == 1L, 0, sbp)

  # Variable transformations (centered and scaled)
  # UKB TC/HDL already in mmol/L — no mg/dL conversion needed
  age_10       <- (age - 55) / 10
  tc_hdlc_mmol <- (tc - hdl) - 3.5
  hdlc_mmol    <- (hdl - 1.3) / 0.3

  # SBP spline at 110
  min_sbp_110 <- fifelse(
    treated_sbp == 0,
    (pmin(untreated_sbp, 110) - 110) / 20,
    (pmin(treated_sbp, 110) - 110) / 20
  )
  max_sbp_110 <- fifelse(
    treated_sbp == 0,
    (pmax(untreated_sbp, 110) - 130) / 20,
    (pmax(treated_sbp, 110) - 130) / 20
  )

  # BMI spline at 30

  min_bmi_30 <- (pmin(bmi, 30) - 25) / 5
  max_bmi_30 <- (pmax(bmi, 30) - 30) / 5

  # eGFR spline at 60
  min_egfr_60 <- (pmin(egfr, 60) - 60) / -15
  max_egfr_60 <- (pmax(egfr, 60) - 90) / -15

  # Composite (interaction) terms
  max_sbp_110_treated_comp <- fifelse(
    treated_sbp == 0, 0,
    (pmax(treated_sbp, 110) - 130) / 20 * antihtn
  )
  tc_hdlc_mmol_treated_comp <- tc_hdlc_mmol * statin
  tc_hdlc_mmol_comp         <- age_10 * tc_hdlc_mmol
  hdlc_mmol_comp            <- age_10 * hdlc_mmol
  age_max_sbp_110_comp <- fifelse(
    treated_sbp == 0,
    age_10 * (pmax(untreated_sbp, 110) - 130) / 20,
    age_10 * (pmax(treated_sbp, 110) - 130) / 20
  )
  diabetes_comp    <- age_10 * has_diabetes
  cursmk_comp      <- age_10 * current_smoker
  max_bmi_30_comp  <- age_10 * max_bmi_30
  min_egfr_60_comp <- age_10 * min_egfr_60

  # Female log-odds
  lo_f <- -3.819975 +
    0.719883  * age_10 +
    0.1176967 * tc_hdlc_mmol +
   -0.151185  * hdlc_mmol +
   -0.0835358 * min_sbp_110 +
    0.3592852 * max_sbp_110 +
    0.8348585 * has_diabetes +
    0.4831078 * current_smoker +
    0.0       * min_bmi_30 +
    0.0       * max_bmi_30 +
    0.4864619 * min_egfr_60 +
    0.0397779 * max_egfr_60 +
    0.2265309 * antihtn +
   -0.0592374 * statin +
   -0.0395762 * max_sbp_110_treated_comp +
    0.0844423 * tc_hdlc_mmol_treated_comp +
   -0.0567839 * tc_hdlc_mmol_comp +
    0.0325692 * hdlc_mmol_comp +
   -0.1035985 * age_max_sbp_110_comp +
   -0.2417542 * diabetes_comp +
   -0.0791142 * cursmk_comp +
    0.0       * max_bmi_30_comp +
    0.1671492 * min_egfr_60_comp

  # Male log-odds
  lo_m <- -3.500655 +
    0.7099847 * age_10 +
    0.1658663 * tc_hdlc_mmol +
   -0.1144285 * hdlc_mmol +
   -0.2837212 * min_sbp_110 +
    0.3239977 * max_sbp_110 +
    0.7189597 * has_diabetes +
    0.3956973 * current_smoker +
    0.0       * min_bmi_30 +
    0.0       * max_bmi_30 +
    0.3690075 * min_egfr_60 +
    0.0203619 * max_egfr_60 +
    0.2036522 * antihtn +
   -0.0865581 * statin +
   -0.0322916 * max_sbp_110_treated_comp +
    0.114563  * tc_hdlc_mmol_treated_comp +
   -0.0300005 * tc_hdlc_mmol_comp +
    0.0232747 * hdlc_mmol_comp +
   -0.0927024 * age_max_sbp_110_comp +
   -0.2018525 * diabetes_comp +
   -0.0970527 * cursmk_comp +
    0.0       * max_bmi_30_comp +
   -0.1217081 * min_egfr_60_comp

  lo <- fifelse(is_female, lo_f, lo_m)
  risk <- exp(lo) / (1 + exp(lo)) * 100
  risk
}

# Diabetes status for PREVENT: on diabetes medication OR HbA1c >= 48 mmol/mol
pheno[, has_diabetes := as.integer(on_diabetes_med == 1L |
                                   (!is.na(HbA1c) & HbA1c >= 48))]

pheno[, PREVENT_ASCVD := compute_prevent_ascvd_10yr(
  age = age,
  tc = total_cholesterol,
  hdl = HDL,
  sbp = SBP,
  egfr = eGFR,
  has_diabetes = has_diabetes,
  current_smoker = current_smoker,
  bmi = BMI,
  on_statin = on_statin,
  on_antihyp = on_antihypertensive,
  sex = sex
)]
cat(sprintf("  PREVENT-ASCVD 10yr: %d computed\n", sum(!is.na(pheno$PREVENT_ASCVD))))

# -- PCE 10-year ASCVD risk (Goff et al. 2014, JACC) --
# Pooled Cohort Equations — Cox survival model, valid for ages 40-75.
# Expects cholesterol in mg/dL; UKB stores mmol/L, so convert × 38.67.
compute_pce_ascvd_10yr <- function(age, tc_mmol, hdl_mmol, sbp,
                                    antihtn, current_smoke, diabetes,
                                    sex) {
  # sex: 0=female, 1=male (UKB field 31)
  # All EUR validation individuals treated as white (race=0)
  tc  <- tc_mmol  * 38.67
  hdl <- hdl_mmol * 38.67

  is_female <- (sex == 0L)
  is_male   <- (sex == 1L)
  in_range  <- (age >= 40 & age <= 75)

  treated_sbp   <- fifelse(antihtn == 1L, sbp, 1)
  untreated_sbp <- fifelse(antihtn == 0L, sbp, 1)
  treat_flag    <- as.numeric(antihtn == 1L)
  untreat_flag  <- as.numeric(antihtn == 0L)

  # White female
  lo_wf <- -29.799 * log(age) +
    4.884 * (log(age))^2 +
    13.54 * log(tc) +
    -3.114 * log(age) * log(tc) +
    -13.578 * log(hdl) +
    3.149 * log(age) * log(hdl) +
    2.019 * log(treated_sbp) * treat_flag +
    1.957 * log(untreated_sbp) * untreat_flag +
    7.574 * current_smoke +
    -1.665 * log(age) * current_smoke +
    0.661 * diabetes - (-29.18)

  # White male
  lo_wm <- 12.344 * log(age) +
    11.853 * log(tc) +
    -2.664 * log(age) * log(tc) +
    -7.990 * log(hdl) +
    1.769 * log(age) * log(hdl) +
    1.797 * log(treated_sbp) * treat_flag +
    1.764 * log(untreated_sbp) * untreat_flag +
    7.837 * current_smoke +
    -1.795 * log(age) * current_smoke +
    0.658 * diabetes - 61.18

  lo <- fifelse(is_female, lo_wf, lo_wm)
  s0 <- fifelse(is_female, 0.9665, 0.9144)

  risk <- (1 - s0^exp(lo)) * 100
  risk <- fifelse(in_range, risk, NA_real_)
  risk
}

pheno[, PCE_ASCVD := compute_pce_ascvd_10yr(
  age = age,
  tc_mmol = total_cholesterol,
  hdl_mmol = HDL,
  sbp = SBP,
  antihtn = on_antihypertensive,
  current_smoke = current_smoker,
  diabetes = has_diabetes,
  sex = sex
)]
cat(sprintf("  PCE-ASCVD 10yr: %d computed\n", sum(!is.na(pheno$PCE_ASCVD))))

# ============================================================
# 8. Prepare output
# ============================================================
cat("\n--- Preparing output ---\n")

# Add FID/IID
pheno[, `:=`(FID = sample_id, IID = sample_id)]

# Define analysis-ready trait columns
analysis_traits <- c("Lpa", "CRP", "HbA1c_corrected", "glucose_corrected",
                     "eGFR", "non_HDL_chol", "triglycerides_corrected",
                     "LDL_corrected", "WHtR", "waist_circumference",
                     "BMI", "SBP_corrected", "PREVENT_ASCVD", "PCE_ASCVD")

# Keep FID, IID, analysis traits, medication flags, and raw age/sex for reference
out_cols <- c("FID", "IID", analysis_traits,
              "on_statin", "on_antihypertensive", "on_diabetes_med",
              "age", "sex")
pheno_out <- pheno[, ..out_cols]

# Write phenotype file
pheno_file <- file.path(results_dir, "biomarker_phenotypes.tsv")
fwrite(pheno_out, pheno_file, sep = "\t")
cat(sprintf("  Phenotypes: %s (%d rows)\n", pheno_file, nrow(pheno_out)))

# Write medication flags
med_cols <- c("FID", "IID", "on_statin", "on_antihypertensive", "on_diabetes_med",
              "has_diabetes", "current_smoker")
med_out <- pheno[, ..med_cols]
med_file <- file.path(results_dir, "medication_flags.tsv")
fwrite(med_out, med_file, sep = "\t")
cat(sprintf("  Medication flags: %s (%d rows)\n", med_file, nrow(med_out)))

# Write phenotype summary
trait_display <- c("Lpa (nmol/L)", "CRP (mg/L)", "HbA1c (mmol/mol, corrected)",
                   "Fasting glucose (mmol/L, corrected)", "eGFR (mL/min/1.73m2)",
                   "non-HDL cholesterol (mmol/L)", "Triglycerides (mmol/L, corrected)",
                   "LDL (mmol/L, corrected)", "WHtR", "Waist circumference (cm)",
                   "BMI (kg/m2)", "SBP (mmHg, corrected)", "PREVENT-ASCVD 10yr (%)",
                   "PCE-ASCVD 10yr (%)")

summary_list <- lapply(seq_along(analysis_traits), function(i) {
  vals <- pheno_out[[analysis_traits[i]]]
  vals <- vals[!is.na(vals)]
  data.table(
    trait = analysis_traits[i],
    display_name = trait_display[i],
    n = length(vals),
    mean = round(mean(vals), 4),
    sd = round(sd(vals), 4),
    median = round(median(vals), 4),
    min = round(min(vals), 4),
    max = round(max(vals), 4)
  )
})
summary_dt <- rbindlist(summary_list)

summary_file <- file.path(results_dir, "phenotype_summary.tsv")
fwrite(summary_dt, summary_file, sep = "\t")
cat(sprintf("  Summary: %s\n", summary_file))

cat("\n=== Phenotype Summary ===\n")
print(summary_dt, nrow = Inf)

cat("\n=== Phenotype extraction complete ===\n")
