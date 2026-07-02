#!/usr/bin/env Rscript
# 00_extract_binary_phenotypes.R
# Extract binary case/control phenotypes for 6 traits from UKB SQLite DB
# using ICD-9 and ICD-10 diagnosis codes (main + secondary).
#
# Traits: T2D, CAD, Angina, MI, Stroke, PAD
#
# ICD code definitions:
#   T2D:    ICD-10 E11*          ICD-9 250* (exclude type 1: 5th digit 1 or 3)
#   CAD:    ICD-10 I25*          ICD-9 414*
#   Angina: ICD-10 I20*          ICD-9 413*
#   MI:     ICD-10 I21*, I22*    ICD-9 410*
#   Stroke: ICD-10 I60*, I61*, I63*, I64*   ICD-9 430-436
#   PAD:    ICD-10 I739          ICD-9 4439
#
# Uses both main (f41202/f41203) and secondary (f41204/f41205) diagnosis tables.
# Excludes withdrawn participants.
#
# Usage:
#   Rscript scripts/a0_analysis/00_extract_binary_phenotypes.R
#   Rscript scripts/a0_analysis/00_extract_binary_phenotypes.R \
#     --db /path/to/ukb.db \
#     --covar /path/to/ukb.covar \
#     --withdrawn /path/to/withdrawn.txt \
#     --pheno-dir phenotypes \
#     --out-dir results/a0_analysis/prs_ct

library(RSQLite)
library(data.table)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)

db_path       <- "/sc/arion/projects/paul_oreilly/data/ukb/phenotype/ukb18177.db"
covar_path    <- "/sc/arion/projects/paul_oreilly/data/ukb/phenotype/ukb18177.covar"
withdrawn_path <- "/sc/arion/projects/paul_oreilly/data/ukb/withdrawn/withdraw18177_398_20240119.txt"
pheno_dir     <- "phenotypes"
out_dir       <- "results/a0_analysis/prs_ct"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--db") {
    db_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--covar") {
    covar_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--withdrawn") {
    withdrawn_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--pheno-dir") {
    pheno_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--out-dir") {
    out_dir <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

dir.create(pheno_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Binary Phenotype Extraction ===\n")
cat(sprintf("  DB:        %s\n", db_path))
cat(sprintf("  Covar:     %s\n", covar_path))
cat(sprintf("  Withdrawn: %s\n", withdrawn_path))
cat(sprintf("  Pheno dir: %s\n", pheno_dir))
cat(sprintf("  Out dir:   %s\n\n", out_dir))

# --- ICD code definitions ---
# Each entry: list of regex patterns for ICD-10 and ICD-9
icd_definitions <- list(
  T2D = list(
    icd10 = "^E11",
    # ICD-9 250*: include all, then exclude type 1 (5th digit 1 or 3)
    icd9  = "^250",
    icd9_exclude = "[13]$"
  ),
  CAD = list(
    icd10 = "^I25",
    icd9  = "^414"
  ),
  Angina = list(
    icd10 = "^I20",
    icd9  = "^413"
  ),
  MI = list(
    icd10 = "^I2[12]",
    icd9  = "^410"
  ),
  Stroke = list(
    icd10 = "^I6[0134]",
    icd9  = "^43[0-6]"
  ),
  PAD = list(
    icd10 = "^I739$",
    icd9  = "^4439$"
  )
)

cat("ICD code definitions:\n")
for (trait in names(icd_definitions)) {
  def <- icd_definitions[[trait]]
  cat(sprintf("  %s: ICD-10=%s, ICD-9=%s\n", trait, def$icd10, def$icd9))
}
cat("\n")

# --- Connect to UKB database ---
con <- dbConnect(SQLite(), dbname = db_path)
on.exit(dbDisconnect(con), add = TRUE)

# --- Load withdrawn IDs ---
cat("--- Loading withdrawn participants ---\n")
withdrawn_ids <- fread(withdrawn_path, header = FALSE)$V1
cat(sprintf("  Withdrawn: %d participants\n", length(withdrawn_ids)))

# --- Get all participant IDs (excluding withdrawn) ---
cat("\n--- Extracting participant IDs ---\n")
all_ids <- as.data.table(dbGetQuery(con,
  "SELECT sample_id FROM participant WHERE withdrawn = 0"))
# Also exclude from external withdrawal file (more up-to-date)
all_ids <- all_ids[!sample_id %in% withdrawn_ids]
cat(sprintf("  Non-withdrawn participants: %d\n", nrow(all_ids)))

# --- Query all ICD codes (vectorized, single pass per table) ---
cat("\n--- Querying ICD diagnosis tables ---\n")

icd_tables <- list(
  icd10_main      = "f41202",
  icd10_secondary = "f41204",
  icd9_main       = "f41203",
  icd9_secondary  = "f41205"
)

# Read each ICD table once, keeping only sample_id and pheno (the code)
icd_data <- lapply(names(icd_tables), function(tbl_name) {
  tbl <- icd_tables[[tbl_name]]
  cat(sprintf("  Reading %s (%s)...", tbl_name, tbl))
  dt <- as.data.table(dbGetQuery(con,
    sprintf("SELECT sample_id, pheno FROM %s", tbl)))
  cat(sprintf(" %d records\n", nrow(dt)))
  dt[, source := tbl_name]
  dt
})
names(icd_data) <- names(icd_tables)

# Combine ICD-10 and ICD-9 separately
icd10_all <- rbind(icd_data$icd10_main, icd_data$icd10_secondary)
icd9_all  <- rbind(icd_data$icd9_main, icd_data$icd9_secondary)

cat(sprintf("\n  Total ICD-10 records: %d\n", nrow(icd10_all)))
cat(sprintf("  Total ICD-9 records:  %d\n", nrow(icd9_all)))

# --- Identify cases for each trait ---
cat("\n--- Identifying cases ---\n")

case_ids <- lapply(names(icd_definitions), function(trait) {
  def <- icd_definitions[[trait]]

  # ICD-10 matches
  ids_10 <- icd10_all[grepl(def$icd10, pheno), unique(sample_id)]

  # ICD-9 matches (with optional exclusion for T2D type 1)
  ids_9 <- icd9_all[grepl(def$icd9, pheno)]
  if (!is.null(def$icd9_exclude)) {
    ids_9 <- ids_9[!grepl(def$icd9_exclude, pheno)]
  }
  ids_9 <- unique(ids_9$sample_id)

  # Union of ICD-10 and ICD-9 cases
  case_set <- unique(c(ids_10, ids_9))

  cat(sprintf("  %s: %d ICD-10 cases, %d ICD-9 cases, %d total unique cases\n",
              trait, length(ids_10), length(ids_9), length(case_set)))
  case_set
})
names(case_ids) <- names(icd_definitions)

# --- Build phenotype data.table ---
cat("\n--- Building phenotype table ---\n")

pheno_dt <- copy(all_ids)
pheno_dt[, FID := sample_id]
pheno_dt[, IID := sample_id]

for (trait in names(case_ids)) {
  pheno_dt[, (trait) := fifelse(sample_id %in% case_ids[[trait]], 1L, 0L)]
}

# Reorder columns
pheno_dt <- pheno_dt[, c("FID", "IID", names(icd_definitions)), with = FALSE]

cat(sprintf("  Total individuals: %d\n", nrow(pheno_dt)))
cat("\n  Case/control counts:\n")
summary_rows <- list()
for (trait in names(icd_definitions)) {
  n_case <- sum(pheno_dt[[trait]] == 1)
  n_ctrl <- sum(pheno_dt[[trait]] == 0)
  cat(sprintf("    %s: %d cases / %d controls (prevalence %.2f%%)\n",
              trait, n_case, n_ctrl, 100 * n_case / nrow(pheno_dt)))
  summary_rows[[trait]] <- data.table(
    trait = trait, n_cases = n_case, n_controls = n_ctrl,
    prevalence_pct = round(100 * n_case / nrow(pheno_dt), 2)
  )
}
summary_dt <- rbindlist(summary_rows)

# --- Extract covariates ---
cat("\n--- Extracting covariates ---\n")

# Age at assessment (instance 0)
get_baseline <- function(field_id, col_name) {
  tbl_name <- paste0("f", field_id)
  query <- sprintf(
    "SELECT sample_id, pheno AS %s FROM %s WHERE instance = 0",
    col_name, tbl_name
  )
  as.data.table(dbGetQuery(con, query))
}

age_dt <- get_baseline(21003, "age")
cat(sprintf("  Age (f21003): %d records\n", nrow(age_dt)))

sex_dt <- get_baseline(31, "sex")
cat(sprintf("  Sex (f31): %d records\n", nrow(sex_dt)))

# Read pre-computed PCs and batch
cat("  Reading covariate file...\n")
pcs <- fread(covar_path)
cat(sprintf("  PC file: %d rows, columns: %s\n",
            nrow(pcs), paste(head(colnames(pcs), 15), collapse = ", ")))

# Merge demographics
covar_dt <- merge(all_ids, age_dt, by = "sample_id", all.x = TRUE)
covar_dt <- merge(covar_dt, sex_dt, by = "sample_id", all.x = TRUE)

# Merge PCs and Batch
covar_dt <- merge(covar_dt,
                  pcs[, .(FID, IID, Batch, PC1, PC2, PC3, PC4, PC5,
                          PC6, PC7, PC8, PC9, PC10)],
                  by.x = "sample_id", by.y = "IID", all.x = FALSE)

covar_dt[, age := as.numeric(age)]
covar_dt[, sex := as.integer(sex)]
covar_dt[, `:=`(FID = sample_id, IID = sample_id, age2 = age^2)]
covar_dt <- covar_dt[, .(FID, IID, age, age2, sex,
                          PC1, PC2, PC3, PC4, PC5,
                          PC6, PC7, PC8, PC9, PC10, Batch)]

# Drop rows with missing covariates
n_before <- nrow(covar_dt)
covar_dt <- covar_dt[complete.cases(covar_dt[, .(age, sex)])]
cat(sprintf("  Covariate table: %d individuals (%d dropped for missing age/sex)\n",
            nrow(covar_dt), n_before - nrow(covar_dt)))

# --- Restrict phenotypes to individuals with covariates ---
pheno_dt <- pheno_dt[IID %in% covar_dt$IID]
cat(sprintf("  Phenotype table (with covariates): %d individuals\n", nrow(pheno_dt)))

# --- Write individual trait CSVs to phenotypes/ ---
cat("\n--- Writing output ---\n")

for (trait in names(icd_definitions)) {
  trait_file <- file.path(pheno_dir, sprintf("%s.csv", trait))
  trait_dt <- pheno_dt[, c("FID", "IID", trait), with = FALSE]
  fwrite(trait_dt, trait_file)
  cat(sprintf("  %s: %s\n", trait, trait_file))
}

# Write combined PLINK-format phenotype file
combined_file <- file.path(out_dir, "phenotypes_combined.txt")
fwrite(pheno_dt, combined_file, sep = "\t")
cat(sprintf("  Combined: %s\n", combined_file))

# Write covariate file
covar_file <- file.path(out_dir, "covariates.txt")
fwrite(covar_dt, covar_file, sep = "\t")
cat(sprintf("  Covariates: %s\n", covar_file))

# Write summary
summary_file <- file.path(out_dir, "phenotype_summary.txt")
fwrite(summary_dt, summary_file, sep = "\t")
cat(sprintf("  Summary: %s\n", summary_file))

cat("\n=== Phenotype extraction complete ===\n")
