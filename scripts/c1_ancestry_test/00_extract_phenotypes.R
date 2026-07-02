#!/usr/bin/env Rscript
# 00_extract_phenotypes.R
# Extract UKB phenotypes and covariates from the SQLite database for the
# C1.1 pipeline. Writes PLINK-format pheno and covar files.
#
# Phenotype fields (baseline, instance 0):
#   LDL:     f30780 (blood biochemistry LDL direct, mmol/L)
#   glucose: f30740 (blood biochemistry glucose, mmol/L)
#   BMI:     f21001 (body mass index, kg/m^2)
#
# Covariate fields:
#   age:  f21022 (age at recruitment)
#   age2: age^2
#   sex:  f31    (1=Male, 2=Female)
#   PCs:  from pre-computed covariate file (PC1-PC10)
#
# Pattern based on:
#   /sc/arion/projects/psychgen/projects/prs/cad_subtype/scripts/phenotype/wrangle_sql_statin_response.r
#
# Usage:
#   Rscript scripts/c1_ancestry_test/00_extract_phenotypes.R
#   Rscript scripts/c1_ancestry_test/00_extract_phenotypes.R \
#     --db /path/to/ukb18177-all.db \
#     --covar /path/to/ukb18177.covar \
#     --outdir results/c1_ancestry_test

library(RSQLite)
library(data.table)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)

db_path    <- "/sc/arion/projects/ukb18177_oreilly/ukb/phenotype/ukb18177-all.db"
covar_path <- "/sc/arion/projects/ukb18177_oreilly/ukb/phenotype/ukb18177.covar"
out_dir    <- "results/c1_ancestry_test"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--db") {
    db_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--covar") {
    covar_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--outdir") {
    out_dir <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== UKB Phenotype Extraction ===\n")
cat(sprintf("  DB: %s\n", db_path))
cat(sprintf("  Covar file: %s\n", covar_path))
cat(sprintf("  Output dir: %s\n\n", out_dir))

# --- Connect to UKB database ---
con <- dbConnect(SQLite(), dbname = db_path)
on.exit(dbDisconnect(con), add = TRUE)

# Helper: pull baseline (instance 0) field
get_baseline <- function(field_id, col_name) {
  tbl_name <- paste0("f", field_id)
  if (!dbExistsTable(con, tbl_name)) {
    stop(sprintf("Table '%s' not found in database", tbl_name))
  }
  query <- sprintf(
    "SELECT sample_id, pheno AS %s FROM %s WHERE instance = 0",
    col_name, tbl_name
  )
  as.data.table(dbGetQuery(con, query))
}

# --- Get non-withdrawn participant IDs ---
cat("--- Extracting participant IDs ---\n")
if (dbExistsTable(con, "Participant")) {
  ids <- as.data.table(dbGetQuery(con,
    "SELECT sample_id FROM Participant WHERE withdrawn = 0"))
  cat(sprintf("  Non-withdrawn participants: %d\n", nrow(ids)))
} else {
  # Fallback: get all unique sample_ids from a common field
  ids <- as.data.table(dbGetQuery(con,
    "SELECT DISTINCT sample_id FROM f31"))
  cat(sprintf("  Participants (from f31): %d\n", nrow(ids)))
}

# --- Extract phenotypes ---
cat("\n--- Extracting phenotypes ---\n")

ldl <- get_baseline(30780, "LDL")
cat(sprintf("  LDL (f30780): %d non-NA values out of %d\n",
            sum(!is.na(ldl$LDL)), nrow(ldl)))

glucose <- get_baseline(30740, "glucose")
cat(sprintf("  Glucose (f30740): %d non-NA values out of %d\n",
            sum(!is.na(glucose$glucose)), nrow(glucose)))

bmi <- get_baseline(21001, "BMI")
cat(sprintf("  BMI (f21001): %d non-NA values out of %d\n",
            sum(!is.na(bmi$BMI)), nrow(bmi)))

# --- Extract covariates ---
cat("\n--- Extracting covariates ---\n")

age <- get_baseline(21022, "age")
cat(sprintf("  Age (f21022): %d non-NA values\n", sum(!is.na(age$age))))

sex <- get_baseline(31, "sex")
cat(sprintf("  Sex (f31): %d non-NA values\n", sum(!is.na(sex$sex))))

# --- Merge phenotypes ---
cat("\n--- Merging ---\n")

pheno <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = FALSE),
                list(ids, ldl, glucose, bmi))
pheno[, `:=`(FID = sample_id, IID = sample_id)]
pheno <- pheno[, .(FID, IID, LDL, glucose, BMI)]

cat(sprintf("  Phenotype table: %d individuals with all 3 traits\n", nrow(pheno)))
cat(sprintf("  LDL non-NA: %d\n", sum(!is.na(pheno$LDL))))
cat(sprintf("  Glucose non-NA: %d\n", sum(!is.na(pheno$glucose))))
cat(sprintf("  BMI non-NA: %d\n", sum(!is.na(pheno$BMI))))

# --- Merge covariates ---
covar_demo <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = FALSE),
                     list(ids, age, sex))

# Read pre-computed PCs
cat("\n--- Reading pre-computed PCs ---\n")
pcs <- fread(covar_path)
cat(sprintf("  PC file: %d rows, columns: %s\n",
            nrow(pcs), paste(head(colnames(pcs), 15), collapse = ", ")))

# Merge demographics with PCs
covar <- merge(covar_demo, pcs[, .(FID, IID, PC1, PC2, PC3, PC4, PC5,
                                    PC6, PC7, PC8, PC9, PC10)],
               by.x = "sample_id", by.y = "IID")
covar[, age := as.numeric(age)]
covar[, `:=`(FID = sample_id, IID = sample_id, age2 = age^2)]
covar <- covar[, .(FID, IID, age, age2, sex, PC1, PC2, PC3, PC4, PC5,
                    PC6, PC7, PC8, PC9, PC10)]

cat(sprintf("  Covariate table: %d individuals\n", nrow(covar)))

# --- Write output ---
cat("\n--- Writing output ---\n")

pheno_file <- file.path(out_dir, "phenotypes.txt")
fwrite(pheno, pheno_file, sep = "\t")
cat(sprintf("  Phenotypes: %s (%d rows)\n", pheno_file, nrow(pheno)))

covar_file <- file.path(out_dir, "covariates.txt")
fwrite(covar, covar_file, sep = "\t")
cat(sprintf("  Covariates: %s (%d rows)\n", covar_file, nrow(covar)))

cat("\n=== Phenotype extraction complete ===\n")
