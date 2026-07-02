#!/usr/bin/env Rscript
# 03_prepare_prset_input.R
# Format METAL meta-analysis output for PRSice/PRSet input.
# Also combine EUR train + validation into full EUR sample list
# and create phenotype file with Synthetic (T2D | CAD) column.
#
# Usage:
#   Rscript scripts/b1_2_analysis/03_prepare_prset_input.R
#   Rscript scripts/b1_2_analysis/03_prepare_prset_input.R --config config/b1_2_config.yaml

library(data.table)
library(yaml)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b1_2_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir
metal_dir <- file.path(results_dir, "metal")
prset_dir <- file.path(results_dir, "prset")
dir.create(prset_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== B1.2 Step 3: Prepare PRSet Input ===\n")

# --- Load METAL output ---
cat("\n--- Loading METAL output ---\n")
meta <- fread(file.path(metal_dir, "meta_t2d_cad.txt"))
cat(sprintf("  METAL output: %d variants\n", nrow(meta)))
cat(sprintf("  Columns: %s\n", paste(colnames(meta), collapse = ", ")))

# --- Load RSID lookup ---
cat("\n--- Loading RSID lookup ---\n")
rsid_lookup <- fread(file.path(metal_dir, "rsid_lookup.tsv"))
cat(sprintf("  RSID lookup: %d mappings\n", nrow(rsid_lookup)))

# --- Parse METAL columns ---
# METAL output columns: MarkerName Allele1 Allele2 Effect StdErr P-value Direction HetISq HetChiSq HetDf HetPVal
# MarkerName = VAR_ID, Allele1/2 are lowercased by METAL

# Derive CHR and BP from MarkerName (VAR_ID format: CHR_POS_REF_ALT)
meta[, c("CHR", "BP") := tstrsplit(MarkerName, "_", keep = 1:2)]
meta[, CHR := as.integer(CHR)]
meta[, BP := as.integer(BP)]

# Uppercase alleles (METAL lowercases them)
meta[, A1 := toupper(Allele1)]
meta[, A2 := toupper(Allele2)]

# Rename columns for PRSice
setnames(meta, c("Effect", "P-value"), c("BETA", "P"), skip_absent = TRUE)

# --- Map VAR_ID -> RSID ---
cat("\n--- Mapping VAR_ID to RSID ---\n")
meta <- merge(meta, rsid_lookup, by.x = "MarkerName", by.y = "VAR_ID", all.x = TRUE)

n_with_rsid <- sum(!is.na(meta$RSID_bim) & grepl("^rs", meta$RSID_bim))
n_without <- sum(is.na(meta$RSID_bim) | !grepl("^rs", meta$RSID_bim))
cat(sprintf("  With RSID: %d\n", n_with_rsid))
cat(sprintf("  Without RSID (will be dropped): %d\n", n_without))

# Filter to variants with valid RSIDs
meta <- meta[!is.na(RSID_bim) & grepl("^rs", RSID_bim)]
meta[, SNP := RSID_bim]

# Remove duplicates by SNP (keep first)
n_before <- nrow(meta)
meta <- meta[!duplicated(SNP)]
cat(sprintf("  After dedup by RSID: %d (removed %d duplicates)\n", nrow(meta), n_before - nrow(meta)))

# --- Filter out problematic variants ---
meta <- meta[!is.na(BETA) & !is.na(P) & P > 0 & P <= 1]
cat(sprintf("  After filtering invalid BETA/P: %d variants\n", nrow(meta)))

# --- Write PRSice-format base file ---
cat("\n--- Writing PRSice base file ---\n")
prsice_dt <- meta[, .(SNP, CHR, BP, A1, A2, P, BETA)]
out_path <- file.path(prset_dir, "meta_t2d_cad_prsice.txt")
fwrite(prsice_dt, out_path, sep = "\t")
cat(sprintf("  Written: %s (%d variants)\n", out_path, nrow(prsice_dt)))

# --- Combine EUR train + validation into full EUR ---
cat("\n--- Combining EUR samples ---\n")
eur_train <- fread(cfg$samples$eur_train, header = FALSE)
eur_val <- fread(cfg$samples$eur_validation, header = FALSE)
eur_all <- rbindlist(list(eur_train, eur_val))
out_keep <- file.path(prset_dir, "eur_all.keep")
fwrite(eur_all, out_keep, sep = "\t", col.names = FALSE)
cat(sprintf("  EUR train: %d, validation: %d, total: %d\n",
            nrow(eur_train), nrow(eur_val), nrow(eur_all)))
cat(sprintf("  Written: %s\n", out_keep))

# --- Create phenotype file with Synthetic column ---
cat("\n--- Creating phenotype file with Synthetic column ---\n")
pheno <- fread(cfg$phenotypes$phenotype_file)
cat(sprintf("  Phenotype file: %d individuals\n", nrow(pheno)))

# Create Synthetic (T2D | CAD)
pheno[, Synthetic := fifelse(T2D == 1 | CAD == 1, 1L,
                     fifelse(T2D == 0 & CAD == 0, 0L, NA_integer_))]

cat(sprintf("  Synthetic: %d cases, %d controls, %d NA\n",
            sum(pheno$Synthetic == 1, na.rm = TRUE),
            sum(pheno$Synthetic == 0, na.rm = TRUE),
            sum(is.na(pheno$Synthetic))))

# Write phenotype file (FID, IID, T2D, CAD, Synthetic)
pheno_out <- pheno[, .(FID, IID, T2D, CAD, Synthetic)]
out_pheno <- file.path(prset_dir, "phenotypes_with_synthetic.txt")
fwrite(pheno_out, out_pheno, sep = "\t")
cat(sprintf("  Written: %s\n", out_pheno))

cat("\n=== Step 3 complete ===\n")
