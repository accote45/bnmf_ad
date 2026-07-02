#!/usr/bin/env Rscript
# 01_prepare_metal_input.R
# Prepare T2D and CAD EUR GWAS files for METAL IVW meta-analysis.
# Also build RSID lookup from UKB .bim files (T2D GWAS has "." for RSIDs).
#
# Usage:
#   Rscript scripts/b1_2_analysis/01_prepare_metal_input.R
#   Rscript scripts/b1_2_analysis/01_prepare_metal_input.R --config config/b1_2_config.yaml

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
dir.create(metal_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== B1.2 Step 1: Prepare METAL Input ===\n")

# --- Load GWAS files ---
cat("\n--- Loading GWAS files ---\n")
t2d <- fread(cfg$gwas$t2d, select = c("VAR_ID", "RSID", "Effect_Allele", "BETA", "SE", "P_VALUE"))
cat(sprintf("  T2D: %d variants\n", nrow(t2d)))

cad <- fread(cfg$gwas$cad, select = c("VAR_ID", "RSID", "Effect_Allele", "BETA", "SE", "P_VALUE"))
cat(sprintf("  CAD: %d variants\n", nrow(cad)))

# --- Derive Other_Allele from VAR_ID (CHR_POS_REF_ALT) ---
cat("\n--- Deriving Other_Allele from VAR_ID ---\n")

derive_other_allele <- function(dt) {
  dt[, c("v_chr", "v_pos", "v_ref", "v_alt") := tstrsplit(VAR_ID, "_", keep = 1:4)]
  dt[, Other_Allele := fifelse(Effect_Allele == v_ref, v_alt, v_ref)]
  dt[, c("v_chr", "v_pos", "v_ref", "v_alt") := NULL]
  dt
}

t2d <- derive_other_allele(t2d)
cad <- derive_other_allele(cad)

# --- Filter to valid variants ---
cat("\n--- Filtering variants ---\n")
t2d <- t2d[!is.na(BETA) & !is.na(SE) & SE > 0]
cad <- cad[!is.na(BETA) & !is.na(SE) & SE > 0]
cat(sprintf("  T2D after filter: %d variants\n", nrow(t2d)))
cat(sprintf("  CAD after filter: %d variants\n", nrow(cad)))

# --- Write METAL input files ---
# METAL format: MarkerName Allele1 Allele2 Effect StdErr
cat("\n--- Writing METAL input files ---\n")

write_metal_input <- function(dt, outpath, label) {
  metal_dt <- data.table(
    MarkerName = dt$VAR_ID,
    Allele1 = toupper(dt$Effect_Allele),
    Allele2 = toupper(dt$Other_Allele),
    Effect = dt$BETA,
    StdErr = dt$SE
  )
  fwrite(metal_dt, outpath, sep = "\t")
  cat(sprintf("  %s: %d variants -> %s\n", label, nrow(metal_dt), outpath))
}

write_metal_input(t2d, file.path(metal_dir, "metal_input_t2d.txt"), "T2D")
write_metal_input(cad, file.path(metal_dir, "metal_input_cad.txt"), "CAD")

# --- Build RSID lookup from .bim files ---
cat("\n--- Building RSID lookup from .bim files ---\n")

geno_prefix <- cfg$genotypes$prefix
rsid_list <- list()

for (chr in 1:22) {
  bim_path <- paste0(gsub("\\{chr\\}", chr, geno_prefix), ".bim")
  if (!file.exists(bim_path)) {
    bim_path <- sprintf("%s/chr%d.bim", dirname(gsub("\\{chr\\}", "1", geno_prefix)), chr)
  }
  bim <- fread(bim_path, col.names = c("CHR", "SNP", "CM", "POS", "A1", "A2"))
  # Build VAR_ID from bim: CHR_POS_A1_A2 and CHR_POS_A2_A1 (both orientations)
  bim[, VAR_ID_1 := paste(CHR, POS, A1, A2, sep = "_")]
  bim[, VAR_ID_2 := paste(CHR, POS, A2, A1, sep = "_")]

  rsid_1 <- bim[grepl("^rs", SNP), .(VAR_ID = VAR_ID_1, RSID_bim = SNP)]
  rsid_2 <- bim[grepl("^rs", SNP), .(VAR_ID = VAR_ID_2, RSID_bim = SNP)]
  rsid_list[[chr]] <- rbindlist(list(rsid_1, rsid_2))

  cat(sprintf("  chr%d: %d RSIDs\n", chr, nrow(rsid_1)))
}

rsid_lookup <- rbindlist(rsid_list)
rsid_lookup <- unique(rsid_lookup, by = "VAR_ID")
fwrite(rsid_lookup, file.path(metal_dir, "rsid_lookup.tsv"), sep = "\t")
cat(sprintf("\n  RSID lookup: %d unique VAR_ID -> RSID mappings\n", nrow(rsid_lookup)))

# --- Write METAL command script ---
cat("\n--- Writing METAL command script ---\n")

metal_script <- '
SCHEME STDERR

MARKER MarkerName
ALLELE Allele1 Allele2
EFFECT Effect
STDERR StdErr

PROCESS metal_input_t2d.txt
PROCESS metal_input_cad.txt

OUTFILE meta_t2d_cad .tbl

ANALYZE HETEROGENEITY

QUIT
'

writeLines(metal_script, file.path(metal_dir, "run_metal.sh"))
cat(sprintf("  Written: %s\n", file.path(metal_dir, "run_metal.sh")))

cat("\n=== Step 1 complete ===\n")
