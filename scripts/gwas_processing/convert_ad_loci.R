#!/usr/bin/env Rscript
# convert_ad_loci.R
# Convert the AD "top loci" table (127 genome-wide-significant lead SNPs from the
# full-cohort PGC3 AD analysis) into the pipeline's harmonized reference format.
#
# This is the REFERENCE GWAS for the AD bNMF run: it defines the variant universe
# (matrix rows). It is a small locus-annotation table, so we convert it directly
# instead of routing it through the genome-wide harmonize_sumstats.py.
#
# Output columns (match harmonize_sumstats.py exactly):
#   VAR_ID  RSID  Effect_Allele  P_VALUE  BETA  SE  N  MAF  EAF
# where VAR_ID = CHR_POS_A1_A2 with A1/A2 sorted ALPHABETICALLY (a1=min, a2=max),
# and Effect_Allele stored separately (the Z-matrix build flips BETA when effect
# alleles disagree). This guarantees the AD leads join to the trait files.
#
# Input: the LeadSNP column is coded  chr:position:effect_allele:other_allele  (hg19).
#
# Usage (on Minerva, from bnmf_ad project root):
#   module load R/4.2.0
#   Rscript scripts/gwas_processing/convert_ad_loci.R \
#     --input /sc/arion/projects/paul_oreilly/data/GWASs/NonBiobanks/raw_data/ad/PGC3_Unpublished/uffleman_top_loci_ad_gwas.xlsx \
#     --output sumstats/harmonized/AD_PGC3_top127.META.GRCh37.processed.txt.gz
#
# Column names default to the observed header but can be overridden via flags
# (--leadsnp-col, --beta-col, --se-col, --z-col, --n-col, --rsid-col, --sheet).

suppressWarnings(suppressMessages({
  library(data.table)
}))

# ---------------------------------------------------------------------------
# CLI args
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  if (flag %in% args) args[which(args == flag) + 1] else default
}

input_path  <- get_arg("--input",
  "/sc/arion/projects/paul_oreilly/data/GWASs/NonBiobanks/raw_data/ad/PGC3_Unpublished/uffleman_top_loci_ad_gwas.xlsx")
output_path <- get_arg("--output", "sumstats/harmonized/AD_PGC3_top127.META.GRCh37.processed.txt.gz")
sheet       <- get_arg("--sheet", 1)

col_leadsnp <- get_arg("--leadsnp-col", "LeadSNP")   # coded chr:pos:effect:other
col_beta    <- get_arg("--beta-col",    "beta")
col_se      <- get_arg("--se-col",      "standard_error")
col_z       <- get_arg("--z-col",       "z_value")   # optional; used for P if present
col_n       <- get_arg("--n-col",       "neff")      # else n_case + n_control
col_rsid    <- get_arg("--rsid-col",    NULL)        # optional rsID column, if any

if (!file.exists(input_path)) stop(sprintf("Input not found: %s", input_path))

# ---------------------------------------------------------------------------
# Read xlsx / csv / tsv
# ---------------------------------------------------------------------------
read_any <- function(path, sheet = 1) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (requireNamespace("readxl", quietly = TRUE)) {
      return(as.data.table(readxl::read_excel(path, sheet = sheet)))
    }
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      return(as.data.table(openxlsx::read.xlsx(path, sheet = as.integer(sheet))))
    }
    stop("Reading .xlsx needs the 'readxl' or 'openxlsx' R package. ",
         "Install one (install.packages('readxl')), or re-save the sheet as .tsv ",
         "and pass that to --input.")
  }
  fread(path)  # csv / tsv / txt
}

cat(sprintf("Reading: %s\n", input_path))
dt <- read_any(input_path, sheet)
setnames(dt, trimws(names(dt)))   # normalize header whitespace
cat(sprintf("  %d rows x %d cols. Columns: %s\n",
            nrow(dt), ncol(dt), paste(names(dt), collapse = ", ")))

need <- c(col_leadsnp, col_beta, col_se)
missing_cols <- setdiff(need, names(dt))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing expected column(s): %s\nOverride with --leadsnp-col/--beta-col/--se-col.",
               paste(missing_cols, collapse = ", ")))
}

# ---------------------------------------------------------------------------
# Parse LeadSNP = chr:pos:effect:other
# ---------------------------------------------------------------------------
lead <- as.character(dt[[col_leadsnp]])
parts <- tstrsplit(lead, ":", fixed = TRUE)
if (length(parts) < 4) {
  stop(sprintf("LeadSNP ('%s') did not split into 4 fields on ':'. Example value: %s",
               col_leadsnp, lead[1]))
}
CHR <- gsub("^chr", "", parts[[1]], ignore.case = TRUE)
POS <- suppressWarnings(as.integer(parts[[2]]))
EA  <- toupper(trimws(parts[[3]]))   # effect allele (field 3)
OA  <- toupper(trimws(parts[[4]]))   # other allele  (field 4)

out <- data.table(CHR = CHR, POS = POS, EA = EA, OA = OA)
out[, BETA := suppressWarnings(as.numeric(dt[[col_beta]]))]
out[, SE   := suppressWarnings(as.numeric(dt[[col_se]]))]

# P: prefer deriving from z_value if present, else from BETA/SE
if (!is.null(col_z) && col_z %in% names(dt)) {
  zc <- suppressWarnings(as.numeric(dt[[col_z]]))
  out[, P_VALUE := 2 * pnorm(-abs(zc))]
} else {
  out[, P_VALUE := 2 * pnorm(-abs(BETA / SE))]
}

# N: neff if present, else n_case + n_control if both present, else NA
if (!is.null(col_n) && col_n %in% names(dt)) {
  out[, N := suppressWarnings(as.numeric(dt[[col_n]]))]
} else if (all(c("n_case", "n_control") %in% names(dt))) {
  out[, N := suppressWarnings(as.numeric(dt[["n_case"]]) + as.numeric(dt[["n_control"]]))]
} else {
  out[, N := NA_real_]
}

# RSID: optional; clumping re-derives rsIDs from the reference panel by position
if (!is.null(col_rsid) && col_rsid %in% names(dt)) {
  out[, RSID := as.character(dt[[col_rsid]])]
} else {
  out[, RSID := NA_character_]
}

# MAF / EAF: not needed for the 127 GW-significant leads (MAF filter is optional)
out[, MAF := NA_real_]
out[, EAF := NA_real_]

# ---------------------------------------------------------------------------
# Build VAR_ID (alleles sorted alphabetically) + Effect_Allele
# ---------------------------------------------------------------------------
out[, A1 := pmin(EA, OA)]
out[, A2 := pmax(EA, OA)]
out[, VAR_ID := sprintf("%s_%d_%s_%s", CHR, POS, A1, A2)]
out[, Effect_Allele := EA]

# ---------------------------------------------------------------------------
# QC / sanity
# ---------------------------------------------------------------------------
bad_pos    <- out[is.na(POS)]
bad_beta   <- out[is.na(BETA) | is.na(SE)]
bad_allele <- out[!(EA %in% c("A","C","G","T") & OA %in% c("A","C","G","T"))]
dups       <- out[duplicated(VAR_ID) | duplicated(VAR_ID, fromLast = TRUE)]

cat("\n--- QC ---\n")
cat(sprintf("  rows in: %d\n", nrow(out)))
if (nrow(bad_pos))    cat(sprintf("  WARNING: %d rows with unparseable POS\n", nrow(bad_pos)))
if (nrow(bad_beta))   cat(sprintf("  WARNING: %d rows missing BETA/SE\n", nrow(bad_beta)))
if (nrow(bad_allele)) cat(sprintf("  NOTE: %d rows are non-SNP (indel/multichar alleles) — QC keeps SNPs only downstream\n", nrow(bad_allele)))
if (nrow(dups))       cat(sprintf("  WARNING: %d duplicated VAR_IDs\n", nrow(dups)))
cat(sprintf("  P range: %.2e .. %.2e (all should be < 5e-8)\n",
            min(out$P_VALUE, na.rm = TRUE), max(out$P_VALUE, na.rm = TRUE)))

# Cross-check LeadSNP effect allele vs an explicit effect_allele column, if present
if ("effect_allele" %in% names(dt)) {
  mism <- sum(toupper(trimws(as.character(dt[["effect_allele"]]))) != out$EA, na.rm = TRUE)
  if (mism > 0) cat(sprintf("  WARNING: %d rows where LeadSNP effect allele != 'effect_allele' column\n", mism))
}

# Drop unusable rows
keep <- !is.na(out$POS) & !is.na(out$BETA) & !is.na(out$SE)
out  <- out[keep]

final <- out[, .(VAR_ID, RSID, Effect_Allele, P_VALUE, BETA, SE, N, MAF, EAF)]
final <- unique(final, by = "VAR_ID")

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
fwrite(final, output_path, sep = "\t")
cat(sprintf("\nWrote %d variants -> %s\n", nrow(final), output_path))
cat("Preview:\n")
print(head(final, 3))
