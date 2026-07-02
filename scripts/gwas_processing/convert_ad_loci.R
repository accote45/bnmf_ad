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
# VAR_ID = CHR_POS_A1_A2 with A1/A2 sorted ALPHABETICALLY (a1=min, a2=max);
# Effect_Allele stored separately (the Z-matrix build flips BETA when effect
# alleles disagree). This guarantees the AD leads join to the trait files.
#
# Alleles + effect direction come from the explicit `effect_allele`/`other_allele`
# columns (authoritative). `LeadSNP` (chr:pos:a:a) is used only for CHR:POS and as
# a cross-check, because its allele ORDER is not reliably effect:other.
#
# Usage (on Minerva, from bnmf_ad project root):
#   module load R/4.2.0
#   Rscript scripts/gwas_processing/convert_ad_loci.R

suppressWarnings(suppressMessages({
  library(data.table)
}))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  if (flag %in% args) args[which(args == flag) + 1] else default
}

input_path  <- get_arg("--input",
  "/sc/arion/projects/paul_oreilly/data/GWASs/NonBiobanks/raw_data/ad/PGC3_Unpublished/uffleman_top_loci_ad_gwas.xlsx")
output_path <- get_arg("--output", "sumstats/harmonized/AD_PGC3_top127.META.GRCh37.processed.txt.gz")
sheet       <- get_arg("--sheet", 1)

col_leadsnp <- get_arg("--leadsnp-col", "LeadSNP")            # chr:pos[:a:a] -> CHR/POS
col_ea      <- get_arg("--ea-col",      "effect_allele")      # authoritative effect allele
col_oa      <- get_arg("--oa-col",      "other_allele")       # authoritative other allele
col_beta    <- get_arg("--beta-col",    "beta")               # w.r.t. effect_allele
col_se      <- get_arg("--se-col",      "standard_error")
col_p       <- get_arg("--p-col",       "p_value")            # preferred P source
col_z       <- get_arg("--z-col",       "z_value")            # fallback P source
col_n       <- get_arg("--n-col",       "neff")
col_rsid    <- get_arg("--rsid-col",    "rsid")
col_eaf     <- get_arg("--eaf-col",     "effect_allele_frequency")

if (!file.exists(input_path)) stop(sprintf("Input not found: %s", input_path))

read_any <- function(path, sheet = 1) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (requireNamespace("readxl", quietly = TRUE))
      return(as.data.table(readxl::read_excel(path, sheet = sheet)))
    if (requireNamespace("openxlsx", quietly = TRUE))
      return(as.data.table(openxlsx::read.xlsx(path, sheet = as.integer(sheet))))
    stop("Reading .xlsx needs the 'readxl' or 'openxlsx' R package. ",
         "Install one (install.packages('readxl')), or re-save as .tsv and pass to --input.")
  }
  fread(path)
}

cat(sprintf("Reading: %s\n", input_path))
dt <- read_any(input_path, sheet)
setnames(dt, trimws(names(dt)))
cat(sprintf("  %d rows x %d cols. Columns: %s\n",
            nrow(dt), ncol(dt), paste(names(dt), collapse = ", ")))

# --- Resolve alleles: prefer explicit columns, fall back to LeadSNP fields 3/4 ---
lead  <- as.character(dt[[col_leadsnp]])
parts <- tstrsplit(lead, ":", fixed = TRUE)
if (length(parts) < 2)
  stop(sprintf("LeadSNP ('%s') did not split into >=2 fields on ':'. Example: %s",
               col_leadsnp, lead[1]))
CHR <- gsub("^chr", "", parts[[1]], ignore.case = TRUE)
POS <- suppressWarnings(as.integer(parts[[2]]))
lead_a3 <- if (length(parts) >= 3) toupper(trimws(parts[[3]])) else NA_character_
lead_a4 <- if (length(parts) >= 4) toupper(trimws(parts[[4]])) else NA_character_

use_explicit <- col_ea %in% names(dt) && col_oa %in% names(dt)
if (use_explicit) {
  EA <- toupper(trimws(as.character(dt[[col_ea]])))
  OA <- toupper(trimws(as.character(dt[[col_oa]])))
  cat(sprintf("  Alleles from explicit '%s'/'%s' columns (authoritative).\n", col_ea, col_oa))
} else {
  EA <- lead_a3; OA <- lead_a4
  cat("  WARNING: no explicit effect/other allele columns; parsed from LeadSNP fields 3/4.\n")
}

for (nm in c(col_beta, col_se)) if (!nm %in% names(dt)) stop(sprintf("Missing column: %s", nm))

out <- data.table(CHR = CHR, POS = POS, EA = EA, OA = OA)
out[, BETA := suppressWarnings(as.numeric(dt[[col_beta]]))]
out[, SE   := suppressWarnings(as.numeric(dt[[col_se]]))]

# P: prefer the reported p_value column, then z_value, then BETA/SE.
if (col_p %in% names(dt)) {
  out[, P_VALUE := suppressWarnings(as.numeric(dt[[col_p]]))]
  cat(sprintf("  P from '%s' column.\n", col_p))
} else if (col_z %in% names(dt)) {
  out[, P_VALUE := 2 * pnorm(-abs(suppressWarnings(as.numeric(dt[[col_z]]))))]
} else {
  out[, P_VALUE := 2 * pnorm(-abs(BETA / SE))]
}

out[, N := if (col_n %in% names(dt)) suppressWarnings(as.numeric(dt[[col_n]]))
          else if (all(c("n_case","n_control") %in% names(dt)))
            suppressWarnings(as.numeric(dt[["n_case"]]) + as.numeric(dt[["n_control"]]))
          else NA_real_]

out[, RSID := if (col_rsid %in% names(dt)) as.character(dt[[col_rsid]]) else NA_character_]

if (col_eaf %in% names(dt)) {
  out[, EAF := suppressWarnings(as.numeric(dt[[col_eaf]]))]
  out[, MAF := pmin(EAF, 1 - EAF)]
} else {
  out[, EAF := NA_real_]; out[, MAF := NA_real_]
}

# VAR_ID with alphabetically sorted alleles; Effect_Allele = EA (matches BETA sign)
out[, A1 := pmin(EA, OA)]
out[, A2 := pmax(EA, OA)]
out[, VAR_ID := sprintf("%s_%d_%s_%s", CHR, POS, A1, A2)]
out[, Effect_Allele := EA]

# --- QC / sanity ---
bad_pos    <- out[is.na(POS)]
bad_beta   <- out[is.na(BETA) | is.na(SE)]
bad_allele <- out[!(EA %in% c("A","C","G","T") & OA %in% c("A","C","G","T"))]
dups       <- out[duplicated(VAR_ID) | duplicated(VAR_ID, fromLast = TRUE)]

cat("\n--- QC ---\n")
cat(sprintf("  rows in: %d\n", nrow(out)))
if (nrow(bad_pos))    cat(sprintf("  WARNING: %d rows with unparseable POS\n", nrow(bad_pos)))
if (nrow(bad_beta))   cat(sprintf("  WARNING: %d rows missing BETA/SE\n", nrow(bad_beta)))
if (nrow(bad_allele)) cat(sprintf("  NOTE: %d non-SNP (indel/multichar) rows — SNP-only kept downstream\n", nrow(bad_allele)))
if (nrow(dups))       cat(sprintf("  WARNING: %d duplicated VAR_IDs\n", nrow(dups)))
maxP <- max(out$P_VALUE, na.rm = TRUE)
cat(sprintf("  P range: %.2e .. %.2e\n", min(out$P_VALUE, na.rm = TRUE), maxP))
if (maxP >= 5e-8)
  cat(sprintf("  NOTE: max P = %.2e is >= 5e-8 -> %d lead(s) would be dropped at p_threshold 5e-8.\n",
              maxP, sum(out$P_VALUE >= 5e-8, na.rm = TRUE)))

# Cross-check explicit alleles vs LeadSNP fields 3/4 (as an unordered set)
if (use_explicit && !all(is.na(lead_a3))) {
  set_ok <- (EA == lead_a3 & OA == lead_a4) | (EA == lead_a4 & OA == lead_a3)
  n_swap <- sum(EA != lead_a3 & set_ok, na.rm = TRUE)
  n_bad  <- sum(!set_ok, na.rm = TRUE)
  cat(sprintf("  Allele cross-check vs LeadSNP: %d order-swapped (harmless), %d true set-mismatches.\n",
              n_swap, n_bad))
  if (n_bad > 0)
    cat("  WARNING: set-mismatch rows have alleles that disagree with LeadSNP — inspect these.\n")
}

keep  <- !is.na(out$POS) & !is.na(out$BETA) & !is.na(out$SE)
final <- unique(out[keep, .(VAR_ID, RSID, Effect_Allele, P_VALUE, BETA, SE, N, MAF, EAF)],
                by = "VAR_ID")

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
fwrite(final, output_path, sep = "\t")
cat(sprintf("\nWrote %d variants -> %s\n", nrow(final), output_path))
cat("Preview:\n"); print(head(final, 3))
