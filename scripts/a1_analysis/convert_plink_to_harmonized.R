#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# convert_plink_to_harmonized.R
# ------------------------------------------------------------------
# Converts a single PLINK .assoc.gz + .frq file pair to the
# harmonized format expected by the a1_analysis pipeline.
#
# Output columns: VAR_ID, RSID, Effect_Allele, P_VALUE, BETA, SE, N, MAF, EAF
#
# VAR_ID is constructed with alphabetically sorted alleles to match
# the convention in harmonize_sumstats.py (lines 53-56).
#
# Usage:
#   Rscript convert_plink_to_harmonized.R \
#     --source-dir /path/to/plink/gwas/ \
#     --output-dir /path/to/output/ \
#     --field-id 21001 \
#     --trait-name BMI
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

# ---- CLI arguments ----
option_list <- list(
  make_option("--source-dir", type = "character",
              help = "Directory containing {field_id}.assoc.gz and {field_id}.frq files"),
  make_option("--output-dir", type = "character",
              help = "Output directory for harmonized files"),
  make_option("--field-id", type = "character",
              help = "UKB field ID (e.g., 21001)"),
  make_option("--trait-name", type = "character",
              help = "Human-readable trait name (e.g., BMI)")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required args
stopifnot(!is.null(opt$`source-dir`), !is.null(opt$`output-dir`),
          !is.null(opt$`field-id`), !is.null(opt$`trait-name`))

source_dir  <- opt$`source-dir`
output_dir  <- opt$`output-dir`
field_id    <- opt$`field-id`
trait_name  <- opt$`trait-name`

assoc_file <- file.path(source_dir, paste0(field_id, ".assoc.gz"))
frq_file   <- file.path(source_dir, paste0(field_id, ".frq"))

stopifnot(file.exists(assoc_file), file.exists(frq_file))

cat(sprintf("Processing field_id=%s (%s)\n", field_id, trait_name))
cat(sprintf("  assoc: %s\n  frq:   %s\n", assoc_file, frq_file))

# ---- Read PLINK association results ----
assoc <- fread(assoc_file,
               select = c("CHR", "SNP", "BP", "NMISS", "BETA", "SE", "P", "A1", "A2"))
cat(sprintf("  Rows in .assoc: %d\n", nrow(assoc)))

# Remove rows with missing P or BETA (PLINK outputs NA for monomorphic/failed)
assoc <- assoc[!is.na(P) & !is.na(BETA)]
cat(sprintf("  Rows after removing NA P/BETA: %d\n", nrow(assoc)))

# ---- Read allele frequency file ----
frq <- fread(frq_file, select = c("SNP", "A1", "A2", "MAF"))
setnames(frq, c("A1", "A2"), c("FRQ_A1", "FRQ_A2"))

# ---- Join MAF onto association results ----
assoc <- merge(assoc, frq, by = "SNP", all.x = TRUE)
cat(sprintf("  Rows after MAF join: %d\n", nrow(assoc)))

# ---- Compute EAF ----
# In PLINK --linear, A1 is the tested (effect) allele.
# In .frq, MAF is reported for FRQ_A1.
# If assoc A1 == frq A1, then EAF = MAF.
# If assoc A1 == frq A2, then EAF = 1 - MAF.
assoc[, EAF := fifelse(A1 == FRQ_A1, MAF, 1 - MAF)]

# ---- Construct VAR_ID with alphabetically sorted alleles ----
# Matches harmonize_sumstats.py convention: sorted alleles for cross-file matching
assoc[, a1_sorted := fifelse(A1 <= A2, A1, A2)]
assoc[, a2_sorted := fifelse(A1 <= A2, A2, A1)]
assoc[, VAR_ID := paste(CHR, BP, a1_sorted, a2_sorted, sep = "_")]

# ---- Build output table ----
out <- assoc[, .(
  VAR_ID       = VAR_ID,
  RSID         = SNP,
  Effect_Allele = A1,
  P_VALUE      = P,
  BETA         = BETA,
  SE           = SE,
  N            = NMISS,
  MAF          = MAF,
  EAF          = EAF
)]

# Remove duplicates on VAR_ID (keep lowest P)
setorder(out, P_VALUE)
out <- out[!duplicated(VAR_ID)]
cat(sprintf("  Rows after dedup on VAR_ID: %d\n", nrow(out)))

# ---- Write output ----
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(output_dir,
                      sprintf("splitUKB.%s.%s.EUR.GRCh37.processed.txt.gz",
                              field_id, trait_name))
fwrite(out, out_file, sep = "\t", compress = "gzip")
cat(sprintf("  Written: %s (%d variants)\n", out_file, nrow(out)))
cat("Done.\n")
