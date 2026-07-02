#!/usr/bin/env Rscript
# 05b_find_overlapping_loci.R
# Identify overlapping GWAS loci between T2D and 5 CVD traits.
# Overlap defined as: lead SNPs within 500kb OR LD r2 > 0.2
#
# Inputs:  results/a0_analysis/loci_overlap/clumped/{trait}_clumped.txt
# Outputs: results/a0_analysis/loci_overlap/overlap_pairwise.csv
#          results/a0_analysis/loci_overlap/overlap_summary.csv
#          results/a0_analysis/loci_overlap/overlap_all5.csv

suppressPackageStartupMessages({
  library(data.table)
})

# --- Configuration ---
PROJECT_DIR <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
CLUMPED_DIR <- file.path(PROJECT_DIR, "results/a0_analysis/loci_overlap/clumped")
OUT_DIR     <- file.path(PROJECT_DIR, "results/a0_analysis/loci_overlap")
REF_PREFIX  <- file.path(PROJECT_DIR, "reference/1kg_eur/1000G.EUR.QC")

TRAITS <- c("t2d", "cad", "mi", "angina", "stroke", "pad")
CVD_TRAITS <- c("cad", "mi", "angina", "stroke", "pad")

DIST_THRESHOLD <- 500000   # 500kb
LD_R2_THRESHOLD <- 0.2
LD_WINDOW_KB <- 10000      # 10Mb window for plink --r2

# --- Step 1: Read clumped lead SNPs ---
cat("=== Loading clumped lead SNPs ===\n")

read_clumped <- function(trait) {
  f <- file.path(CLUMPED_DIR, paste0(trait, "_clumped.txt"))
  if (!file.exists(f)) {
    cat("  WARNING: No clumped file for", trait, "\n")
    return(data.table(SNP = character(), CHR = integer(), BP = integer(), P = numeric()))
  }
  # Plink .clumped is whitespace-delimited; SP2 column can wrap lines, so use fill=TRUE
  dt <- fread(f, header = TRUE, fill = TRUE)
  # Keep only relevant columns
  dt <- dt[, .(SNP, CHR, BP, P)]
  # Remove empty rows
  dt <- dt[!is.na(SNP) & SNP != ""]
  cat(sprintf("  %s: %d lead SNPs\n", trait, nrow(dt)))
  dt
}

leads <- lapply(setNames(TRAITS, TRAITS), read_clumped)

# --- Step 2: Distance-based overlap ---
cat("\n=== Distance-based overlap (within", DIST_THRESHOLD / 1000, "kb) ===\n")

t2d_leads <- leads[["t2d"]]

find_distance_overlap <- function(cvd_trait) {
  cvd_leads <- leads[[cvd_trait]]
  if (nrow(t2d_leads) == 0 || nrow(cvd_leads) == 0) {
    return(data.table())
  }

  # Non-equi join: same CHR, |BP_diff| <= threshold
  setnames(cvd_leads, c("SNP", "CHR", "BP", "P"),
           c("CVD_SNP", "CHR", "CVD_BP", "CVD_P"))

  cvd_leads[, `:=`(BP_lo = CVD_BP - DIST_THRESHOLD, BP_hi = CVD_BP + DIST_THRESHOLD)]

  overlap <- t2d_leads[cvd_leads,
    on = .(CHR = CHR, BP >= BP_lo, BP <= BP_hi),
    .(T2D_SNP = x.SNP, CHR = x.CHR, T2D_BP = x.BP,
      T2D_P = x.P, CVD_SNP = i.CVD_SNP, CVD_BP = i.CVD_BP, CVD_P = i.CVD_P),
    nomatch = NULL, allow.cartesian = TRUE
  ]

  # Clean up temp columns
  cvd_leads[, `:=`(BP_lo = NULL, BP_hi = NULL)]
  setnames(cvd_leads, c("CVD_SNP", "CHR", "CVD_BP", "CVD_P"),
           c("SNP", "CHR", "BP", "P"))

  if (nrow(overlap) > 0) {
    overlap[, distance_bp := abs(T2D_BP - CVD_BP)]
    overlap[, CVD_trait := cvd_trait]
    overlap[, overlap_type := "distance"]
  }

  overlap
}

dist_overlaps <- rbindlist(lapply(CVD_TRAITS, find_distance_overlap), fill = TRUE)
cat(sprintf("  Total distance-based overlaps: %d\n", nrow(dist_overlaps)))

# --- Step 3: LD-based overlap ---
cat("\n=== LD-based overlap (r2 >", LD_R2_THRESHOLD, ") ===\n")

# Collect all unique lead SNPs across all traits for LD computation
all_t2d_snps <- t2d_leads$SNP
all_cvd_snps <- unique(unlist(lapply(CVD_TRAITS, function(tr) leads[[tr]]$SNP)))

# Already-overlapping pairs (by distance) — skip these for LD
if (nrow(dist_overlaps) > 0) {
  dist_pairs <- dist_overlaps[, paste(T2D_SNP, CVD_SNP, sep = "_")]
} else {
  dist_pairs <- character(0)
}

# Compute LD per chromosome using plink
ld_results <- list()

for (chr in 1:22) {
  t2d_chr <- t2d_leads[CHR == chr, SNP]
  cvd_chr <- unique(unlist(lapply(CVD_TRAITS, function(tr) leads[[tr]][CHR == chr, SNP])))

  if (length(t2d_chr) == 0 || length(cvd_chr) == 0) next

  # Write SNP lists to temp files
  t2d_snp_file <- tempfile(pattern = "t2d_snps_", fileext = ".txt")
  cvd_snp_file <- tempfile(pattern = "cvd_snps_", fileext = ".txt")
  writeLines(t2d_chr, t2d_snp_file)
  writeLines(cvd_chr, cvd_snp_file)

  ld_out <- tempfile(pattern = paste0("ld_chr", chr, "_"))

  # Use plink to compute r2 between T2D leads and CVD leads
  cmd <- sprintf(
    "plink --bfile %s.%d --r2 --ld-snp-list %s --ld-window-kb %d --ld-window 99999 --ld-window-r2 %g --out %s 2>/dev/null",
    REF_PREFIX, chr, t2d_snp_file, LD_WINDOW_KB, LD_R2_THRESHOLD, ld_out
  )
  system(cmd)

  ld_file <- paste0(ld_out, ".ld")
  if (file.exists(ld_file)) {
    ld_dt <- fread(ld_file)
    if (nrow(ld_dt) > 0) {
      # Filter: SNP_A is T2D lead, SNP_B is CVD lead (or vice versa)
      ld_dt <- ld_dt[
        (SNP_A %in% t2d_chr & SNP_B %in% cvd_chr) |
        (SNP_B %in% t2d_chr & SNP_A %in% cvd_chr)
      ]
      # Normalize so T2D SNP is always first
      swap <- ld_dt$SNP_A %in% cvd_chr
      if (any(swap)) {
        ld_dt[swap, c("SNP_A", "SNP_B", "BP_A", "BP_B") :=
               .(SNP_B, SNP_A, BP_B, BP_A)]
      }
      ld_results[[length(ld_results) + 1]] <- ld_dt[, .(T2D_SNP = SNP_A, CVD_SNP = SNP_B, r2 = R2, CHR = CHR_A)]
    }
  }

  # Cleanup temp files
  unlink(c(t2d_snp_file, cvd_snp_file, paste0(ld_out, c(".ld", ".log", ".nosex"))))
}

if (length(ld_results) > 0) {
  ld_all <- rbindlist(ld_results)
  # Remove pairs already found by distance
  ld_all[, pair_key := paste(T2D_SNP, CVD_SNP, sep = "_")]
  ld_new <- ld_all[!pair_key %in% dist_pairs]
  ld_new[, pair_key := NULL]
  cat(sprintf("  LD-based pairs (not already distance-overlapping): %d\n", nrow(ld_new)))
} else {
  ld_new <- data.table(T2D_SNP = character(), CVD_SNP = character(), r2 = numeric(), CHR = integer())
  cat("  No LD-based overlaps found beyond distance overlaps.\n")
}

# Map LD overlaps back to traits
if (nrow(ld_new) > 0) {
  # For each CVD SNP, find which trait(s) it belongs to
  cvd_snp_to_trait <- rbindlist(lapply(CVD_TRAITS, function(tr) {
    data.table(CVD_SNP = leads[[tr]]$SNP, CVD_trait = tr,
               CVD_BP = leads[[tr]]$BP, CVD_P = leads[[tr]]$P)
  }))

  ld_overlap_rows <- merge(ld_new, cvd_snp_to_trait, by = "CVD_SNP", allow.cartesian = TRUE)

  # Add T2D info
  ld_overlap_rows <- merge(ld_overlap_rows,
    t2d_leads[, .(T2D_SNP = SNP, T2D_BP = BP, T2D_P = P)],
    by = "T2D_SNP"
  )
  ld_overlap_rows[, distance_bp := abs(T2D_BP - CVD_BP)]
  ld_overlap_rows[, overlap_type := "LD"]
} else {
  ld_overlap_rows <- data.table(
    T2D_SNP = character(), CHR = integer(), T2D_BP = integer(), T2D_P = numeric(),
    CVD_SNP = character(), CVD_BP = integer(), CVD_P = numeric(),
    CVD_trait = character(), distance_bp = integer(), overlap_type = character(), r2 = numeric()
  )
}

# --- Step 4: Merge all overlaps ---
cat("\n=== Merging overlaps ===\n")

# Add r2 = NA to distance overlaps
if (nrow(dist_overlaps) > 0) {
  dist_overlaps[, r2 := NA_real_]
}

all_overlaps <- rbindlist(list(dist_overlaps, ld_overlap_rows), fill = TRUE, use.names = TRUE)

# Reorder columns
cols <- c("T2D_SNP", "CHR", "T2D_BP", "T2D_P", "CVD_SNP", "CVD_trait", "CVD_BP", "CVD_P",
          "overlap_type", "distance_bp", "r2")
all_overlaps <- all_overlaps[, ..cols]

# Deduplicate (same T2D-CVD SNP pair for same trait)
all_overlaps <- unique(all_overlaps, by = c("T2D_SNP", "CVD_SNP", "CVD_trait"))

cat(sprintf("  Total unique overlapping pairs: %d\n", nrow(all_overlaps)))

# --- Step 5: Summary per CVD lead SNP ---
# Count from the CVD side to avoid inflating counts when many T2D lead SNPs
# in a dense region (e.g. TCF7L2) all match a single CVD lead SNP.
cat("\n=== Building summary tables (CVD-side counts) ===\n")

if (nrow(all_overlaps) > 0) {
  # Unique CVD lead SNPs that overlap at least one T2D locus, per trait
  cvd_summary <- unique(all_overlaps, by = c("CVD_SNP", "CVD_trait"))
  cvd_summary <- cvd_summary[, .(CVD_SNP, CHR, CVD_BP, CVD_P, CVD_trait,
                                  overlap_type, distance_bp, r2)]

  # Per-trait counts
  cvd_counts <- cvd_summary[, .(n_shared = .N), by = CVD_trait]

  # Per CVD lead SNP: which traits does it belong to, and does T2D overlap it?
  # (already guaranteed by being in all_overlaps)
  cvd_wide <- dcast(
    unique(all_overlaps[, .(CVD_SNP, CVD_trait, overlaps_T2D = TRUE)]),
    CVD_SNP ~ CVD_trait, value.var = "overlaps_T2D", fill = FALSE,
    fun.aggregate = function(x) TRUE
  )

  # T2D lead SNPs shared across all 5 CVD traits (T2D-side, for the all5 output)
  t2d_trait_counts <- unique(all_overlaps[, .(T2D_SNP, CVD_trait)])[
    , .(n_CVD_traits = .N), by = T2D_SNP
  ]
  all5_snps <- t2d_trait_counts[n_CVD_traits == 5, T2D_SNP]
  all5 <- t2d_leads[SNP %in% all5_snps, .(T2D_SNP = SNP, CHR, BP, T2D_P = P)]
  setorder(all5, CHR, BP)
} else {
  cvd_summary <- data.table(CVD_SNP = character(), CHR = integer(),
                            CVD_BP = integer(), CVD_P = numeric(),
                            CVD_trait = character(), overlap_type = character(),
                            distance_bp = integer(), r2 = numeric())
  cvd_counts <- data.table(CVD_trait = character(), n_shared = integer())
  cvd_wide <- data.table(CVD_SNP = character())
  all5 <- data.table(T2D_SNP = character(), CHR = integer(), BP = integer(),
                     T2D_P = numeric())
}

# --- Step 6: Print summary and write outputs ---
cat("\n=== Results Summary (CVD lead SNPs overlapping T2D) ===\n")
for (tr in CVD_TRAITS) {
  n <- cvd_counts[CVD_trait == tr, n_shared]
  if (length(n) == 0) n <- 0
  total <- nrow(leads[[tr]])
  cat(sprintf("  T2D - %s: %d / %d lead SNPs shared\n", toupper(tr), n, total))
}
cat(sprintf("  T2D loci shared with ALL 5 CVD: %d\n", nrow(all5)))

fwrite(all_overlaps, file.path(OUT_DIR, "overlap_pairwise.csv"))
fwrite(cvd_summary, file.path(OUT_DIR, "overlap_summary_cvd_side.csv"))
fwrite(all5, file.path(OUT_DIR, "overlap_all5.csv"))

cat(sprintf("\nOutputs written to %s/\n", OUT_DIR))
cat("  overlap_pairwise.csv\n")
cat("  overlap_summary_cvd_side.csv\n")
cat("  overlap_all5.csv\n")
cat("\n=== Done ===\n")
