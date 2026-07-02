#!/usr/bin/env Rscript
# 05c_write_loci_overlap_xlsx.R
# Build the supplementary workbook t2d_cvd_loci_overlap.xlsx summarising how
# T2D lead SNPs overlap each CVD trait's lead SNPs.
#
# Definitions (self-contained — computed directly from the 05a clumped lead
# files, NOT from 05b's pairwise output):
#   Lead SNPs        — plink --clump index SNPs (P<5e-8, r2>0.05 absorbed, 500kb;
#                      MHC excluded). One row per independent association signal.
#   Shared lead SNPs — T2D lead SNPs within SHARE_KB of any lead SNP of the CVD
#                      trait (same chromosome, |BP difference| <= SHARE_KB).
#   Shared loci      — those qualifying T2D lead SNPs collapsed into loci by
#                      merging SNPs within LOCUS_KB of each other (same chrom).
#
# Summary columns (one row per CVD trait):
#   n_T2D_lead_SNPs | CVD_trait | n_CVD_lead_SNPs | n_shared_lead_SNPs | n_shared_loci
#
# Inputs:  results/a0_analysis/loci_overlap/clumped/{trait}_clumped.txt
# Output:  results/supplementary_tables/t2d_cvd_loci_overlap.xlsx
#
# Usage:
#   Rscript scripts/a0_analysis/05c_write_loci_overlap_xlsx.R
#   Rscript scripts/a0_analysis/05c_write_loci_overlap_xlsx.R \
#     --overlap-dir results/a0_analysis/loci_overlap \
#     --out results/supplementary_tables/t2d_cvd_loci_overlap.xlsx \
#     --share-kb 100 --locus-kb 500

suppressMessages({
  library(data.table)
  library(openxlsx)
})

# --- Args ---
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1L]]
}
overlap_dir <- get_arg("--overlap-dir", "results/a0_analysis/loci_overlap")
out_file    <- get_arg("--out", "results/supplementary_tables/t2d_cvd_loci_overlap.xlsx")
SHARE_KB    <- as.numeric(get_arg("--share-kb", "100"))   # T2D-CVD lead proximity for "shared"
LOCUS_KB    <- as.numeric(get_arg("--locus-kb", "500"))   # window to merge shared leads into loci
SHARE_BP    <- SHARE_KB * 1000
LOCUS_BP    <- LOCUS_KB * 1000

CVD_TRAITS  <- c("cad", "mi", "angina", "stroke", "pad")
clumped_dir <- file.path(overlap_dir, "clumped")

# --- Read lead SNPs from a *_clumped.txt (plink --clump output) ---
read_leads <- function(trait) {
  f <- file.path(clumped_dir, paste0(trait, "_clumped.txt"))
  if (!file.exists(f)) stop("Missing clumped file: ", f)
  dt <- fread(f, header = TRUE, fill = TRUE)
  dt <- dt[!is.na(SNP) & SNP != "", .(SNP, CHR = as.integer(CHR), BP = as.integer(BP),
                                       P = as.numeric(P))]
  setorder(dt, CHR, BP)
  dt
}

# --- Merge lead SNPs into loci: new locus when chrom changes or gap > window ---
count_loci <- function(dt, window_bp) {
  if (nrow(dt) == 0) return(0L)
  setorder(dt, CHR, BP)
  newloc <- c(TRUE, dt$CHR[-1] != dt$CHR[-nrow(dt)] |
                    (dt$BP[-1] - dt$BP[-nrow(dt)]) > window_bp)
  sum(newloc)
}

# --- Assign lead SNPs to loci (integer id); same merge rule as count_loci ---
assign_loci <- function(dt, window_bp) {
  setorder(dt, CHR, BP)
  if (nrow(dt) == 0) { dt[, locus := integer()]; return(dt) }
  dt[, locus := cumsum(c(TRUE, CHR[-1] != CHR[-.N] | (BP[-1] - BP[-.N]) > window_bp))]
  dt
}

# --- Nearest protein-coding gene (GENCODE v19 / GRCh37) ---
gtf_path <- get_arg("--gtf", "reference/gencode.v19.annotation.gtf.gz")
message("Parsing protein-coding genes from ", gtf_path, " ...")
g <- fread(cmd = sprintf("zcat %s | awk -F'\\t' '$3==\"gene\"'", gtf_path),
           header = FALSE, sep = "\t", quote = "")
g <- g[grepl('gene_type "protein_coding"', V9)]
g[, `:=`(chr = sub("^chr", "", V1), start = as.integer(V4), end = as.integer(V5),
         gene = sub('.*gene_name "([^"]+)".*', "\\1", V9))]
gene_dt <- g[chr %in% as.character(1:22), .(chr, start, end, gene)]
# Drop uninformative clone/uncharacterized "protein_coding" symbols so labels are
# recognisable gene names. Catches dot-versioned clone IDs (e.g. AC022431.2,
# RP11-145E5.5) and the common GENCODE clone-library prefixes.
clone_re <- paste0("[.][0-9]+$|",
  "^(RP[0-9]+-|CT[ABCD]-|AC[0-9]{6}|AL[0-9]{6}|AP[0-9]{6}|GS1-|KB-|KC[0-9]{6}|",
  "WI2-|XX|LL[0-9]{5}|CH[0-9]{2,}-|Z[0-9]{5}|BX[0-9]{6}|FP[0-9]{6}|LINC[0-9])")
gene_dt <- gene_dt[!grepl(clone_re, gene)]
message(sprintf("  %d autosomal protein-coding genes (clone-named symbols excluded)", nrow(gene_dt)))

nearest_gene <- function(chr_i, pos_i) {
  gg <- gene_dt[chr == as.character(chr_i)]
  if (nrow(gg) == 0) return(NA_character_)
  d <- pmax(gg$start - pos_i, pos_i - gg$end, 0)  # 0 if inside gene body
  gg$gene[which.min(d)]
}

t2d   <- read_leads("t2d")
leads <- lapply(setNames(CVD_TRAITS, CVD_TRAITS), read_leads)
n_t2d_loci <- count_loci(t2d, LOCUS_BP)

# --- Per-CVD-trait: shared pairs -> loci -> joint strength + nearest gene ---
# For each CVD trait, find all (T2D lead, CVD lead) pairs within SHARE_BP, group
# the shared T2D leads into loci (merge within LOCUS_BP), label each locus by the
# nearest gene to its top (lowest-T2D-P) lead SNP, and rank loci by joint signal
# strength = min over pairs of max(T2D_P, CVD_P) (smaller = both more significant).
per_trait <- lapply(setNames(CVD_TRAITS, CVD_TRAITS), function(tr) {
  cvd <- leads[[tr]]
  jn <- t2d[cvd, on = .(CHR), allow.cartesian = TRUE,
            .(T2D_SNP = x.SNP, CHR, T2D_BP = x.BP, T2D_P = x.P,
              CVD_SNP = i.SNP, CVD_BP = i.BP, CVD_P = i.P)]
  jn <- jn[!is.na(T2D_BP) & abs(T2D_BP - CVD_BP) <= SHARE_BP]
  if (nrow(jn) == 0)
    return(list(n_shared_snps = 0L, n_loci = 0L, strongest = NA_character_,
                loci_list = "", detail = NULL))
  jn[, joint_p := pmax(T2D_P, CVD_P)]
  jn[, distance_bp := abs(T2D_BP - CVD_BP)]

  # Locus id from the unique shared T2D leads, mapped back onto the pairs
  shared <- assign_loci(unique(jn[, .(SNP = T2D_SNP, CHR, BP = T2D_BP)]), LOCUS_BP)
  jn <- merge(jn, shared, by.x = c("T2D_SNP", "CHR", "T2D_BP"),
              by.y = c("SNP", "CHR", "BP"), all.x = TRUE)

  # Per-locus: index = lowest-T2D-P shared lead; joint = min max(T2D_P,CVD_P)
  loci <- jn[, {
    idx <- which.min(T2D_P)
    .(index_SNP = T2D_SNP[idx], CHR = CHR[idx], index_BP = T2D_BP[idx],
      index_T2D_P = T2D_P[idx], locus_joint_p = min(joint_p))
  }, by = locus]
  loci[, gene := mapply(nearest_gene, CHR, index_BP)]
  setorder(loci, CHR, index_BP)               # genomic order for the list
  loci_list <- paste(unique(loci$gene), collapse = ", ")
  strongest <- loci[which.min(locus_joint_p), gene]

  list(n_shared_snps = uniqueN(jn$T2D_SNP), n_loci = nrow(loci),
       strongest = strongest, loci_list = loci_list,
       detail = jn[order(distance_bp)][, .SD[1], by = T2D_SNP][
         , .(T2D_SNP, CHR, T2D_BP, T2D_P, CVD_trait = toupper(tr),
             nearest_CVD_SNP = CVD_SNP, CVD_BP, CVD_P, distance_bp,
             locus, gene = loci$gene[match(locus, loci$locus)])])
})

# --- Summary table ---
summary_rows <- rbindlist(lapply(CVD_TRAITS, function(tr) {
  p <- per_trait[[tr]]
  data.table(
    n_T2D_lead_SNPs       = nrow(t2d),
    n_T2D_loci            = n_t2d_loci,
    CVD_trait             = toupper(tr),
    n_CVD_lead_SNPs       = nrow(leads[[tr]]),
    n_CVD_loci            = count_loci(leads[[tr]], LOCUS_BP),
    n_shared_lead_SNPs    = p$n_shared_snps,
    n_shared_loci         = p$n_loci,
    strongest_shared_locus = p$strongest,
    shared_loci           = p$loci_list
  )
}))

# --- Detail sheet ---
detail <- rbindlist(lapply(CVD_TRAITS, function(tr) per_trait[[tr]]$detail), fill = TRUE)
if (!is.null(detail) && nrow(detail) > 0) setorder(detail, CHR, T2D_BP, CVD_trait)

# --- Write workbook ---
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold")
add_sheet <- function(name, df) {
  addWorksheet(wb, name); writeData(wb, name, df, headerStyle = hs)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
}
add_sheet("Summary", summary_rows)
if (!is.null(detail) && nrow(detail) > 0) add_sheet("Shared_T2D_leads", detail)

# Parameter provenance sheet
params <- data.frame(parameter = c(
  "lead SNP clumping", "MHC", "LD reference",
  "shared lead SNP rule", "shared locus merge window",
  "strongest shared locus", "locus label",
  "T2D GWAS", "CVD GWAS"),
  value = c(
  "plink --clump P1=5e-8, P2=1, r2>0.05, 500kb",
  "excluded (chr6:25-35Mb, GRCh37)", "1000G EUR (n=503)",
  sprintf("T2D lead within %g kb of a CVD lead", SHARE_KB),
  sprintf("merge shared T2D leads within %g kb", LOCUS_KB),
  "locus minimising max(T2D_P, CVD_P) over its shared lead pairs",
  "nearest protein-coding gene (GENCODE v19/GRCh37) to the locus index SNP (lowest-T2D-P shared lead)",
  "Suzuki Nature 2024 EUR",
  "Aragam 2022 (CAD); Verma 2024 (MI/Angina/PAD); Mishra 2022 (Stroke)"))
add_sheet("Parameters", params)

saveWorkbook(wb, out_file, overwrite = TRUE)
cat(sprintf("Wrote: %s\n", out_file))
cat(sprintf("\nParameters: lead clumping r2>0.05/500kb; shared = T2D lead within %g kb of CVD lead; loci merged within %g kb\n",
            SHARE_KB, LOCUS_KB))
cat("\n=== Summary ===\n"); print(summary_rows, row.names = FALSE)
