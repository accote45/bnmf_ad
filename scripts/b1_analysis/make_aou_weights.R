#!/usr/bin/env Rscript
# make_aou_weights.R
# Build a single, self-contained cluster-PGS weights file for external
# replication of B1 (e.g. in All of Us). Reproduces the weight-building block of
# 01_compute_cluster_prs.R exactly: per variant BETA = max(|BETA_T2D|,|BETA_CAD|),
# then Wbeta_Kj = BETA * W_Kj. Adds GRCh38 coordinates via UCSC liftOver so the
# weights (native GRCh37) can be matched against GRCh38 genotypes by rsID first
# and lifted chr:pos:allele as a fallback.
#
# Usage:
#   ml liftover            # provides liftOver on PATH (rocky9)
#   Rscript scripts/b1_analysis/make_aou_weights.R
#   Rscript scripts/b1_analysis/make_aou_weights.R --config config/b1_config.yaml

library(data.table)
library(yaml)

# strongest_beta() — per-variant strongest-signal BETA across a trait's GWAS list
source("scripts/b1_analysis/prs_scoring_helpers.R")

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b1_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}

cfg <- read_yaml(config_path)

w_matrix_path <- cfg$a1_results$w_matrix
filt_var_path <- cfg$a1_results$filtered_variants
t2d_files     <- unlist(cfg$gwas$t2d_files)
cad_files     <- unlist(cfg$gwas$cad_files)
results_dir   <- cfg$results_dir

out_dir   <- file.path(results_dir, "aou")
chain_gz  <- file.path(out_dir, "hg19ToHg38.over.chain.gz")
out_path  <- file.path(out_dir, "b1_cluster_weights.tsv")
liftover_bin <- "liftOver"   # from `ml liftover`
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Report all parameters before doing any work (CLAUDE.md) ---
cat("=== make_aou_weights: parameters ===\n")
cat(sprintf("  W matrix:        %s\n", w_matrix_path))
cat(sprintf("  Filtered var:    %s\n", filt_var_path))
cat(sprintf("  T2D GWAS (%d):   %s\n", length(t2d_files), paste(basename(t2d_files), collapse = ", ")))
cat(sprintf("  CAD GWAS (%d):   %s\n", length(cad_files), paste(basename(cad_files), collapse = ", ")))
cat(sprintf("  liftOver chain:  %s\n", chain_gz))
cat(sprintf("  liftOver binary: %s\n", liftover_bin))
cat(sprintf("  Output:          %s\n", out_path))

stopifnot(file.exists(chain_gz))

# --- Load W matrix ---
w_dt <- fread(w_matrix_path)
cluster_cols <- setdiff(colnames(w_dt), "VAR_ID")
cat(sprintf("\nW matrix: %d variants x %d clusters (%s)\n",
            nrow(w_dt), length(cluster_cols), paste(cluster_cols, collapse = ", ")))

# --- Load filtered variants (RSID + Effect_Allele) ---
filt_dt <- fread(filt_var_path)
cat(sprintf("Filtered variants: %d rows\n", nrow(filt_dt)))

# --- Strongest-signal BETA per variant across each trait's GWAS list ---
# (matches 01_compute_cluster_prs.R; restrict to the W-matrix universe)
cat("\n--- Loading GWAS for BETA lookup (strongest-signal across list) ---\n")
universe_ids <- w_dt$VAR_ID
t2d <- strongest_beta(universe_ids, t2d_files)[, .(VAR_ID, BETA_T2D = BETA)]
cat(sprintf("  T2D: %d/%d universe variants with BETA\n", nrow(t2d), length(universe_ids)))

cad <- strongest_beta(universe_ids, cad_files)[, .(VAR_ID, BETA_CAD = BETA)]
cat(sprintf("  CAD: %d/%d universe variants with BETA\n", nrow(cad), length(universe_ids)))

# --- Merge W matrix with variant info and GWAS betas (same as script 01) ---
merged <- merge(w_dt, filt_dt[, .(VAR_ID, RSID, Effect_Allele)], by = "VAR_ID")
merged <- merge(merged, t2d, by = "VAR_ID", all.x = TRUE)
merged <- merge(merged, cad, by = "VAR_ID", all.x = TRUE)

cat(sprintf("\nMerged: %d variants\n", nrow(merged)))
cat(sprintf("  With T2D BETA: %d\n", sum(!is.na(merged$BETA_T2D))))
cat(sprintf("  With CAD BETA: %d\n", sum(!is.na(merged$BETA_CAD))))

# --- Select BETA with max |BETA| (identical cascade to 01_compute_cluster_prs.R) ---
merged[, BETA := fifelse(
  is.na(BETA_T2D) & is.na(BETA_CAD), NA_real_,
  fifelse(is.na(BETA_T2D), BETA_CAD,
  fifelse(is.na(BETA_CAD), BETA_T2D,
  fifelse(abs(BETA_T2D) >= abs(BETA_CAD), BETA_T2D, BETA_CAD)))
)]

merged[, beta_source := fifelse(
  is.na(BETA_T2D) & is.na(BETA_CAD), "none",
  fifelse(is.na(BETA_T2D), "CAD",
  fifelse(is.na(BETA_CAD), "T2D",
  fifelse(abs(BETA_T2D) >= abs(BETA_CAD), "T2D", "CAD")))
)]

cat(sprintf("  BETA source: T2D=%d, CAD=%d, none=%d\n",
            sum(merged$beta_source == "T2D"),
            sum(merged$beta_source == "CAD"),
            sum(merged$beta_source == "none")))

# Drop variants with no BETA
merged <- merged[!is.na(BETA)]
cat(sprintf("  Variants with BETA: %d\n", nrow(merged)))

# --- Compute weighted scores: Wbeta_Kj = BETA * W_Kj ---
wbeta_cols <- paste0("Wbeta_", cluster_cols)
for (j in seq_along(cluster_cols)) {
  merged[[wbeta_cols[j]]] <- merged$BETA * merged[[cluster_cols[j]]]
}

# Keep variants with at least one non-zero weighted score
keep <- rowSums(abs(merged[, ..wbeta_cols])) > 0
merged <- merged[keep]
cat(sprintf("  Variants with non-zero weights: %d\n", nrow(merged)))

# --- Parse CHR37/POS37/REF/ALT from VAR_ID (CHR_POS_REF_ALT, GRCh37) ---
merged[, c("CHR37", "POS37", "REF", "ALT") := tstrsplit(VAR_ID, "_", fixed = TRUE)]
merged[, POS37 := as.integer(POS37)]

# --- liftOver GRCh37 -> GRCh38 ---------------------------------------------
# Build a 0-based BED (chr, start, end, name=VAR_ID), run liftOver, read back
# the mapped/unmapped coordinates. liftOver needs "chr"-prefixed contigs.
cat("\n--- liftOver GRCh37 -> GRCh38 ---\n")
bed_in   <- file.path(out_dir, ".lift_in.bed")
bed_out  <- file.path(out_dir, ".lift_out.bed")
bed_unmap<- file.path(out_dir, ".lift_unmapped.bed")

bed_dt <- merged[, .(chr = paste0("chr", CHR37),
                     start = POS37 - 1L,
                     end = POS37,
                     name = VAR_ID)]
fwrite(bed_dt, bed_in, sep = "\t", col.names = FALSE)

lift_cmd <- sprintf("%s %s %s %s %s",
                    liftover_bin, bed_in, chain_gz, bed_out, bed_unmap)
cat(sprintf("  %s\n", lift_cmd))
ret <- system(lift_cmd)
stopifnot(ret == 0, file.exists(bed_out))

lifted <- fread(bed_out, header = FALSE,
                col.names = c("chr38", "start38", "end38", "VAR_ID"))
lifted[, CHR38 := sub("^chr", "", chr38)]
lifted[, POS38 := end38]                       # 1-based GRCh38 position
lifted <- lifted[, .(VAR_ID, CHR38, POS38)]

n_unmapped <- nrow(merged) - nrow(lifted)
cat(sprintf("  liftOver: %d mapped, %d unmapped (dropped)\n",
            nrow(lifted), n_unmapped))

merged <- merge(merged, lifted, by = "VAR_ID")  # inner: keep only lifted variants

# --- Assemble and write the consolidated weights file ----------------------
out_cols <- c("VAR_ID", "CHR37", "POS37", "CHR38", "POS38",
              "RSID", "REF", "ALT", "Effect_Allele", "BETA", "beta_source",
              wbeta_cols)
out_dt <- merged[, ..out_cols]
setnames(out_dt, "VAR_ID", "VAR_ID37")
# Stable order by genomic position (GRCh38)
out_dt[, .chr_n := suppressWarnings(as.integer(CHR38))]
setorder(out_dt, .chr_n, POS38, na.last = TRUE)
out_dt[, .chr_n := NULL]

fwrite(out_dt, out_path, sep = "\t")
cat(sprintf("\nOutput: %s (%d variants x %d cluster weights)\n",
            out_path, nrow(out_dt), length(wbeta_cols)))

# --- Sanity summary ---
cat("\n--- Weighted-score summary ---\n")
for (wc in wbeta_cols) {
  v <- out_dt[[wc]]
  cat(sprintf("  %s: mean=%.6g sd=%.6g NAs=%d\n",
              wc, mean(v, na.rm = TRUE), sd(v, na.rm = TRUE), sum(is.na(v))))
}
cat(sprintf("\n  rsID present: %d / %d\n",
            sum(out_dt$RSID != "." & !is.na(out_dt$RSID)), nrow(out_dt)))

# --- Clean temp liftOver files ---
file.remove(c(bed_in, bed_out, bed_unmap))

cat("\n=== Done ===\n")
