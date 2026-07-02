#!/usr/bin/env Rscript
# 01b_compute_genomewide_prs.R
# Genome-wide PRS comparator for B1, built the SAME way as the bNMF variant
# universe so the comparison with the cluster PRS is fair: same input GWAS
# (4 CAD + 5 T2D for META), same QC + LD clumping (5e-8, r2=0.05, kb=500 on the
# 1000G EUR panel), via prep_bnmf.R::qc_variants_multi(). Two flavours:
#
#   (A) Trait-specific  — union-clump each trait's GWAS list on its own:
#         GW_T2D  = sum(BETA_T2D * dosage) over T2D-significant variants
#         GW_CAD  = sum(BETA_CAD * dosage) over CAD-significant variants
#
#   (B) Combined-variant (supplementary) — fix the variant set to the full
#       pooled T2D+CAD union universe (= the bNMF cluster variant set,
#       filtered_variants_<ancestry>.tsv) and vary only the weighting:
#         GW_T2D_combined = sum(BETA_T2D * dosage) over the shared universe
#         GW_CAD_combined = sum(BETA_CAD * dosage) over the shared universe
#       Labelled "GW T2D (combined variants)" / "GW CAD (combined variants)".
#
# Per-variant BETA = strongest-signal (min-P) effect across the trait's GWAS
# list (strongest_beta() in prs_scoring_helpers.R), matching the cluster-PRS
# BETAs in 01_compute_cluster_prs.R.
#
# Output: results/b1_analysis/prs/genomewide_prs_all.tsv
#         (FID, IID, GW_T2D, GW_CAD, GW_T2D_combined, GW_CAD_combined, group)
#
# Usage:
#   Rscript scripts/b1_analysis/01b_compute_genomewide_prs.R
#   Rscript scripts/b1_analysis/01b_compute_genomewide_prs.R --config config/b1_config.yaml

library(data.table)
library(yaml)

# bNMF variant pipeline (qc_variants_multi) + shared B1 scorer / BETA helpers.
source("scripts/a1_analysis/prep_bnmf.R")
source("scripts/b1_analysis/prs_scoring_helpers.R")

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b1_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}

cfg <- read_yaml(config_path)

t2d_files   <- unlist(cfg$gwas$t2d_files)
cad_files   <- unlist(cfg$gwas$cad_files)
filt_var    <- cfg$a1_results$filtered_variants
geno_prefix <- cfg$genotypes$prefix
results_dir <- cfg$results_dir

# Clumping / variant-selection params (mirror config/a1_config.yaml for META).
p_threshold   <- as.numeric(cfg$ld_clump$p_threshold)
maf_threshold <- as.numeric(cfg$ld_clump$maf_threshold)
clump_r2      <- as.numeric(cfg$ld_clump$r2)
clump_kb      <- as.integer(cfg$ld_clump$kb)
ref_panel     <- cfg$ld_clump$ref_panel
union_clump   <- isTRUE(cfg$ld_clump$union_clump)
plink_bin     <- cfg$ld_clump$plink_bin
hapmap3_path  <- cfg$ld_clump$hapmap3_snp_file
if (!is.null(hapmap3_path) && !file.exists(hapmap3_path)) {
  cat(sprintf("WARNING: HapMap3 file not found at %s, proceeding without restriction.\n",
              hapmap3_path))
  hapmap3_path <- NULL
}

plink2_bin <- "/hpc/packages/minerva-centos7/plink2/2.3/plink2"

cat("=== B1 Step 1b: Compute Genome-Wide PRS ===\n")
cat(sprintf("  T2D GWAS (%d): %s\n", length(t2d_files), paste(basename(t2d_files), collapse = ", ")))
cat(sprintf("  CAD GWAS (%d): %s\n", length(cad_files), paste(basename(cad_files), collapse = ", ")))
cat(sprintf("  Combined universe: %s\n", filt_var))
cat(sprintf("  P<%.0e | MAF>=%.4f | clump r2<%g, %dkb | ref %s | union_clump=%s\n",
            p_threshold, maf_threshold, clump_r2, clump_kb, ref_panel, union_clump))

prs_dir    <- file.path(results_dir, "prs")
sample_dir <- if (!is.null(cfg$samples_dir)) cfg$samples_dir else file.path(results_dir, "samples")
dir.create(prs_dir, recursive = TRUE, showWarnings = FALSE)

keep_files <- list(
  eur_train      = file.path(sample_dir, "eur_train.keep"),
  eur_validation = file.path(sample_dir, "eur_validation.keep"),
  afr_validation = file.path(sample_dir, "afr_validation.keep"),
  eas_validation = file.path(sample_dir, "eas_validation.keep"),
  sas_validation = file.path(sample_dir, "sas_validation.keep")
)

# --- Build a trait's union-clumped, genome-wide-significant variant set -------
# Mirrors 01_run_bnmf.R: per-file QC + clump, union (clump for META), then a
# post-hoc HapMap3 restriction for non-union ancestries.
build_universe <- function(files) {
  rp <- if (!is.null(ref_panel) && file.exists(paste0(ref_panel, ".1.bed"))) ref_panel else NULL
  if (is.null(rp)) cat(sprintf("WARNING: LD reference panel not found at %s\n", ref_panel))
  qc <- qc_variants_multi(files, p_threshold = p_threshold,
                          maf_threshold = maf_threshold,
                          ref_panel_prefix = rp,
                          clump_r2 = clump_r2, clump_kb = clump_kb,
                          plink_bin = plink_bin,
                          union_clump = union_clump,
                          union_hapmap3_file = if (union_clump) hapmap3_path else NULL)
  filtered <- qc$data
  if (!is.null(hapmap3_path) && !union_clump) {
    hm3 <- fread(hapmap3_path, select = "VAR_ID")
    filtered <- filtered[VAR_ID %in% hm3$VAR_ID]
  }
  filtered
}

# --- Score one variant set with a single weight column -----------------------
# vars must have VAR_ID, Effect_Allele, and a numeric column named `weight_col`.
# Each set scores into its own subdirectory so the per-chr .sscore/.cachekey
# files never collide with the cluster PRS or the other GW sets.
score_set <- function(vars, weight_col, out_name, tag) {
  vars <- copy(vars)
  vars[, c("v_chr", "v_pos") := tstrsplit(VAR_ID, "_", keep = 1:2)]
  vars[, v_chr := as.integer(v_chr)]
  vars[, v_pos := as.integer(v_pos)]
  vars <- vars[!is.na(get(weight_col))]

  set_dir <- file.path(prs_dir, paste0("gw_", tag))
  dir.create(set_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("\n--- Scoring set '%s' (%s): %d variants ---\n",
              tag, out_name, nrow(vars)))
  fwrite(vars[, .(VAR_ID, Effect_Allele, weight = get(weight_col))],
         file.path(set_dir, sprintf("gw_variant_weights_%s.tsv", tag)), sep = "\t")

  chr_score_paths <- write_chr_score_files(vars, weight_col, geno_prefix, set_dir, tag = tag)
  score_by_group(chr_score_paths, 1L, out_name, keep_files, geno_prefix, plink2_bin, set_dir)
}

# ============================================================
# (A) Trait-specific genome-wide PRS
# ============================================================
cat("\n=== Building trait-specific variant universes ===\n")
t2d_uni <- build_universe(t2d_files)
cad_uni <- build_universe(cad_files)
cat(sprintf("\n  GW T2D universe: %d variants | GW CAD universe: %d variants\n",
            nrow(t2d_uni), nrow(cad_uni)))

# Weight each universe by the strongest-signal trait BETA (consistent with the
# cluster PRS); keep the harmonized Effect_Allele from the universe table.
t2d_w <- strongest_beta(t2d_uni$VAR_ID, t2d_files)[, .(VAR_ID, BETA_T2D = BETA)]
cad_w <- strongest_beta(cad_uni$VAR_ID, cad_files)[, .(VAR_ID, BETA_CAD = BETA)]
t2d_uni <- merge(t2d_uni[, .(VAR_ID, Effect_Allele)], t2d_w, by = "VAR_ID")
cad_uni <- merge(cad_uni[, .(VAR_ID, Effect_Allele)], cad_w, by = "VAR_ID")

d_t2d <- score_set(t2d_uni, "BETA_T2D", "GW_T2D", "t2d")
d_cad <- score_set(cad_uni, "BETA_CAD", "GW_CAD", "cad")

# ============================================================
# (B) Combined-variant (shared universe) genome-wide PRS
# ============================================================
cat("\n=== Building combined-variant PRS over the shared bNMF universe ===\n")
universe <- fread(filt_var)[, .(VAR_ID, Effect_Allele)]
cat(sprintf("  Shared universe: %d variants (%s)\n", nrow(universe), basename(filt_var)))

t2d_cw <- strongest_beta(universe$VAR_ID, t2d_files)[, .(VAR_ID, BETA_T2D = BETA)]
cad_cw <- strongest_beta(universe$VAR_ID, cad_files)[, .(VAR_ID, BETA_CAD = BETA)]
uni_t2d <- merge(universe, t2d_cw, by = "VAR_ID")
uni_cad <- merge(universe, cad_cw, by = "VAR_ID")
cat(sprintf("  Shared universe with T2D BETA: %d | with CAD BETA: %d\n",
            nrow(uni_t2d), nrow(uni_cad)))

d_t2d_c <- score_set(uni_t2d, "BETA_T2D", "GW_T2D_combined", "t2d_combined")
d_cad_c <- score_set(uni_cad, "BETA_CAD", "GW_CAD_combined", "cad_combined")

# ============================================================
# Merge all four genome-wide scores per person
# ============================================================
cat("\n--- Combining genome-wide PRS scores ---\n")
score_tables <- list(d_t2d, d_cad, d_t2d_c, d_cad_c)
gw_all <- Reduce(function(a, b) merge(a, b, by = c("FID", "IID", "group"), all = TRUE),
                 score_tables)

# Column order: identifiers, four GW scores, group
setcolorder(gw_all, c("FID", "IID",
                      "GW_T2D", "GW_CAD", "GW_T2D_combined", "GW_CAD_combined",
                      "group"))

out_path <- file.path(prs_dir, "genomewide_prs_all.tsv")
fwrite(gw_all, out_path, sep = "\t")
cat(sprintf("\nOutput: %s\n", out_path))
cat(sprintf("Total individuals: %d | Groups: %s\n",
            nrow(gw_all), paste(unique(gw_all$group), collapse = ", ")))

# --- Summary ---
gw_cols <- c("GW_T2D", "GW_CAD", "GW_T2D_combined", "GW_CAD_combined")
cat("\n--- GW PRS summary (all samples) ---\n")
for (pc in gw_cols) {
  vals <- gw_all[[pc]]
  cat(sprintf("  %s: mean=%.6f, sd=%.6f, NAs=%d\n",
              pc, mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE), sum(is.na(vals))))
}

# --- Clean up per-chromosome plink2 logs (keep score/sscore/cachekey caches) ---
for (tag in c("t2d", "cad", "t2d_combined", "cad_combined")) {
  set_dir <- file.path(prs_dir, paste0("gw_", tag))
  chr_logs <- list.files(set_dir, pattern = "_chr\\d+\\.log$", full.names = TRUE)
  if (length(chr_logs) > 0) file.remove(chr_logs)
}

cat("\n=== Done ===\n")
