#!/usr/bin/env Rscript
# ============================================================================
# a1_3: Cross-ancestry GWAS signal consistency (EUR vs EAS) for CAD and T2D
# ============================================================================
# Compares ancestry-specific GWAS arms to assess consistency of signals between
# EUR and EAS. Per trait (EUR = reference):
#   (1) Lead-SNP definition per arm (genome-wide sig, ancestry-matched LD clump).
#   (2) Locus sharing within +/-500 kb -> EUR-only / shared / EAS-only bar chart.
#       (Only 2 ancestries -> a 3-category bar, not an UpSet; UpSet would only
#        pay off with >=3 ancestries.)
#   (3) Effect-size concordance: lm(beta_EUR ~ beta_EAS) at union lead SNPs,
#       with R2 + 95% CI; scatter, EAS-unique-significant SNPs highlighted.
#   (4) Allele-frequency concordance: harmonized EAF(EUR) vs EAF(EAS) scatter.
#
# Comparisons (single-ancestry arms only; trans-ancestry META_ files excluded):
#   CAD: Aragam_NatureGenetics_2022.CAD.EUR  vs  Sakaue_NatGenet_2021.CAD.EAS  (cross-study)
#   T2D: Suzuki_Nature_2024.t2d.EUR          vs  Suzuki_Nature_2024.t2d.EAS    (single consortium)
#
# Framing: validation + mechanism illustration (shared biology -> concordant beta;
# allele-freq divergence / power / LD differences -> ancestry-unique loci),
# NOT discovery of ancestry-specific biology. CAD additionally has cross-study
# heterogeneity; unique-locus counts are confounded by EUR's larger sample size.
#
# Usage:
#   Rscript scripts/a1_3_analysis/01_cross_ancestry_consistency.R [--trait CAD|T2D|both]
# (Run under module R/4.2.0 + plink, or the rocky9 R 4.2.0 recipe.)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
setwd(base_dir)

# Reuse existing QC / clumping / IO helpers (do not reimplement)
source(file.path(base_dir, "scripts/a1_analysis/prep_bnmf.R"))

# -- Config -------------------------------------------------------------------
sm_dir   <- file.path(base_dir, "sumstats/harmonized")
out_dir  <- file.path(base_dir, "results/a1_3_analysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Analysis parameters
P_THRESH   <- 5e-8
MAF_THRESH <- 0.01
CLUMP_R2   <- 0.05
CLUMP_KB   <- 500
OVERLAP_KB <- 500            # locus-sharing distance window (each side)
PANEL_EUR  <- file.path(base_dir, "reference/1kg_eur/1000G.EUR.QC")
PANEL_EAS  <- file.path(base_dir, "reference/1kg_eas/1000G.EAS.QC")

# Trait pairs: EUR (reference) vs EAS
comparisons <- list(
  CAD = list(
    EUR = list(file = file.path(sm_dir, "Aragam_NatureGenetics_2022.CAD.EUR.GRCh37.processed.txt.gz"),
               panel = PANEL_EUR, study = "Aragam 2022"),
    EAS = list(file = file.path(sm_dir, "Sakaue_NatGenet_2021.CAD.EAS.GRCh37.processed.txt.gz"),
               panel = PANEL_EAS, study = "Sakaue 2021"),
    note = "cross-study (Aragam EUR vs Sakaue EAS); Sakaue MAF/EAF may be sparse"
  ),
  T2D = list(
    EUR = list(file = file.path(sm_dir, "Suzuki_Nature_2024.t2d.EUR.GRCh37.processed.txt.gz"),
               panel = PANEL_EUR, study = "Suzuki 2024"),
    EAS = list(file = file.path(sm_dir, "Suzuki_Nature_2024.t2d.EAS.GRCh37.processed.txt.gz"),
               panel = PANEL_EAS, study = "Suzuki 2024"),
    note = "single consortium (Suzuki 2024 EUR vs EAS)"
  )
)

# -- CLI ----------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
trait_arg <- "both"
if ("--trait" %in% args) trait_arg <- args[which(args == "--trait") + 1]
traits <- if (toupper(trait_arg) == "BOTH") names(comparisons) else toupper(trait_arg)

message("=== a1_3: Cross-ancestry GWAS consistency (EUR vs EAS) ===\n")
message("Parameters:")
message(sprintf("  P threshold        : %.0e", P_THRESH))
message(sprintf("  MAF threshold      : %.2f", MAF_THRESH))
message(sprintf("  LD clump           : r2<%.2f, %d kb", CLUMP_R2, CLUMP_KB))
message(sprintf("  Locus overlap      : +/-%d kb", OVERLAP_KB))
message(sprintf("  LD panels          : EUR=%s ; EAS=%s", basename(PANEL_EUR), basename(PANEL_EAS)))
message(sprintf("  HapMap3 restriction: NONE; MHC + strand-ambiguous excluded (qc_variants defaults)"))
message(sprintf("  Traits             : %s\n", paste(traits, collapse = ", ")))

# -- Helpers ------------------------------------------------------------------

# Parse CHR/POS/REF/ALT from VAR_ID (chr_pos_ref_alt)
parse_varid <- function(dt) {
  p <- tstrsplit(dt$VAR_ID, "_", fixed = TRUE)
  dt[, `:=`(CHR = as.integer(p[[1]]), POS = as.integer(p[[2]]),
            REF = p[[3]], ALT = p[[4]])]
  dt[]
}

# Lead SNPs for one arm (full sumstat rows for clumped genome-wide-sig variants)
get_leads <- function(arm) {
  res <- qc_variants(arm$file, p_threshold = P_THRESH, maf_threshold = MAF_THRESH,
                     ref_panel_prefix = arm$panel, clump_r2 = CLUMP_R2,
                     clump_kb = CLUMP_KB, hapmap3_file = NULL)
  parse_varid(res$data)
}

# Classify locus sharing between two lead sets by distance (+/-OVERLAP_KB).
# A lead is "shared" if a lead on the other side sits within the window.
classify_sharing <- function(eur, eas) {
  win <- OVERLAP_KB * 1000L
  e <- eur[, .(VAR_ID, CHR, POS)]
  a <- eas[, .(VAR_ID, CHR, POS)]
  # EUR leads with an EAS lead within window
  a2 <- copy(a)[, `:=`(lo = POS - win, hi = POS + win)]
  eur_shared <- e[a2, on = .(CHR, POS >= lo, POS <= hi), nomatch = NULL,
                  .(VAR_ID = x.VAR_ID), allow.cartesian = TRUE]
  e2 <- copy(e)[, `:=`(lo = POS - win, hi = POS + win)]
  eas_shared <- a[e2, on = .(CHR, POS >= lo, POS <= hi), nomatch = NULL,
                  .(VAR_ID = x.VAR_ID), allow.cartesian = TRUE]
  list(
    eur_only  = setdiff(e$VAR_ID, unique(eur_shared$VAR_ID)),
    eas_only  = setdiff(a$VAR_ID, unique(eas_shared$VAR_ID)),
    shared_eur = unique(eur_shared$VAR_ID),
    shared_eas = unique(eas_shared$VAR_ID)
  )
}

# Read full sumstats (selected cols) and subset to a VAR_ID set
read_subset <- function(file, varids) {
  dt <- fread(file, select = c("VAR_ID", "Effect_Allele", "BETA", "SE", "P_VALUE", "MAF", "EAF"))
  dt[VAR_ID %in% varids]
}

# -- Plotting (loaded after data.table work to limit dplyr masking) -----------
suppressPackageStartupMessages({
  library(ggplot2)
})
source(file.path(base_dir, "scripts/a1_analysis/figure_utils.R"))  # theme_big_text()

PAL <- c("EUR only" = "#4682B4", "Shared" = "#6A3D9A", "EAS only" = "#F4845F")

summary_rows <- list()

for (tr in traits) {
  cmp <- comparisons[[tr]]
  message(sprintf("\n===== %s : %s vs %s =====", tr, cmp$EUR$study, cmp$EAS$study))
  message(sprintf("  (%s)", cmp$note))

  # (1) Lead SNPs per arm
  leads_eur <- get_leads(cmp$EUR)
  leads_eas <- get_leads(cmp$EAS)
  fwrite(leads_eur, file.path(out_dir, sprintf("lead_snps_%s_EUR.csv", tr)))
  fwrite(leads_eas, file.path(out_dir, sprintf("lead_snps_%s_EAS.csv", tr)))
  message(sprintf("  Lead SNPs: EUR=%d, EAS=%d", nrow(leads_eur), nrow(leads_eas)))

  # (2) Locus sharing (+/-500 kb)
  sh <- classify_sharing(leads_eur, leads_eas)
  n_eur_only <- length(sh$eur_only)
  n_eas_only <- length(sh$eas_only)
  n_shared   <- length(sh$shared_eur)   # EUR leads with EAS partner
  share_dt <- data.table(
    category = factor(c("EUR only", "Shared", "EAS only"),
                      levels = c("EUR only", "Shared", "EAS only")),
    n = c(n_eur_only, n_shared, n_eas_only)
  )
  fwrite(share_dt, file.path(out_dir, sprintf("locus_sharing_%s.csv", tr)))
  message(sprintf("  Locus sharing (+/-%dkb): EUR-only=%d, shared=%d, EAS-only=%d",
                  OVERLAP_KB, n_eur_only, n_shared, n_eas_only))

  p_bar <- ggplot(share_dt, aes(category, n, fill = category)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n), vjust = -0.3, size = 5) +
    scale_fill_manual(values = PAL, guide = "none") +
    labs(x = NULL, y = "Lead loci",
         title = sprintf("%s lead loci: EUR vs EAS sharing (+/-%dkb)", tr, OVERLAP_KB),
         subtitle = sprintf("%s | EUR n_lead=%d, EAS n_lead=%d (EUR power >> EAS)",
                            cmp$note, nrow(leads_eur), nrow(leads_eas))) +
    theme_big_text(base_size = 14)
  ggsave(file.path(out_dir, sprintf("barplot_locus_sharing_%s.png", tr)),
         p_bar, width = 7, height = 6, dpi = 300)

  # (3)+(4) Concordance on the ASCERTAINED set: EAS lead SNPs with EUR effects
  # looked up (the published design — "select lead SNPs in the non-EUR ancestry,
  # extract beta/SE from EUR"). Regressing over the union of leads instead dilutes
  # the relationship because EUR-only leads are null/noisy in underpowered EAS.
  eur_at_eas <- read_subset(cmp$EUR$file, leads_eas$VAR_ID)
  setnames(eur_at_eas, setdiff(names(eur_at_eas), "VAR_ID"),
           paste0(setdiff(names(eur_at_eas), "VAR_ID"), "_EUR"))
  eas_sub <- leads_eas[, .(VAR_ID, Effect_Allele_EAS = Effect_Allele,
                           BETA_EAS = BETA, SE_EAS = SE, EAF_EAS = EAF)]
  m <- merge(eas_sub, eur_at_eas, by = "VAR_ID")   # EAS leads present in EUR sumstats
  message(sprintf("  EAS leads present in EUR sumstats: %d / %d", nrow(m), nrow(leads_eas)))

  # Harmonize EAS to EUR effect allele: flip beta sign + EAF complement when alleles differ
  flip <- m$Effect_Allele_EUR != m$Effect_Allele_EAS
  m[flip, `:=`(BETA_EAS = -BETA_EAS, EAF_EAS = 1 - EAF_EAS)]
  m[, also_eur_lead := VAR_ID %in% leads_eur$VAR_ID]
  m[, lead_in := factor(fifelse(also_eur_lead, "Shared", "EAS only"),
                        levels = c("Shared", "EAS only"))]

  # Allele-orientation QC: after harmonizing, effect alleles agree but betas should
  # be positively concordant. If they systematically anti-correlate, the two files'
  # BETA columns are oriented to opposite alleles despite matching Effect_Allele
  # labels (a harmonization inconsistency) -> the beta comparison is NOT interpretable.
  mb <- m[is.finite(BETA_EUR) & is.finite(BETA_EAS)]
  sign_conc <- if (nrow(mb) > 0) mean(sign(mb$BETA_EUR) == sign(mb$BETA_EAS)) else NA_real_
  allele_ok <- is.na(sign_conc) || sign_conc >= 0.5
  if (!allele_ok) {
    message(sprintf("  *** ALLELE-ORIENTATION WARNING: sign concordance %.0f%% (n=%d) — ",
                    100 * sign_conc, nrow(mb)),
            "BETA columns appear oriented to opposite alleles between the two studies. ",
            "Beta concordance for this trait is NOT interpretable (likely cross-study harmonization mismatch).")
  }

  fit <- lm(BETA_EUR ~ BETA_EAS, data = mb)
  r2 <- summary(fit)$r.squared
  slope <- coef(fit)[["BETA_EAS"]]
  ci <- tryCatch(confint(fit)["BETA_EAS", ], error = function(e) c(NA, NA))
  message(sprintf("  Beta concordance (EAS leads, n=%d): slope=%.3f [%.3f, %.3f], R2=%.3f, sign-conc=%.0f%%",
                  nrow(mb), slope, ci[1], ci[2], r2, 100 * sign_conc))

  flag_txt <- if (allele_ok) "" else " | ALLELE-ORIENTATION MISMATCH — not interpretable"
  lim <- max(abs(c(mb$BETA_EUR, mb$BETA_EAS)), na.rm = TRUE)
  p_beta <- ggplot(mb, aes(BETA_EAS, BETA_EUR)) +
    geom_hline(yintercept = 0, color = "grey80") +
    geom_vline(xintercept = 0, color = "grey80") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = lead_in), alpha = 0.8, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
    scale_color_manual(values = PAL[c("Shared", "EAS only")], name = "Lead in") +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(x = sprintf("EAS beta (%s, EAS lead SNPs)", cmp$EAS$study),
         y = sprintf("EUR beta (%s, looked up)", cmp$EUR$study),
         title = sprintf("%s effect-size concordance: EAS leads vs EUR", tr),
         subtitle = sprintf("slope=%.2f [%.2f, %.2f], R2=%.2f, n=%d%s",
                            slope, ci[1], ci[2], r2, nrow(mb), flag_txt)) +
    theme_big_text(base_size = 14)
  ggsave(file.path(out_dir, sprintf("scatter_beta_%s.png", tr)),
         p_beta, width = 8, height = 8, dpi = 300)

  # Allele-frequency concordance (harmonized EAF) at EAS lead SNPs
  ma <- m[is.finite(EAF_EUR) & is.finite(EAF_EAS)]
  af_cov <- nrow(ma)
  if (af_cov >= 5) {
    p_af <- ggplot(ma, aes(EAF_EUR, EAF_EAS)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(aes(color = lead_in), alpha = 0.8, size = 2) +
      scale_color_manual(values = PAL[c("Shared", "EAS only")], name = "Lead in") +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(x = "EUR effect-allele freq", y = "EAS effect-allele freq",
           title = sprintf("%s allele-frequency divergence: EUR vs EAS", tr),
           subtitle = sprintf("%d EAS lead SNPs with freq in both arms (dashed=identity)", af_cov)) +
      theme_big_text(base_size = 14)
    ggsave(file.path(out_dir, sprintf("scatter_AF_%s.png", tr)),
           p_af, width = 8, height = 8, dpi = 300)
  } else {
    message(sprintf("  AF scatter SKIPPED: only %d EAS lead SNPs have freq in both arms (sparse).", af_cov))
  }

  summary_rows[[tr]] <- data.table(
    trait = tr, eur_study = cmp$EUR$study, eas_study = cmp$EAS$study,
    n_lead_EUR = nrow(leads_eur), n_lead_EAS = nrow(leads_eas),
    n_EUR_only = n_eur_only, n_shared = n_shared, n_EAS_only = n_eas_only,
    n_beta = nrow(mb), beta_slope = slope, beta_ci_lo = ci[1], beta_ci_hi = ci[2],
    beta_r2 = r2, beta_sign_conc = sign_conc, beta_interpretable = allele_ok,
    n_AF = af_cov
  )

  rm(eur_at_eas, m); gc()
}

# -- Combined summary ---------------------------------------------------------
summ <- rbindlist(summary_rows, use.names = TRUE)
fwrite(summ, file.path(out_dir, "cross_ancestry_summary.csv"))
message("\nSummary:")
print(summ)
message(sprintf("\nOutputs in %s", out_dir))
message("Done.")
