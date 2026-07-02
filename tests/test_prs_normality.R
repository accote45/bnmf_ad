#!/usr/bin/env Rscript
# ============================================================================
# test_prs_normality.R
# Quick PRS sanity check: compute PRS for EUR, apply weight cutoff,
# standardize, and verify approximate normality with histograms + QQ plots.
#
# Run on AoU:
#   Rscript tests/test_prs_normality.R              # chr16-19 only (fast)
#   Rscript tests/test_prs_normality.R --all-chr    # all 22 chr (full test)
#   Rscript tests/test_prs_normality.R --no-weight-cutoff  # skip cutoff
#
# Output:
#   tests/prs_normality_check/prs_normality_histograms.png
#   tests/prs_normality_check/prs_normality_qq.png
#   tests/prs_normality_check/weight_cutoff_eur.png  (diagnostic)
#   Console: summary stats + Shapiro-Wilk p-values
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# Parse simple flags
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
USE_ALL_CHR        <- "--all-chr" %in% args
APPLY_WEIGHT_CUTOFF <- !("--no-weight-cutoff" %in% args)

TEST_CHRS <- if (USE_ALL_CHR) 1:22 else 16:19
ANC       <- "EUR"

# ---------------------------------------------------------------------------
# Paths (AoU workspace)
# ---------------------------------------------------------------------------
ws  <- file.path(Sys.getenv("HOME"), "workspaces", "duplicateofstatinresponse")
mp  <- file.path(ws, "multiancestry_polygenic")

geno_dir      <- file.path(ws, "genotype", "plink")
geno_prefix   <- "chr"
results_dir   <- file.path(mp, "genotypes")
liftover_path <- file.path(ws, "genome_useful", "liftOver")
chain_path    <- file.path(ws, "genome_useful", "hg19ToHg38.over.chain.gz")
plink2_path   <- "plink2"

out_dir <- file.path(mp, "tests", "prs_normality_check")
dir.create(file.path(out_dir, "score_files"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "output"),       recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tmp"),          recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "diagnostics"),  recursive = TRUE, showWarnings = FALSE)

cat("=============================================\n")
cat(sprintf("PRS Normality Test — %s chr%d:%d\n", ANC, min(TEST_CHRS), max(TEST_CHRS)))
cat(sprintf("Weight cutoff: %s\n", if (APPLY_WEIGHT_CUTOFF) "enabled" else "disabled"))
cat("=============================================\n\n")

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
w_path  <- file.path(results_dir, ANC, sprintf("W_matrix_%s.tsv", ANC))
fv_path <- file.path(results_dir, ANC, sprintf("filtered_variants_%s.tsv", ANC))

stopifnot(
  "W_matrix not found"           = file.exists(w_path),
  "filtered_variants not found"  = file.exists(fv_path),
  "Genotype directory not found" = dir.exists(geno_dir),
  "liftOver binary not found"    = file.exists(liftover_path),
  "Chain file not found"         = file.exists(chain_path)
)

# ===========================================================================
# Weight cutoff helper (mirrors compute_prs.R::compute_weight_cutoff)
# ===========================================================================
compute_weight_cutoff <- function(w_dt, cluster_cols,
                                  top_pct = 0.01, bot_pct = 0.80) {

  all_weights <- unlist(w_dt[, cluster_cols, with = FALSE], use.names = FALSE)
  all_weights <- all_weights[all_weights > 0]
  n_total <- length(all_weights)

  empty_result <- list(
    cutoff = 0, n_total = n_total, n_above = n_total, n_below = 0,
    sorted_weights = sort(all_weights, decreasing = TRUE),
    signal_fit = NULL, noise_fit = NULL, crossover_rank = NA_integer_
  )

  if (n_total < 10) {
    warning("Too few non-zero weights (", n_total, "); skipping cutoff")
    return(empty_result)
  }

  sorted_w <- sort(all_weights, decreasing = TRUE)
  ranks    <- seq_along(sorted_w)

  n_top   <- max(2L, ceiling(top_pct * n_total))
  top_idx <- seq_len(n_top)
  signal_fit <- lm(sorted_w[top_idx] ~ ranks[top_idx])

  bot_start <- max(1L, floor((1 - bot_pct) * n_total) + 1L)
  bot_idx   <- seq(bot_start, n_total)
  noise_fit <- lm(sorted_w[bot_idx] ~ ranks[bot_idx])

  perp_dist <- function(x, y, fit) {
    a <- coef(fit)[[1]]
    b <- coef(fit)[[2]]
    abs(b * x - y + a) / sqrt(b^2 + 1)
  }

  d_signal <- perp_dist(ranks, sorted_w, signal_fit)
  d_noise  <- perp_dist(ranks, sorted_w, noise_fit)

  closer_to_signal <- which(d_signal < d_noise)

  if (length(closer_to_signal) == 0) {
    warning("No crossover found; all points closer to noise line. Setting cutoff = 0")
    empty_result$signal_fit <- signal_fit
    empty_result$noise_fit  <- noise_fit
    return(empty_result)
  }

  crossover_rank <- max(closer_to_signal)
  cutoff <- sorted_w[crossover_rank]

  list(
    cutoff         = cutoff,
    n_total        = n_total,
    n_above        = sum(sorted_w >= cutoff),
    n_below        = sum(sorted_w <  cutoff),
    sorted_weights = sorted_w,
    signal_fit     = signal_fit,
    noise_fit      = noise_fit,
    crossover_rank = crossover_rank
  )
}

plot_weight_cutoff <- function(cutoff_result, ancestry, out_dir) {
  sorted_w   <- cutoff_result$sorted_weights
  ranks      <- seq_along(sorted_w)
  cutoff     <- cutoff_result$cutoff
  cross_rank <- cutoff_result$crossover_rank
  signal_fit <- cutoff_result$signal_fit
  noise_fit  <- cutoff_result$noise_fit

  out_path <- file.path(out_dir, sprintf("weight_cutoff_%s.png", tolower(ancestry)))

  png(out_path, width = 800, height = 600, res = 150)
  plot(ranks, sorted_w, pch = 16, cex = 0.5,
       xlab = "Rank (descending weight)", ylab = "Raw W_k value",
       main = sprintf("bNMF Weight Cutoff - %s (cutoff = %.4f)", ancestry, cutoff))
  if (!is.null(signal_fit)) abline(signal_fit, col = "blue", lwd = 2, lty = 2)
  if (!is.null(noise_fit))  abline(noise_fit, col = "red", lwd = 2, lty = 2)
  abline(h = cutoff, col = "darkgreen", lwd = 2, lty = 3)
  if (!is.na(cross_rank)) {
    points(cross_rank, sorted_w[cross_rank], pch = 4, cex = 2, col = "darkgreen", lwd = 2)
  }
  legend("topright",
    legend = c("Signal line (top 1%)", "Noise line (bottom 80%)",
               sprintf("Cutoff = %.4f", cutoff)),
    col = c("blue", "red", "darkgreen"),
    lty = c(2, 2, 3), lwd = 2, cex = 0.6)
  dev.off()

  cat(sprintf("  Cutoff diagnostic plot saved: %s\n", out_path))
}

# ---------------------------------------------------------------------------
# 1. Load W_matrix + filtered_variants, join on VAR_ID
# ---------------------------------------------------------------------------
w_mat <- fread(w_path)
fv    <- fread(fv_path)

cluster_cols <- setdiff(names(w_mat), "VAR_ID")
cat(sprintf("W_matrix:           %d variants x %d clusters (%s)\n",
            nrow(w_mat), length(cluster_cols), paste(cluster_cols, collapse = ", ")))
cat(sprintf("filtered_variants:  %d variants\n", nrow(fv)))

dt <- merge(w_mat, fv[, .(VAR_ID, Effect_Allele, BETA)], by = "VAR_ID")
cat(sprintf("After join:         %d variants\n\n", nrow(dt)))

# Show BETA distribution (positive vs negative)
n_pos <- sum(dt$BETA > 0)
n_neg <- sum(dt$BETA < 0)
cat(sprintf("BETA distribution:  %d positive, %d negative (%.0f%% negative)\n",
            n_pos, n_neg, 100 * n_neg / nrow(dt)))
cat("  — All should be included in the score file after the abs() fix.\n\n")

# ---------------------------------------------------------------------------
# 2. LiftOver GRCh37 → GRCh38 (cached)
# ---------------------------------------------------------------------------
cache_path <- file.path(out_dir, "score_files", sprintf("lifted_%s.tsv", tolower(ANC)))

if (file.exists(cache_path)) {
  cached <- fread(cache_path, nrows = 1)
  if ("VAR_ID_38_alt" %in% names(cached)) {
    cat("Loading cached liftover...\n")
    dt <- fread(cache_path)
  } else {
    file.remove(cache_path)
  }
}

if (!"VAR_ID_38" %in% names(dt)) {
  cat("Running liftOver GRCh37 → GRCh38...\n")

  dt[, c("CHR", "POS", "A1", "A2") := {
    parts <- tstrsplit(VAR_ID, "_", fixed = TRUE)
    list(parts[[1]], as.integer(parts[[2]]), parts[[3]], parts[[4]])
  }]

  bed_in       <- file.path(out_dir, "tmp", "input.bed")
  bed_out      <- file.path(out_dir, "tmp", "output.bed")
  bed_unmapped <- file.path(out_dir, "tmp", "unmapped.bed")

  bed_dt <- dt[, .(chrom = paste0("chr", CHR), start = POS - 1L, end = POS, name = VAR_ID)]
  fwrite(bed_dt, bed_in, sep = "\t", col.names = FALSE)

  cmd <- sprintf("%s %s %s %s %s",
                 shQuote(liftover_path), shQuote(bed_in), shQuote(chain_path),
                 shQuote(bed_out), shQuote(bed_unmapped))
  exit_code <- system(cmd, ignore.stdout = TRUE)
  if (exit_code != 0) stop("liftOver failed with exit code ", exit_code)

  lifted <- fread(bed_out, header = FALSE,
                  col.names = c("chrom_38", "start_38", "end_38", "VAR_ID"))
  lifted[, POS_38 := end_38]
  lifted[, CHR_38 := sub("^chr", "", chrom_38)]

  cat(sprintf("  Lifted: %d variants, unmapped: %d\n", nrow(lifted), nrow(dt) - nrow(lifted)))
  dt <- merge(dt, lifted[, .(VAR_ID, CHR_38, POS_38)], by = "VAR_ID", all.x = FALSE)

  dt[, VAR_ID_38     := paste0("chr", CHR_38, ":", POS_38, ":", A1, ":", A2)]
  dt[, VAR_ID_38_alt := paste0("chr", CHR_38, ":", POS_38, ":", A2, ":", A1)]

  fwrite(dt, cache_path, sep = "\t")
  file.remove(bed_in, bed_out)
  if (file.exists(bed_unmapped)) file.remove(bed_unmapped)
}

cat(sprintf("Variants after liftover: %d\n\n", nrow(dt)))

# ---------------------------------------------------------------------------
# 3. Apply bNMF weight cutoff (on raw W_k, before BETA multiplication)
# ---------------------------------------------------------------------------
# Operate on a copy so the liftover cache stays unmodified
dt_scored <- copy(dt)

if (APPLY_WEIGHT_CUTOFF) {
  cutoff_result <- compute_weight_cutoff(dt_scored, cluster_cols)
  cutoff_val <- cutoff_result$cutoff

  if (cutoff_val > 0) {
    cat(sprintf("Weight cutoff: %.6f (retaining %d / %d non-zero weights)\n",
      cutoff_val, cutoff_result$n_above, cutoff_result$n_total))

    for (k in cluster_cols) {
      n_before <- sum(dt_scored[[k]] > 0)
      dt_scored[get(k) < cutoff_val & get(k) > 0, (k) := 0]
      n_after <- sum(dt_scored[[k]] > 0)
      cat(sprintf("  %s: %d -> %d variants (removed %d)\n",
        k, n_before, n_after, n_before - n_after))
    }

    plot_weight_cutoff(cutoff_result, ANC, file.path(out_dir, "diagnostics"))
  } else {
    cat("Weight cutoff: not applied (cutoff = 0 or insufficient data)\n")
  }
  cat("\n")
} else {
  cat("Weight cutoff: DISABLED (--no-weight-cutoff flag)\n\n")
}

# ---------------------------------------------------------------------------
# 4. Weight clusters by BETA + create score file
# ---------------------------------------------------------------------------
for (k in cluster_cols) {
  dt_scored[, (k) := get(k) * BETA]
}
dt_scored[, grs_total := BETA]

all_score_cols <- c(cluster_cols, "grs_total")

# Write score file — with the abs() fix
keep <- rowSums(abs(dt_scored[, all_score_cols, with = FALSE])) > 0
dt_keep <- dt_scored[keep]
cat(sprintf("Variants in score file: %d (of %d — %d excluded for all-zero weights)\n",
            nrow(dt_keep), nrow(dt_scored), nrow(dt_scored) - nrow(dt_keep)))

primary <- dt_keep[, c("VAR_ID_38", "Effect_Allele", all_score_cols), with = FALSE]
setnames(primary, c("VAR_ID_38", "Effect_Allele"), c("ID", "A1"))

alt <- dt_keep[, c("VAR_ID_38_alt", "Effect_Allele", all_score_cols), with = FALSE]
setnames(alt, c("VAR_ID_38_alt", "Effect_Allele"), c("ID", "A1"))

score_dt <- rbind(primary, alt)
score_path <- file.path(out_dir, "score_files", sprintf("score_%s_test.tsv", tolower(ANC)))
fwrite(score_dt, score_path, sep = "\t")
cat(sprintf("Score file written: %s (%d rows)\n\n", basename(score_path), nrow(score_dt)))

# ---------------------------------------------------------------------------
# 5. Run plink2 --score
# ---------------------------------------------------------------------------
cat(sprintf("Running plink2 --score for chr %d-%d...\n", min(TEST_CHRS), max(TEST_CHRS)))

col_nums <- paste(seq(3, 3 + length(all_score_cols) - 1), collapse = ",")
sscore_files <- c()

for (chr in TEST_CHRS) {
  bfile   <- file.path(geno_dir, paste0(geno_prefix, chr))
  chr_out <- file.path(out_dir, "output", sprintf("test_chr%d", chr))

  if (!file.exists(paste0(bfile, ".bed"))) {
    cat(sprintf("  WARNING: %s.bed not found, skipping chr%d\n", bfile, chr))
    next
  }

  cmd <- sprintf(
    "%s --bfile %s --score %s 1 2 header cols=+scoresums ignore-dup-ids --score-col-nums %s --out %s 2>&1",
    shQuote(plink2_path), shQuote(bfile), shQuote(score_path),
    col_nums, shQuote(chr_out)
  )
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

  sscore_path <- paste0(chr_out, ".sscore")
  if (file.exists(sscore_path)) {
    sscore_files <- c(sscore_files, sscore_path)
    cat(sprintf("  chr%d done\n", chr))
  } else {
    cat(sprintf("  chr%d — no sscore output\n", chr))
  }
}

cat(sprintf("\nCompleted %d / %d chromosomes\n\n", length(sscore_files), length(TEST_CHRS)))

if (length(sscore_files) == 0) stop("No sscore files produced. Check plink2 + genotype paths.")

# ---------------------------------------------------------------------------
# 6. Merge per-chromosome results
# ---------------------------------------------------------------------------
combined <- fread(sscore_files[1])
score_sum_cols <- grep("SCORE.*SUM", names(combined), value = TRUE)
id_col <- ifelse("#IID" %in% names(combined), "#IID",
          ifelse("IID" %in% names(combined), "IID", names(combined)[1]))

if (length(sscore_files) > 1) {
  for (f in sscore_files[-1]) {
    chr_dt <- fread(f)
    for (sc in score_sum_cols) {
      if (sc %in% names(chr_dt)) {
        combined[, (sc) := get(sc) + chr_dt[[sc]][match(combined[[id_col]], chr_dt[[id_col]])]]
      }
    }
  }
}

# Rename to friendly names
prs_names <- c(
  sprintf("prs_%s_%s", tolower(cluster_cols), tolower(ANC)),
  sprintf("grs_total_%s", tolower(ANC))
)

if (length(score_sum_cols) != length(prs_names)) {
  stop(sprintf("Column count mismatch: %d score columns vs %d expected names",
               length(score_sum_cols), length(prs_names)))
}

prs_dt <- combined[, c(id_col, score_sum_cols), with = FALSE]
setnames(prs_dt, c("person_id", prs_names))

cat(sprintf("PRS computed for %d individuals\n", nrow(prs_dt)))
cat(sprintf("Score columns: %s\n\n", paste(prs_names, collapse = ", ")))

# ---------------------------------------------------------------------------
# 7. Summary stats (raw, before standardization)
# ---------------------------------------------------------------------------
cat("=== RAW PRS Summary (before standardization) ===\n\n")
for (col in prs_names) {
  vals <- prs_dt[[col]]
  vals <- vals[!is.na(vals)]
  cat(sprintf("--- %s ---\n", col))
  print(summary(vals))
  cat(sprintf("  SD: %.6f   Skewness: %.4f\n\n",
              sd(vals),
              mean(((vals - mean(vals)) / sd(vals))^3)))
}

# ---------------------------------------------------------------------------
# 8. Standardize (z-score)
# ---------------------------------------------------------------------------
prs_z <- copy(prs_dt)
for (col in prs_names) {
  vals <- prs_z[[col]]
  prs_z[, (col) := (vals - mean(vals, na.rm = TRUE)) / sd(vals, na.rm = TRUE)]
}

cat("=== STANDARDIZED PRS Summary ===\n\n")
for (col in prs_names) {
  vals <- prs_z[[col]]
  vals <- vals[!is.na(vals)]
  cat(sprintf("--- %s (z-scored) ---\n", col))
  print(summary(vals))
  cat(sprintf("  SD: %.4f   Skewness: %.4f\n", sd(vals),
              mean(((vals - mean(vals)) / sd(vals))^3)))

  # Shapiro-Wilk on a subsample (max 5000 for the test)
  n_test <- min(length(vals), 5000)
  set.seed(42)
  sw <- shapiro.test(sample(vals, n_test))
  cat(sprintf("  Shapiro-Wilk (n=%d): W = %.6f, p = %.4g\n\n", n_test, sw$statistic, sw$p.value))
}

# ---------------------------------------------------------------------------
# 9. Histogram (density overlay with normal curve)
# ---------------------------------------------------------------------------
plot_df <- melt(prs_z, id.vars = "person_id", measure.vars = prs_names,
                variable.name = "score", value.name = "z")

chr_label <- if (USE_ALL_CHR) "all chr" else sprintf("chr%d-%d", min(TEST_CHRS), max(TEST_CHRS))
cutoff_label <- if (APPLY_WEIGHT_CUTOFF) "with cutoff" else "no cutoff"

p_hist <- ggplot(plot_df, aes(x = z)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80,
                 fill = "steelblue", color = "white", alpha = 0.7) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", linewidth = 0.8, linetype = "dashed") +
  facet_wrap(~ score, scales = "free_y") +
  labs(
    title = sprintf("PRS Normality Check — %s %s (z-scored, %s)",
                    ANC, chr_label, cutoff_label),
    subtitle = "Red dashed line = standard normal N(0,1)",
    x = "Standardized PRS (z-score)",
    y = "Density"
  ) +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(face = "bold"))

hist_path <- file.path(out_dir, "prs_normality_histograms.png")
ggsave(hist_path, p_hist, width = 10, height = 6, dpi = 150)
cat(sprintf("Histogram saved: %s\n", hist_path))

# ---------------------------------------------------------------------------
# 10. QQ plot (observed vs theoretical normal quantiles)
# ---------------------------------------------------------------------------
qq_list <- lapply(prs_names, function(col) {
  vals <- sort(prs_z[[col]][!is.na(prs_z[[col]])])
  n <- length(vals)
  data.table(
    score       = col,
    theoretical = qnorm(ppoints(n)),
    observed    = vals
  )
})
qq_df <- rbindlist(qq_list)

p_qq <- ggplot(qq_df, aes(x = theoretical, y = observed)) +
  geom_point(alpha = 0.2, size = 0.5, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ score, scales = "free") +
  labs(
    title = sprintf("QQ Plot — %s %s (z-scored, %s)",
                    ANC, chr_label, cutoff_label),
    subtitle = "Points should follow the red line if normally distributed",
    x = "Theoretical Quantiles (Normal)",
    y = "Observed Quantiles"
  ) +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(face = "bold"))

qq_path <- file.path(out_dir, "prs_normality_qq.png")
ggsave(qq_path, p_qq, width = 10, height = 6, dpi = 150)
cat(sprintf("QQ plot saved:   %s\n", qq_path))

# ---------------------------------------------------------------------------
# Cleanup temp sscore files
# ---------------------------------------------------------------------------
for (f in sscore_files) {
  file.remove(f)
  log_f <- sub("\\.sscore$", ".log", f)
  if (file.exists(log_f)) file.remove(log_f)
}

cat("\n=============================================\n")
cat("Test complete! Check the PNG files.\n")
cat("=============================================\n")
