#!/usr/bin/env Rscript
# 08_bnmf_vs_hclust_comparison.R
# Compare bNMF soft clusters to hierarchical hard clusters via logistic regression.
# Produces a heatmap: rows = bNMF clusters, cols = hierarchical clusters,
# intensity = -log10(p), color = sign(beta).
#
# For singleton/rare-event hierarchical clusters, uses Firth's penalized
# logistic regression (implemented inline to avoid external dependencies).
#
# Usage:
#   Rscript scripts/a1_analysis/08_bnmf_vs_hclust_comparison.R \
#     --config config/a1_config.yaml --ancestry META

library(data.table)
library(ggplot2)
library(yaml)

# --- Parse CLI args ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/a1_config.yaml"
ancestry    <- "META"
if ("--config"   %in% args) config_path <- args[which(args == "--config") + 1]
if ("--ancestry" %in% args) ancestry    <- args[which(args == "--ancestry") + 1]

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir
figures_dir <- file.path(results_dir, "figures")
sens_dir    <- file.path(results_dir, ancestry, "sensitivity")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("=== bNMF vs Hierarchical Comparison: %s ===\n", ancestry))

# --- Load data ---
w_dt <- fread(file.path(results_dir, ancestry, sprintf("W_matrix_%s.tsv", ancestry)))
hc_dt <- fread(file.path(sens_dir, sprintf("hclust_dynamictreecut_assignments_%s.tsv", ancestry)))

# Merge on VAR_ID
merged <- merge(w_dt, hc_dt, by = "VAR_ID")
setnames(merged, "cluster", "hclust_cluster")
bnmf_cols <- setdiff(colnames(w_dt), "VAR_ID")
n_bnmf <- length(bnmf_cols)
n_hclust <- length(unique(merged$hclust_cluster))

cat(sprintf("Variants: %d | bNMF clusters: %d | Hierarchical clusters: %d\n",
            nrow(merged), n_bnmf, n_hclust))
cat(sprintf("Hierarchical cluster sizes: %s\n",
            paste(sprintf("K%d=%d", seq_len(n_hclust),
                          table(merged$hclust_cluster)), collapse = ", ")))

# --- Firth's penalized logistic regression (single-predictor) ---
# Uses iteratively reweighted least squares with Jeffreys prior penalty.
# Returns list(beta, se, p_value) using penalized likelihood ratio test.
firth_logistic <- function(y, x, max_iter = 100, tol = 1e-6) {
  n <- length(y)
  # Initialize with standard glm if possible, else zeros
  beta <- tryCatch({
    fit0 <- glm(y ~ x, family = binomial)
    unname(coef(fit0))
  }, warning = function(w) c(0, 0),
     error = function(e) c(0, 0))

  for (iter in seq_len(max_iter)) {
    eta <- beta[1] + beta[2] * x
    mu <- 1 / (1 + exp(-eta))
    mu <- pmin(pmax(mu, 1e-10), 1 - 1e-10)
    w <- mu * (1 - mu)

    # Design matrix and hat matrix diagonal
    X <- cbind(1, x)
    W <- diag(as.numeric(w))
    XWX <- t(X) %*% W %*% X
    XWX_inv <- tryCatch(solve(XWX), error = function(e) {
      solve(XWX + diag(1e-6, 2))
    })
    H <- sqrt(W) %*% X %*% XWX_inv %*% t(X) %*% sqrt(W)
    h <- diag(H)

    # Firth-corrected score: U* = X'(y - mu + h*(0.5 - mu))
    y_star <- y - mu + h * (0.5 - mu)
    score <- t(X) %*% y_star

    # Update
    delta <- XWX_inv %*% score
    beta <- beta + as.numeric(delta)

    if (max(abs(delta)) < tol) break
  }

  # Penalized log-likelihood at MLE
  eta_mle <- beta[1] + beta[2] * x
  mu_mle <- 1 / (1 + exp(-eta_mle))
  mu_mle <- pmin(pmax(mu_mle, 1e-10), 1 - 1e-10)
  X <- cbind(1, x)
  W_mle <- diag(as.numeric(mu_mle * (1 - mu_mle)))
  info_mle <- t(X) %*% W_mle %*% X
  log_det_mle <- 0.5 * log(max(det(info_mle), 1e-300))
  pl_full <- sum(y * log(mu_mle) + (1 - y) * log(1 - mu_mle)) + log_det_mle

  # Null model (intercept only) via Firth
  beta0 <- log(mean(y) / (1 - mean(y)))
  mu0 <- rep(mean(y), n)
  mu0 <- pmin(pmax(mu0, 1e-10), 1 - 1e-10)
  X0 <- matrix(1, nrow = n, ncol = 1)
  W0 <- diag(as.numeric(mu0 * (1 - mu0)))
  info0 <- t(X0) %*% W0 %*% X0
  log_det_0 <- 0.5 * log(max(as.numeric(info0), 1e-300))
  pl_null <- sum(y * log(mu0) + (1 - y) * log(1 - mu0)) + log_det_0

  # Penalized LRT
  plrt <- 2 * (pl_full - pl_null)
  p_value <- pchisq(max(plrt, 0), df = 1, lower.tail = FALSE)

  # SE from penalized information
  se <- sqrt(diag(XWX_inv))

  list(beta = beta[2], intercept = beta[1], se = se[2], p_value = p_value)
}

# --- Run all 81 regressions ---
cat("\n--- Running logistic regressions (Firth-penalized) ---\n")

results_list <- list()
hclust_ids <- sort(unique(merged$hclust_cluster))

for (h_id in hclust_ids) {
  y <- as.integer(merged$hclust_cluster == h_id)
  n_pos <- sum(y)
  for (b_col in bnmf_cols) {
    x <- merged[[b_col]]
    fit <- firth_logistic(y, x)
    results_list[[length(results_list) + 1]] <- data.table(
      hclust = paste0("K", h_id),
      bnmf = b_col,
      n_pos = n_pos,
      beta = fit$beta,
      se = fit$se,
      p_value = fit$p_value
    )
  }
}

results <- rbindlist(results_list)
results[, neg_log10_p := -log10(pmax(p_value, 1e-300))]
results[, sign_beta := sign(beta)]

cat(sprintf("Completed %d models\n", nrow(results)))

# Save results table
fwrite(results, file.path(sens_dir, sprintf("logistic_bnmf_vs_hclust_%s.tsv", ancestry)),
       sep = "\t")
cat(sprintf("Results table: %s\n",
            file.path(sens_dir, sprintf("logistic_bnmf_vs_hclust_%s.tsv", ancestry))))

# --- Cluster labels ---
# Canonical META labels from scripts/a1_analysis/06_ancestry_trait_barplots.R
bnmf_labels <- c(
  K1 = "Lpa", K2 = "Adiponectin", K3 = "Platelet",
  K4 = "SHBG", K5 = "Blood Pressure-Stature", K6 = "Metabolic",
  K7 = "Triglycerides-HDL", K8 = "ALP-LDL",
  K9 = "Obesity", K10 = "Glycemic"
)
hclust_labels <- setNames(paste0("hc", seq_along(hclust_ids)), paste0("K", hclust_ids))

results[, bnmf_label := bnmf_labels[bnmf]]
results[, hclust_label := hclust_labels[hclust]]

# Cap -log10(p) for visualization
cap_val <- 20
results[, neg_log10_p_cap := pmin(neg_log10_p, cap_val)]

# Signed intensity: positive beta = positive value, negative beta = negative
results[, signed_intensity := neg_log10_p_cap * sign_beta]

# Order hclust columns so the strongest positive association per bNMF row
# appears leftmost, producing an approximate diagonal
bnmf_rank_map <- setNames(seq_along(bnmf_labels), bnmf_labels)
peak_per_hc <- results[signed_intensity > 0,
  .SD[which.max(signed_intensity)], by = hclust_label]
peak_per_hc[, bnmf_rank := bnmf_rank_map[bnmf_label]]
hclust_col_order <- peak_per_hc[order(bnmf_rank, -signed_intensity), hclust_label]
remaining <- setdiff(unique(results$hclust_label), hclust_col_order)
hclust_col_order <- c(hclust_col_order, remaining)

results[, bnmf_label := factor(bnmf_label, levels = bnmf_labels)]
results[, hclust_label := factor(hclust_label, levels = hclust_col_order)]

# Significance annotations
results[, sig := ifelse(p_value < 0.001, "***",
                 ifelse(p_value < 0.01, "**",
                 ifelse(p_value < 0.05, "*", "")))]

# --- Heatmap ---
cat("\n--- Generating heatmap ---\n")

p <- ggplot(results, aes(x = hclust_label, y = bnmf_label, fill = signed_intensity)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig), size = 4, vjust = 0.5) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-cap_val, cap_val),
    name = expression(paste("Signed ", -log[10], "(p)"))
  ) +
  labs(
    title = sprintf("bNMF vs Hierarchical Clustering — %s", ancestry),
    subtitle = "Firth logistic regression: hierarchical membership ~ bNMF weight",
    x = "Hierarchical Cluster",
    y = "bNMF Cluster"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    panel.grid = element_blank()
  )

ggsave(file.path(figures_dir, sprintf("heatmap_bnmf_vs_hclust_%s.png", ancestry)),
       p, width = 10, height = 8, dpi = 300)
cat(sprintf("Saved: %s\n", file.path(figures_dir,
            sprintf("heatmap_bnmf_vs_hclust_%s.png", ancestry))))

# --- Print summary ---
cat("\n--- Top associations (p < 0.05) ---\n")
sig_results <- results[p_value < 0.05][order(p_value)]
if (nrow(sig_results) > 0) {
  for (i in seq_len(min(20, nrow(sig_results)))) {
    r <- sig_results[i]
    cat(sprintf("  bNMF %s -> HClust %s: beta=%.3f, p=%.2e\n",
                r$bnmf_label, r$hclust_label, r$beta, r$p_value))
  }
}

cat("\nDone.\n")
