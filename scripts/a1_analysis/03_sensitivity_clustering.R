#!/usr/bin/env Rscript
# 03_sensitivity_clustering.R
# Sensitivity analysis: compare bNMF clusters to k-means and hierarchical
# clustering on the same Z-score matrix.
#
# Usage:
#   Rscript scripts/a1_analysis/03_sensitivity_clustering.R \
#     --config config/a1_config.yaml --ancestry EUR

library(data.table)
library(yaml)
library(cluster)          # silhouette
library(dynamicTreeCut)   # cutreeHybrid
library(ggplot2)

# --- Parse CLI args ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/a1_config.yaml"
ancestry    <- "EUR"
if ("--config"   %in% args) config_path <- args[which(args == "--config") + 1]
if ("--ancestry" %in% args) ancestry    <- args[which(args == "--ancestry") + 1]

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir

cat(sprintf("=== Sensitivity Clustering Analysis: %s ===\n", ancestry))

# --- Source figure utilities ---
project_root <- getwd()
source(file.path(project_root, "scripts", "a1_analysis", "figure_utils.R"))

# --- Load prepared matrix (non-negative expanded) ---
prepared_path <- file.path(results_dir, ancestry,
                           sprintf("prepared_matrix_%s.tsv", ancestry))
if (!file.exists(prepared_path)) {
  stop(sprintf("Prepared matrix not found: %s", prepared_path))
}
prep_dt <- fread(prepared_path)
var_ids <- prep_dt$VAR_ID

# Reconstruct signed Z-score matrix: Z = trait_pos - trait_neg
trait_cols <- setdiff(colnames(prep_dt), "VAR_ID")
base_traits <- unique(sub("_(pos|neg)$", "", trait_cols))

z_mat <- matrix(0, nrow = nrow(prep_dt), ncol = length(base_traits),
                dimnames = list(var_ids, base_traits))
for (tr in base_traits) {
  pc <- paste0(tr, "_pos")
  nc <- paste0(tr, "_neg")
  pv <- if (pc %in% colnames(prep_dt)) prep_dt[[pc]] else 0
  nv <- if (nc %in% colnames(prep_dt)) prep_dt[[nc]] else 0
  z_mat[, tr] <- pv - nv
}

cat(sprintf("Z-score matrix: %d variants x %d traits\n", nrow(z_mat), ncol(z_mat)))

# --- Load bNMF W matrix for comparison ---
w_path <- file.path(results_dir, ancestry, sprintf("W_matrix_%s.tsv", ancestry))
bnmf_assignments <- NULL
if (file.exists(w_path)) {
  w_dt <- fread(w_path)
  w_var_ids <- w_dt$VAR_ID
  w_mat <- as.matrix(w_dt[, -"VAR_ID"])
  # Assign each variant to its max-weight cluster
  bnmf_clusters <- apply(w_mat, 1, which.max)
  bnmf_assignments <- data.table(VAR_ID = w_var_ids, cluster = bnmf_clusters)
  cat(sprintf("bNMF: K=%d clusters loaded from W matrix\n", ncol(w_mat)))
}

# --- Output directory ---
out_dir <- file.path(results_dir, ancestry, "sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- K range ---
K_range <- 2:10
# Cap at number of observations - 1
K_range <- K_range[K_range < nrow(z_mat)]

cat(sprintf("Sweeping K = %d to %d\n", min(K_range), max(K_range)))

# ---------------------------------------------------------------------------
# K-means
# ---------------------------------------------------------------------------
cat("\n--- K-means ---\n")

kmeans_results <- list()
kmeans_sil <- numeric(length(K_range))

set.seed(42)
for (i in seq_along(K_range)) {
  k <- K_range[i]
  km <- kmeans(z_mat, centers = k, nstart = 25, iter.max = 100)
  kmeans_results[[as.character(k)]] <- km

  sil <- silhouette(km$cluster, dist(z_mat))
  kmeans_sil[i] <- mean(sil[, "sil_width"])
  cat(sprintf("  K=%d: avg silhouette = %.4f\n", k, kmeans_sil[i]))
}

kmeans_optimal_k <- K_range[which.max(kmeans_sil)]
cat(sprintf("K-means optimal K = %d (silhouette = %.4f)\n",
            kmeans_optimal_k, max(kmeans_sil)))

# Save k-means results for optimal K
km_best <- kmeans_results[[as.character(kmeans_optimal_k)]]
kmeans_assign_dt <- data.table(VAR_ID = var_ids, cluster = km_best$cluster)
fwrite(kmeans_assign_dt,
       file.path(out_dir, sprintf("kmeans_assignments_%s.tsv", ancestry)),
       sep = "\t")

# Centroids = trait profiles
kmeans_centroids <- as.data.table(km_best$centers)
kmeans_centroids[, Cluster := paste0("K", seq_len(.N))]
setcolorder(kmeans_centroids, c("Cluster", base_traits))
fwrite(kmeans_centroids,
       file.path(out_dir, sprintf("kmeans_centroids_%s.tsv", ancestry)),
       sep = "\t")

# ---------------------------------------------------------------------------
# Hierarchical clustering (Ward's method)
# ---------------------------------------------------------------------------
cat("\n--- Hierarchical clustering (Ward.D2) ---\n")

dist_mat <- dist(z_mat, method = "euclidean")
hc <- hclust(dist_mat, method = "ward.D2")

hclust_sil <- numeric(length(K_range))
for (i in seq_along(K_range)) {
  k <- K_range[i]
  cuts <- cutree(hc, k = k)
  sil <- silhouette(cuts, dist_mat)
  hclust_sil[i] <- mean(sil[, "sil_width"])
  cat(sprintf("  K=%d: avg silhouette = %.4f\n", k, hclust_sil[i]))
}

hclust_optimal_k <- K_range[which.max(hclust_sil)]
cat(sprintf("Hierarchical optimal K = %d (silhouette = %.4f)\n",
            hclust_optimal_k, max(hclust_sil)))

# Half-max-height K selection: cut dendrogram at half the max merge height
cut_h <- max(hc$height) / 2
halfmax_cuts <- cutree(hc, h = cut_h)
halfmax_k <- length(unique(halfmax_cuts))
halfmax_sil_obj <- silhouette(halfmax_cuts, dist_mat)
halfmax_sil_avg <- mean(halfmax_sil_obj[, "sil_width"])
cat(sprintf("Half-max-height: cut at h=%.2f -> K=%d (avg silhouette = %.4f)\n",
            cut_h, halfmax_k, halfmax_sil_avg))

# Save half-max-height results
halfmax_assign_dt <- data.table(VAR_ID = var_ids, cluster = halfmax_cuts)
fwrite(halfmax_assign_dt,
       file.path(out_dir, sprintf("hclust_halfmax_assignments_%s.tsv", ancestry)),
       sep = "\t")

halfmax_profiles <- do.call(rbind, lapply(seq_len(halfmax_k), function(k) {
  idx <- which(halfmax_cuts == k)
  colMeans(z_mat[idx, , drop = FALSE])
}))
halfmax_profiles_dt <- as.data.table(halfmax_profiles)
halfmax_profiles_dt[, Cluster := paste0("K", seq_len(.N))]
setcolorder(halfmax_profiles_dt, c("Cluster", base_traits))
fwrite(halfmax_profiles_dt,
       file.path(out_dir, sprintf("hclust_halfmax_profiles_%s.tsv", ancestry)),
       sep = "\t")

# ---------------------------------------------------------------------------
# Dynamic Tree Cut (hybrid method)
# ---------------------------------------------------------------------------
cat("\n--- Dynamic Tree Cut (hybrid) ---\n")

dist_matrix <- as.matrix(dist_mat)

# Sweep deepSplit 0-4
cat("deepSplit sweep (minClusterSize = 5):\n")
dtc_results <- list()
for (ds in 0:4) {
  dtc <- cutreeHybrid(dendro = hc, distM = dist_matrix,
                      minClusterSize = 5, deepSplit = ds,
                      pamStage = TRUE, verbose = 0)
  labels <- dtc$labels
  n_unassigned <- sum(labels == 0)
  n_clusters <- length(unique(labels[labels > 0]))
  sizes <- sort(table(labels[labels > 0]), decreasing = TRUE)
  cat(sprintf("  deepSplit=%d: K=%d clusters, %d unassigned, sizes: %s\n",
              ds, n_clusters, n_unassigned,
              paste(sprintf("%d", sizes), collapse = ", ")))
  dtc_results[[as.character(ds)]] <- list(labels = labels, k = n_clusters,
                                          n_unassigned = n_unassigned)
}

# deepSplit=0 recovers K=10 matching bNMF; higher values over-split
dtc_default_ds <- 0
dtc_labels <- dtc_results[[as.character(dtc_default_ds)]]$labels
dtc_k <- dtc_results[[as.character(dtc_default_ds)]]$k
cat(sprintf("\nUsing deepSplit=%d: K=%d clusters\n", dtc_default_ds, dtc_k))

# Reassign cluster-0 (unassigned) variants to nearest cluster centroid
if (any(dtc_labels == 0)) {
  assigned_idx <- which(dtc_labels > 0)
  unassigned_idx <- which(dtc_labels == 0)
  centroids <- do.call(rbind, lapply(seq_len(dtc_k), function(k) {
    colMeans(z_mat[dtc_labels == k, , drop = FALSE])
  }))
  for (ui in unassigned_idx) {
    dists <- apply(centroids, 1, function(ctr) sum((z_mat[ui, ] - ctr)^2))
    dtc_labels[ui] <- which.min(dists)
  }
  cat(sprintf("Reassigned %d unassigned variants to nearest cluster\n",
              length(unassigned_idx)))
}

# Silhouette for Dynamic Tree Cut
dtc_sil_obj <- silhouette(dtc_labels, dist_mat)
dtc_sil_avg <- mean(dtc_sil_obj[, "sil_width"])
cat(sprintf("Dynamic Tree Cut: K=%d (avg silhouette = %.4f)\n", dtc_k, dtc_sil_avg))

# Save Dynamic Tree Cut assignments
dtc_assign_dt <- data.table(VAR_ID = var_ids, cluster = dtc_labels)
fwrite(dtc_assign_dt,
       file.path(out_dir, sprintf("hclust_dynamictreecut_assignments_%s.tsv", ancestry)),
       sep = "\t")

# Compute trait profiles
dtc_profiles <- do.call(rbind, lapply(seq_len(dtc_k), function(k) {
  idx <- which(dtc_labels == k)
  colMeans(z_mat[idx, , drop = FALSE])
}))
dtc_profiles_dt <- as.data.table(dtc_profiles)
dtc_profiles_dt[, Cluster := paste0("K", seq_len(.N))]
setcolorder(dtc_profiles_dt, c("Cluster", base_traits))
fwrite(dtc_profiles_dt,
       file.path(out_dir, sprintf("hclust_dynamictreecut_profiles_%s.tsv", ancestry)),
       sep = "\t")

# Save hierarchical results for optimal K
hc_cuts <- cutree(hc, k = hclust_optimal_k)
hclust_assign_dt <- data.table(VAR_ID = var_ids, cluster = hc_cuts)
fwrite(hclust_assign_dt,
       file.path(out_dir, sprintf("hclust_assignments_%s.tsv", ancestry)),
       sep = "\t")

# Compute trait profiles (mean Z-score per cluster)
hclust_profiles <- do.call(rbind, lapply(seq_len(hclust_optimal_k), function(k) {
  idx <- which(hc_cuts == k)
  colMeans(z_mat[idx, , drop = FALSE])
}))
hclust_profiles_dt <- as.data.table(hclust_profiles)
hclust_profiles_dt[, Cluster := paste0("K", seq_len(.N))]
setcolorder(hclust_profiles_dt, c("Cluster", base_traits))
fwrite(hclust_profiles_dt,
       file.path(out_dir, sprintf("hclust_profiles_%s.tsv", ancestry)),
       sep = "\t")

# ---------------------------------------------------------------------------
# Save dendrogram
# ---------------------------------------------------------------------------
png(file.path(out_dir, sprintf("dendrogram_%s.png", ancestry)),
    width = 1200, height = 600, res = 150)
plot(hc, labels = FALSE, main = sprintf("Hierarchical Clustering — %s", ancestry),
     xlab = "Variants", sub = "Ward.D2 linkage")
rect.hclust(hc, k = hclust_optimal_k, border = "red")
abline(h = cut_h, col = "blue", lty = 2, lwd = 1.5)
legend("topright",
       legend = c(sprintf("Silhouette K=%d", hclust_optimal_k),
                  sprintf("Half-max-height K=%d (h=%.1f)", halfmax_k, cut_h)),
       col = c("red", "blue"), lty = c(1, 2), lwd = 1.5, cex = 0.8)
dev.off()

# ---------------------------------------------------------------------------
# Silhouette comparison plot
# ---------------------------------------------------------------------------
sil_dt <- data.table(
  K = rep(K_range, 2),
  silhouette = c(kmeans_sil, hclust_sil),
  method = rep(c("K-means", "Hierarchical"), each = length(K_range))
)

p_sil <- ggplot(sil_dt, aes(x = K, y = silhouette, color = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = kmeans_optimal_k, linetype = "dashed",
             color = "#E41A1C", alpha = 0.5) +
  geom_vline(xintercept = hclust_optimal_k, linetype = "dashed",
             color = "#377EB8", alpha = 0.5) +
  scale_x_continuous(breaks = K_range) +
  scale_color_manual(values = c("K-means" = "#E41A1C", "Hierarchical" = "#377EB8")) +
  labs(title = sprintf("K Selection — %s", ancestry),
       x = "Number of clusters (K)",
       y = "Average silhouette width",
       color = "Method") +
  theme_big_text(base_size = 14)

ggsave(file.path(out_dir, sprintf("k_selection_%s.png", ancestry)),
       p_sil, width = 8, height = 5, dpi = 150)

# ---------------------------------------------------------------------------
# Trait-profile heatmaps (comparable to bNMF Panel A)
# ---------------------------------------------------------------------------
plot_trait_heatmap <- function(profiles_dt, method_label, ancestry, out_dir) {
  long <- melt(profiles_dt, id.vars = "Cluster",
               variable.name = "Trait", value.name = "Z_score")

  p <- ggplot(long, aes(x = Trait, y = Cluster, fill = Z_score)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "Mean Z") +
    labs(title = sprintf("%s — %s (K=%d)", method_label, ancestry,
                         nrow(profiles_dt))) +
    theme_big_text(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank())

  fname <- sprintf("heatmap_%s_%s.png",
                    tolower(gsub("[- ]", "", method_label)), ancestry)
  ggsave(file.path(out_dir, fname), p,
         width = max(8, 0.5 * ncol(profiles_dt)), height = 3 + 0.5 * nrow(profiles_dt),
         dpi = 150)
  cat(sprintf("Wrote heatmap: %s\n", fname))
}

plot_trait_heatmap(kmeans_centroids, "K-means", ancestry, out_dir)
plot_trait_heatmap(hclust_profiles_dt, "Hierarchical", ancestry, out_dir)
plot_trait_heatmap(halfmax_profiles_dt, "Hierarchical-HalfMax", ancestry, out_dir)
plot_trait_heatmap(dtc_profiles_dt, "DynamicTreeCut", ancestry, out_dir)

# ---------------------------------------------------------------------------
# Adjusted Rand Index comparison
# ---------------------------------------------------------------------------
cat("\n--- Cross-method comparison (ARI) ---\n")

compute_ari <- function(labels1, labels2) {
  if (requireNamespace("mclust", quietly = TRUE)) {
    return(mclust::adjustedRandIndex(labels1, labels2))
  }
  # Fallback: manual ARI computation
  tab <- table(labels1, labels2)
  n <- sum(tab)
  sum_comb_a <- sum(choose(rowSums(tab), 2))
  sum_comb_b <- sum(choose(colSums(tab), 2))
  sum_comb_ab <- sum(choose(tab, 2))
  expected <- sum_comb_a * sum_comb_b / choose(n, 2)
  max_idx <- 0.5 * (sum_comb_a + sum_comb_b)
  if (max_idx == expected) return(0)
  (sum_comb_ab - expected) / (max_idx - expected)
}

# K-means vs hierarchical (using optimal K for each — need matching K for ARI)
# Recut hierarchical at kmeans_optimal_k for fair comparison
hc_at_km_k <- cutree(hc, k = kmeans_optimal_k)
ari_km_hc <- compute_ari(km_best$cluster, hc_at_km_k)
cat(sprintf("ARI(K-means K=%d, Hierarchical K=%d) = %.4f\n",
            kmeans_optimal_k, kmeans_optimal_k, ari_km_hc))

# Half-max vs silhouette-optimal hierarchical
ari_halfmax_sil <- compute_ari(halfmax_cuts, hc_cuts)
cat(sprintf("ARI(HalfMax K=%d, Silhouette-Hierarchical K=%d) = %.4f\n",
            halfmax_k, hclust_optimal_k, ari_halfmax_sil))

# Compare with bNMF if available
ari_km_bnmf <- NA_real_
ari_hc_bnmf <- NA_real_
ari_halfmax_bnmf <- NA_real_
ari_dtc_bnmf <- NA_real_
bnmf_k <- NA_integer_

if (!is.null(bnmf_assignments)) {
  # Match variant IDs
  common_vars <- intersect(var_ids, bnmf_assignments$VAR_ID)
  if (length(common_vars) > 10) {
    km_sub <- km_best$cluster[match(common_vars, var_ids)]
    bnmf_sub <- bnmf_assignments$cluster[match(common_vars, bnmf_assignments$VAR_ID)]
    bnmf_k <- max(bnmf_sub)

    # Rerun k-means and hclust at bNMF K for fair comparison
    set.seed(42)
    km_at_bnmf_k <- kmeans(z_mat, centers = bnmf_k, nstart = 25, iter.max = 100)
    hc_at_bnmf_k <- cutree(hc, k = bnmf_k)

    km_sub_bnmf_k <- km_at_bnmf_k$cluster[match(common_vars, var_ids)]
    hc_sub_bnmf_k <- hc_at_bnmf_k[match(common_vars, var_ids)]

    ari_km_bnmf <- compute_ari(km_sub_bnmf_k, bnmf_sub)
    ari_hc_bnmf <- compute_ari(hc_sub_bnmf_k, bnmf_sub)

    halfmax_sub <- halfmax_cuts[match(common_vars, var_ids)]
    ari_halfmax_bnmf <- compute_ari(halfmax_sub, bnmf_sub)

    cat(sprintf("ARI(K-means K=%d, bNMF K=%d) = %.4f\n",
                bnmf_k, bnmf_k, ari_km_bnmf))
    cat(sprintf("ARI(Hierarchical K=%d, bNMF K=%d) = %.4f\n",
                bnmf_k, bnmf_k, ari_hc_bnmf))
    cat(sprintf("ARI(HalfMax K=%d, bNMF K=%d) = %.4f\n",
                halfmax_k, bnmf_k, ari_halfmax_bnmf))

    dtc_sub <- dtc_labels[match(common_vars, var_ids)]
    ari_dtc_bnmf <- compute_ari(dtc_sub, bnmf_sub)
    cat(sprintf("ARI(DynamicTreeCut K=%d, bNMF K=%d) = %.4f\n",
                dtc_k, bnmf_k, ari_dtc_bnmf))
  } else {
    ari_dtc_bnmf <- NA_real_
    cat("Too few shared variants for bNMF comparison\n")
  }
}

# ---------------------------------------------------------------------------
# Summary report
# ---------------------------------------------------------------------------
summary_lines <- c(
  sprintf("=== Sensitivity Clustering Summary: %s ===", ancestry),
  sprintf("Date: %s", Sys.time()),
  "",
  sprintf("Input: %d variants x %d traits", nrow(z_mat), ncol(z_mat)),
  sprintf("K range tested: %d-%d", min(K_range), max(K_range)),
  "",
  "--- K-means ---",
  sprintf("Optimal K: %d (avg silhouette = %.4f)", kmeans_optimal_k, max(kmeans_sil)),
  sprintf("Cluster sizes: %s",
          paste(sprintf("K%d=%d", seq_along(table(km_best$cluster)),
                        table(km_best$cluster)), collapse = ", ")),
  "",
  "--- Hierarchical (Ward.D2, silhouette) ---",
  sprintf("Optimal K: %d (avg silhouette = %.4f)", hclust_optimal_k, max(hclust_sil)),
  sprintf("Cluster sizes: %s",
          paste(sprintf("K%d=%d", seq_along(table(hc_cuts)),
                        table(hc_cuts)), collapse = ", ")),
  "",
  "--- Hierarchical (Ward.D2, half-max-height) ---",
  sprintf("Cut height: %.2f (max = %.2f)", cut_h, max(hc$height)),
  sprintf("Resulting K: %d (avg silhouette = %.4f)", halfmax_k, halfmax_sil_avg),
  sprintf("Cluster sizes: %s",
          paste(sprintf("K%d=%d", seq_along(table(halfmax_cuts)),
                        table(halfmax_cuts)), collapse = ", ")),
  "",
  sprintf("--- Dynamic Tree Cut (deepSplit=%d) ---", dtc_default_ds),
  sprintf("Resulting K: %d (avg silhouette = %.4f)", dtc_k, dtc_sil_avg),
  sprintf("Cluster sizes: %s",
          paste(sprintf("K%d=%d", seq_along(table(dtc_labels)),
                        table(dtc_labels)), collapse = ", ")),
  "",
  "--- Cross-method ARI ---",
  sprintf("K-means vs Hierarchical (at K=%d): %.4f", kmeans_optimal_k, ari_km_hc),
  sprintf("HalfMax (K=%d) vs Silhouette-Hierarchical (K=%d): %.4f",
          halfmax_k, hclust_optimal_k, ari_halfmax_sil)
)

if (!is.na(ari_km_bnmf)) {
  summary_lines <- c(summary_lines,
    sprintf("K-means vs bNMF (at K=%d): %.4f", bnmf_k, ari_km_bnmf),
    sprintf("Hierarchical vs bNMF (at K=%d): %.4f", bnmf_k, ari_hc_bnmf),
    sprintf("HalfMax (K=%d) vs bNMF (K=%d): %.4f", halfmax_k, bnmf_k, ari_halfmax_bnmf),
    sprintf("DynamicTreeCut (K=%d) vs bNMF (K=%d): %.4f", dtc_k, bnmf_k, ari_dtc_bnmf)
  )
}

summary_lines <- c(summary_lines, "",
  "--- Silhouette scores ---",
  sprintf("K\tK-means\tHierarchical"))
for (i in seq_along(K_range)) {
  summary_lines <- c(summary_lines,
    sprintf("%d\t%.4f\t%.4f", K_range[i], kmeans_sil[i], hclust_sil[i]))
}

summary_path <- file.path(out_dir, sprintf("sensitivity_summary_%s.txt", ancestry))
writeLines(summary_lines, summary_path)
cat(sprintf("\nSummary written to: %s\n", summary_path))
cat(sprintf("Output directory: %s\n", out_dir))
