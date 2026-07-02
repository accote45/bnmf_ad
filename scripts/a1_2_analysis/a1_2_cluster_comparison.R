#!/usr/bin/env Rscript
#
# A1.2: Compare META bNMF clusters against published clustering studies
#
# Constructs trait weight profiles (weighted-average z-scores) for each set
# of clusters, then computes Pearson correlation between the user's 10 META
# clusters and comparator clusters from Suzuki 2024, Smith 2024, and
# Pascat 2026.
#
# Usage:
#   Rscript scripts/a1_2_analysis/a1_2_cluster_comparison.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(readxl)
  library(yaml)
  library(cowplot)
})

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"

theme_big_text <- function(base_size = 18, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title    = element_text(size = base_size * 1.2, face = "bold"),
      plot.subtitle = element_text(size = base_size, face = "italic"),
      axis.title    = element_text(size = base_size),
      axis.text     = element_text(size = base_size * 0.9),
      legend.title  = element_text(size = base_size),
      legend.text   = element_text(size = base_size * 0.9),
      strip.text    = element_text(size = base_size),
      panel.grid.minor = element_blank()
    )
}

# -- Config --------------------------------------------------------------------

out_dir       <- file.path(base_dir, "results/a1_2_analysis")
results_dir   <- file.path(base_dir, "results/a1_analysis/META")
published_dir <- file.path(base_dir, "data/published_clusters")
config_path   <- file.path(base_dir, "config/a1_config.yaml")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Canonical META cluster labels (from 07_gene_bars_metric_comparison.R)
meta_labels <- c(
  K1 = "Lpa",                            K2 = "Adiponectin",
  K3 = "Platelet",                       K4 = "SHBG",
  K5 = "Blood Pressure-Stature",         K6 = "Metabolic",
  K7 = "Triglycerides-HDL",              K8 = "ALP-LDL",
  K9 = "Obesity",                        K10 = "Glycemic"
)

message("=== A1.2: Published Cluster Comparison ===\n")

# -- Step 1: Load user's bNMF data --------------------------------------------

message("Step 1: Loading META bNMF results...")

W <- fread(file.path(results_dir, "W_matrix_META.tsv"))
w_mat <- as.matrix(W[, -1, with = FALSE])
rownames(w_mat) <- W$VAR_ID
message(sprintf("  W matrix: %d SNVs x %d clusters", nrow(w_mat), ncol(w_mat)))

prep <- fread(file.path(results_dir, "prepared_matrix_META.tsv"))
prep_mat <- as.matrix(prep[, -1, with = FALSE])
rownames(prep_mat) <- prep$VAR_ID
message(sprintf("  Prepared matrix: %d SNVs x %d columns", nrow(prep_mat), ncol(prep_mat)))

# Convert pos/neg to signed z-scores
pos_cols <- grep("_pos$", colnames(prep_mat), value = TRUE)
trait_names_full <- sub("_pos$", "", pos_cols)

z_signed <- matrix(NA_real_, nrow = nrow(prep_mat), ncol = length(trait_names_full),
                   dimnames = list(rownames(prep_mat), trait_names_full))
for (tr in trait_names_full) {
  z_signed[, tr] <- prep_mat[, paste0(tr, "_pos")] - prep_mat[, paste0(tr, "_neg")]
}
message(sprintf("  Signed z-score matrix: %d SNVs x %d traits", nrow(z_signed), ncol(z_signed)))

# Compute user cluster profiles: weighted average z-score
user_profiles <- matrix(NA_real_, nrow = ncol(w_mat), ncol = ncol(z_signed),
                        dimnames = list(colnames(w_mat), colnames(z_signed)))
for (k in seq_len(ncol(w_mat))) {
  weights <- w_mat[, k]
  wsum <- sum(weights)
  if (wsum > 0) {
    user_profiles[k, ] <- colSums(weights * z_signed) / wsum
  }
}
message(sprintf("  User profiles: %d clusters x %d traits", nrow(user_profiles), ncol(user_profiles)))

# -- Step 2: Trait name harmonization ------------------------------------------

message("Step 2: Harmonizing trait names...")

strip_source <- function(x) sub("_[A-Za-z]+[0-9]{4}$", "", x)

# Deduplicate: for each base trait, keep version with highest mean |profile|
clean_names <- strip_source(colnames(user_profiles))
if (any(duplicated(clean_names))) {
  mean_abs <- colMeans(abs(user_profiles))
  keep_idx <- tapply(seq_along(clean_names), clean_names, function(idx) idx[which.max(mean_abs[idx])])
  keep_idx <- unlist(keep_idx)
  user_profiles_dedup <- user_profiles[, keep_idx, drop = FALSE]
  colnames(user_profiles_dedup) <- clean_names[keep_idx]
} else {
  user_profiles_dedup <- user_profiles
  colnames(user_profiles_dedup) <- clean_names
}
message(sprintf("  Deduplicated user profiles: %d unique base traits", ncol(user_profiles_dedup)))

# Also keep mapping from full trait name to file path (needed for z-score lookup)
trait_name_to_full <- setNames(colnames(user_profiles), strip_source(colnames(user_profiles)))

# -- Step 3: Load comparator study SNV lists -----------------------------------

message("Step 3: Loading comparator studies...")

# Suzuki 2024: hard clustering, 8 clusters
suzuki_raw <- read_excel(file.path(published_dir, "Suzuki2024_Nature_SuppTables.xlsx"),
                         sheet = 6, skip = 2)
suzuki <- data.table(
  rsid    = suzuki_raw[[3]],
  cluster = suzuki_raw[[5]]
)
suzuki <- suzuki[!is.na(rsid) & !is.na(cluster)]
message(sprintf("  Suzuki 2024: %d SNVs, %d clusters", nrow(suzuki), length(unique(suzuki$cluster))))

# Smith 2024: soft clustering, 12 clusters
smith_raw <- read_excel(file.path(published_dir, "Smith2024_NatMed_SuppTables.xlsx"),
                        sheet = 5, skip = 3)
smith_raw <- smith_raw[, 1:15]
smith_ids <- data.table(
  var_id = as.character(smith_raw[[1]]),
  rsid   = as.character(smith_raw[[2]])
)
# Cluster weight columns (4-15)
smith_clusters <- names(smith_raw)[4:15]
smith_clusters_clean <- sub("\\.\\.\\.[0-9]+$", "", smith_clusters)
smith_w <- as.matrix(smith_raw[, 4:15])
smith_w[is.na(smith_w)] <- 0
colnames(smith_w) <- smith_clusters_clean
smith_ids <- smith_ids[!is.na(var_id)]
smith_w <- smith_w[!is.na(smith_raw[[1]]), , drop = FALSE]
message(sprintf("  Smith 2024: %d SNVs, %d clusters (soft)", nrow(smith_ids), ncol(smith_w)))

# Pascat 2026: hard clustering, 5 clusters
pascat_raw <- read_excel(file.path(published_dir, "Pascat2026_NatComms_SuppTables.xlsx"),
                         sheet = 4, skip = 3)
pascat <- data.table(
  rsid    = as.character(pascat_raw[[1]]),
  cluster = as.character(pascat_raw[[4]])
)
pascat <- pascat[!is.na(rsid) & !is.na(cluster)]
message(sprintf("  Pascat 2026: %d SNVs, %d clusters", nrow(pascat), length(unique(pascat$cluster))))

# -- Step 4: Build z-score lookup from harmonized sumstats ---------------------

message("Step 4: Building z-score lookup from harmonized GWAS sumstats...")

# Collect all needed RSIDs and VAR_IDs
all_rsids  <- unique(c(suzuki$rsid, smith_ids$rsid, pascat$rsid))
all_varids <- unique(smith_ids$var_id)
all_rsids  <- all_rsids[!is.na(all_rsids)]
all_varids <- all_varids[!is.na(all_varids)]
message(sprintf("  Total unique RSIDs to lookup: %d", length(all_rsids)))
message(sprintf("  Total unique VAR_IDs to lookup: %d", length(all_varids)))

# Parse config to get trait → sumstat file paths
cfg <- read_yaml(config_path)
fallback <- cfg$trait_gwas$fallback_order
trait_files <- list()
for (anc in fallback) {
  anc_traits <- cfg$trait_gwas[[anc]]
  if (!is.null(anc_traits)) {
    for (tn in names(anc_traits)) {
      if (!(tn %in% names(trait_files))) {
        trait_files[[tn]] <- file.path(base_dir, anc_traits[[tn]])
      }
    }
  }
}

# Only look up traits that are in the prepared matrix
trait_files <- trait_files[names(trait_files) %in% trait_names_full]
message(sprintf("  Trait sumstat files to read: %d", length(trait_files)))

# Read each sumstat file and extract z-scores for needed variants
zscore_list <- list()
for (i in seq_along(trait_files)) {
  tn <- names(trait_files)[i]
  fp <- trait_files[[tn]]

  if (!file.exists(fp)) {
    message(sprintf("    [%d/%d] MISSING: %s", i, length(trait_files), basename(fp)))
    next
  }

  ss <- fread(fp, select = c("VAR_ID", "RSID", "BETA", "SE"))
  ss <- ss[!is.na(BETA) & !is.na(SE) & SE > 0]
  ss[, z := BETA / SE]

  # Filter to needed variants
  ss_filt <- ss[RSID %in% all_rsids | VAR_ID %in% all_varids]

  if (nrow(ss_filt) > 0) {
    zscore_list[[tn]] <- ss_filt[, .(VAR_ID, RSID, z)]
  }

  if (i %% 10 == 0 || i == length(trait_files)) {
    message(sprintf("    [%d/%d] %s: %d variants matched",
                    i, length(trait_files), tn, nrow(ss_filt)))
  }
}

# Build wide z-score lookup tables keyed by RSID and VAR_ID
build_wide_lookup <- function(zscore_list, key_col) {
  long <- rbindlist(lapply(names(zscore_list), function(tn) {
    dt <- zscore_list[[tn]][, .(key = get(key_col), z)]
    dt[, trait := tn]
    dt[!is.na(key)]
  }))
  if (nrow(long) == 0) return(NULL)
  dcast(long, key ~ trait, value.var = "z", fun.aggregate = mean)
}

lookup_rsid  <- build_wide_lookup(zscore_list, "RSID")
lookup_varid <- build_wide_lookup(zscore_list, "VAR_ID")
message(sprintf("  RSID lookup: %d variants x %d traits",
                if (!is.null(lookup_rsid)) nrow(lookup_rsid) else 0,
                if (!is.null(lookup_rsid)) ncol(lookup_rsid) - 1 else 0))
message(sprintf("  VAR_ID lookup: %d variants x %d traits",
                if (!is.null(lookup_varid)) nrow(lookup_varid) else 0,
                if (!is.null(lookup_varid)) ncol(lookup_varid) - 1 else 0))

# -- Step 5: Compute comparator trait profiles ---------------------------------

message("Step 5: Computing comparator trait profiles...")

# Helper: compute hard-clustering profiles (mean z per cluster)
compute_hard_profiles <- function(rsids, clusters, lookup, study_name) {
  matched <- lookup[key %in% rsids]
  n_matched <- nrow(matched)
  n_total   <- length(unique(rsids))
  message(sprintf("  %s: %d / %d variants matched (%.1f%%)",
                  study_name, n_matched, n_total, 100 * n_matched / n_total))

  # Map matched RSIDs back to cluster assignments
  id_cluster <- as.data.table(list(key = rsids, cluster = clusters))
  id_cluster <- id_cluster[key %in% matched$key]

  trait_cols <- setdiff(names(matched), "key")
  z_mat <- as.matrix(matched[, ..trait_cols])
  rownames(z_mat) <- matched$key

  cluster_names <- sort(unique(id_cluster$cluster))
  profiles <- matrix(NA_real_, nrow = length(cluster_names), ncol = length(trait_cols),
                     dimnames = list(cluster_names, trait_cols))
  for (cl in cluster_names) {
    cl_rsids <- id_cluster[cluster == cl]$key
    cl_z <- z_mat[rownames(z_mat) %in% cl_rsids, , drop = FALSE]
    if (nrow(cl_z) > 0) {
      profiles[cl, ] <- colMeans(cl_z, na.rm = TRUE)
    }
  }
  profiles
}

# Helper: compute soft-clustering profiles (weighted average z)
compute_soft_profiles <- function(ids, w_matrix, lookup, id_col, study_name) {
  if (id_col == "VAR_ID") {
    matched <- lookup_varid[key %in% ids]
  } else {
    matched <- lookup_rsid[key %in% ids]
  }
  n_matched <- nrow(matched)
  n_total   <- length(unique(ids))
  message(sprintf("  %s: %d / %d variants matched (%.1f%%)",
                  study_name, n_matched, n_total, 100 * n_matched / n_total))

  # Align w_matrix rows to matched variants
  match_idx <- match(matched$key, ids)
  match_idx <- match_idx[!is.na(match_idx)]
  w_aligned <- w_matrix[match_idx, , drop = FALSE]

  trait_cols <- setdiff(names(matched), "key")
  z_mat <- as.matrix(matched[, ..trait_cols])
  # Ensure row alignment
  z_keys <- matched$key
  w_keys <- ids[match_idx]
  common <- intersect(z_keys, w_keys)
  z_mat <- z_mat[match(common, z_keys), , drop = FALSE]
  w_aligned <- w_aligned[match(common, w_keys), , drop = FALSE]

  cluster_names <- colnames(w_aligned)
  profiles <- matrix(NA_real_, nrow = length(cluster_names), ncol = length(trait_cols),
                     dimnames = list(cluster_names, trait_cols))
  for (k in seq_along(cluster_names)) {
    weights <- w_aligned[, k]
    wsum <- sum(weights, na.rm = TRUE)
    if (wsum > 0) {
      profiles[k, ] <- colSums(weights * z_mat, na.rm = TRUE) / wsum
    }
  }
  profiles
}

# Suzuki: hard clustering
suzuki_profiles <- compute_hard_profiles(
  suzuki$rsid, suzuki$cluster, lookup_rsid, "Suzuki 2024"
)

# Smith: soft clustering via VAR_ID, with rsID fallback
smith_profiles <- compute_soft_profiles(
  smith_ids$var_id, smith_w, lookup_varid, "VAR_ID", "Smith 2024"
)

# Pascat: hard clustering
pascat_profiles <- compute_hard_profiles(
  pascat$rsid, pascat$cluster, lookup_rsid, "Pascat 2026"
)

# -- Step 6: Trait matching and Pearson correlation ----------------------------

message("Step 6: Computing Pearson correlations...")

# Harmonize comparator trait names and find shared traits with user
harmonize_and_correlate <- function(user_prof, comp_prof, study_name) {
  # Strip source suffixes from comparator trait names
  comp_base <- strip_source(colnames(comp_prof))
  colnames(comp_prof) <- comp_base

  # Deduplicate comparator traits (keep highest mean |profile|)
  if (any(duplicated(comp_base))) {
    mean_abs <- colMeans(abs(comp_prof), na.rm = TRUE)
    keep_idx <- tapply(seq_along(comp_base), comp_base, function(idx) idx[which.max(mean_abs[idx])])
    keep_idx <- unlist(keep_idx)
    comp_prof <- comp_prof[, keep_idx, drop = FALSE]
    colnames(comp_prof) <- comp_base[keep_idx]
  }

  shared <- intersect(colnames(user_prof), colnames(comp_prof))
  message(sprintf("  %s: %d shared traits", study_name, length(shared)))

  if (length(shared) < 3) {
    warning(sprintf("  %s: fewer than 3 shared traits — skipping", study_name))
    return(list(cor_mat = NULL, shared_traits = shared))
  }

  u <- user_prof[, shared, drop = FALSE]
  c <- comp_prof[, shared, drop = FALSE]

  # Pearson correlation: each row of u vs each row of c
  cor_mat <- cor(t(u), t(c), use = "pairwise.complete.obs")
  rownames(cor_mat) <- rownames(u)
  colnames(cor_mat) <- rownames(c)

  # P-values
  p_mat <- matrix(NA_real_, nrow = nrow(u), ncol = nrow(c),
                  dimnames = dimnames(cor_mat))
  for (i in seq_len(nrow(u))) {
    for (j in seq_len(nrow(c))) {
      complete <- complete.cases(u[i, ], c[j, ])
      if (sum(complete) >= 3) {
        ct <- cor.test(u[i, complete], c[j, complete])
        p_mat[i, j] <- ct$p.value
      }
    }
  }

  list(cor_mat = cor_mat, p_mat = p_mat, shared_traits = shared)
}

res_suzuki <- harmonize_and_correlate(user_profiles_dedup, suzuki_profiles, "Suzuki 2024")
res_smith  <- harmonize_and_correlate(user_profiles_dedup, smith_profiles, "Smith 2024")
res_pascat <- harmonize_and_correlate(user_profiles_dedup, pascat_profiles, "Pascat 2026")

# -- Step 7: Visualization ----------------------------------------------------

message("Step 7: Generating heatmaps...")

plot_correlation_heatmap <- function(cor_mat, user_labels, study_name, study_short) {
  if (is.null(cor_mat)) return(NULL)

  # Build semantic row labels for user clusters
  row_ids <- rownames(cor_mat)
  row_labs <- ifelse(row_ids %in% names(user_labels),
                     paste0(row_ids, ": ", user_labels[row_ids]),
                     row_ids)

  col_labs <- colnames(cor_mat)

  long <- as.data.table(expand.grid(
    Row = rownames(cor_mat), Col = colnames(cor_mat), stringsAsFactors = FALSE
  ))
  long[, Value := as.vector(cor_mat)]
  long[, Row_label := factor(row_labs[match(Row, rownames(cor_mat))],
                             levels = rev(row_labs))]
  long[, Col_label := factor(col_labs[match(Col, colnames(cor_mat))],
                             levels = col_labs)]

  p <- ggplot(long, aes(x = Col_label, y = Row_label, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", Value)), size = 3.2) +
    scale_fill_gradient2(low = "steelblue4", mid = "white", high = "firebrick3",
                         midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
    labs(title = sprintf("META vs %s", study_name),
         x = study_name, y = "META Clusters") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white"))

  nc <- ncol(cor_mat)
  nr <- nrow(cor_mat)
  ggsave(file.path(out_dir, sprintf("heatmap_pearson_%s.png", study_short)), p,
         width = max(8, nc * 1.4), height = max(7, nr * 1.0), dpi = 600)
  message(sprintf("  Saved heatmap_pearson_%s.png", study_short))
  p
}

p_suzuki <- plot_correlation_heatmap(res_suzuki$cor_mat, meta_labels, "Suzuki 2024", "suzuki2024")
p_smith  <- plot_correlation_heatmap(res_smith$cor_mat,  meta_labels, "Smith 2024",  "smith2024")
p_pascat <- plot_correlation_heatmap(res_pascat$cor_mat, meta_labels, "Pascat 2026", "pascat2026")

# Combined multi-panel figure
panels <- Filter(Negate(is.null), list(p_suzuki, p_smith, p_pascat))
if (length(panels) > 0) {
  # Relative widths proportional to number of comparator clusters
  widths <- sapply(list(res_suzuki$cor_mat, res_smith$cor_mat, res_pascat$cor_mat), function(m) {
    if (is.null(m)) NULL else ncol(m)
  })
  widths <- unlist(widths[!sapply(widths, is.null)])

  p_combined <- plot_grid(plotlist = panels, nrow = 1, rel_widths = widths,
                          labels = "AUTO", label_size = 14)
  ggsave(file.path(out_dir, "combined_heatmap_published_comparison.png"), p_combined,
         width = sum(widths) * 1.4 + 4, height = 10, dpi = 600)
  message("  Saved combined_heatmap_published_comparison.png")
}

# -- Step 8: Summary outputs --------------------------------------------------

message("Step 8: Writing summary outputs...")

# Correlation summary CSV
build_summary <- function(cor_mat, p_mat, study_name, user_labels, shared_traits) {
  if (is.null(cor_mat)) return(NULL)
  long <- as.data.table(expand.grid(
    User_Cluster = rownames(cor_mat),
    Comparator_Cluster = colnames(cor_mat),
    stringsAsFactors = FALSE
  ))
  long[, Pearson_r := as.vector(cor_mat)]
  long[, P_value := as.vector(p_mat)]
  long[, User_Label := user_labels[User_Cluster]]
  long[, Study := study_name]
  long[, N_Shared_Traits := length(shared_traits)]
  long[, .(User_Cluster, User_Label, Study, Comparator_Cluster, Pearson_r, P_value, N_Shared_Traits)]
}

summary_all <- rbindlist(Filter(Negate(is.null), list(
  build_summary(res_suzuki$cor_mat, res_suzuki$p_mat, "Suzuki 2024", meta_labels, res_suzuki$shared_traits),
  build_summary(res_smith$cor_mat,  res_smith$p_mat,  "Smith 2024",  meta_labels, res_smith$shared_traits),
  build_summary(res_pascat$cor_mat, res_pascat$p_mat, "Pascat 2026", meta_labels, res_pascat$shared_traits)
)))
fwrite(summary_all, file.path(out_dir, "correlation_summary.csv"))
message(sprintf("  Saved correlation_summary.csv (%d rows)", nrow(summary_all)))

# Best match per user cluster per study
best_match <- summary_all[, .SD[which.max(abs(Pearson_r))], by = .(User_Cluster, Study)]
best_match[, Rank := 1]
fwrite(best_match, file.path(out_dir, "best_match_summary.csv"))
message(sprintf("  Saved best_match_summary.csv (%d rows)", nrow(best_match)))

# Per-study correlation matrices as TSV
save_cor_tsv <- function(cor_mat, study_short) {
  if (is.null(cor_mat)) return()
  dt <- as.data.table(cor_mat, keep.rownames = "User_Cluster")
  fwrite(dt, file.path(out_dir, sprintf("cor_matrix_%s.tsv", study_short)), sep = "\t")
}
save_cor_tsv(res_suzuki$cor_mat, "suzuki2024")
save_cor_tsv(res_smith$cor_mat,  "smith2024")
save_cor_tsv(res_pascat$cor_mat, "pascat2026")

# Text summary
sink(file.path(out_dir, "comparison_summary.txt"))
cat("=== A1.2: Published Cluster Comparison Summary ===\n\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

cat("--- Comparator Studies ---\n")
cat(sprintf("Suzuki 2024: %d SNVs, %d clusters (hard)\n", nrow(suzuki), length(unique(suzuki$cluster))))
cat(sprintf("Smith 2024:  %d SNVs, %d clusters (soft)\n", nrow(smith_ids), ncol(smith_w)))
cat(sprintf("Pascat 2026: %d SNVs, %d clusters (hard)\n", nrow(pascat), length(unique(pascat$cluster))))

cat("\n--- Shared Traits ---\n")
cat(sprintf("Suzuki: %d\n", length(res_suzuki$shared_traits)))
cat(sprintf("Smith:  %d\n", length(res_smith$shared_traits)))
cat(sprintf("Pascat: %d\n", length(res_pascat$shared_traits)))

cat("\n--- Best Matches ---\n")
for (i in seq_len(nrow(best_match))) {
  r <- best_match[i]
  cat(sprintf("  %s (%s) <-> %s %s (r = %.3f, p = %.2g)\n",
              r$User_Cluster, r$User_Label, r$Study, r$Comparator_Cluster,
              r$Pearson_r, r$P_value))
}
sink()
message("  Saved comparison_summary.txt")

message("\n=== A1.2 Analysis Complete ===")
