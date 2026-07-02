#!/usr/bin/env Rscript
# 05_cluster_trait_barplots.R
# Barplots of top 20 H-matrix trait loadings per cluster for each arm.
#
# Usage:
#   Rscript scripts/a1_analysis/05_cluster_trait_barplots.R

library(data.table)
library(ggplot2)

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
setwd(base_dir)

out_dir <- "results/a1_comparison"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Helper: process H matrix into signed loadings ---
process_h_matrix <- function(h_file, strip_source = FALSE) {
  h <- fread(h_file)
  clusters <- h$Cluster
  h[, Cluster := NULL]

  cols <- names(h)
  # Pair pos/neg columns by base trait name
  pos_cols <- grep("_pos$", cols, value = TRUE)
  neg_cols <- grep("_neg$", cols, value = TRUE)

  # Extract base trait names (everything before _pos/_neg)
  pos_bases <- sub("_pos$", "", pos_cols)
  neg_bases <- sub("_neg$", "", neg_cols)

  # Match pairs
  shared <- intersect(pos_bases, neg_bases)

  # Compute signed loadings: pos - neg
  result <- rbindlist(lapply(seq_along(clusters), function(i) {
    rbindlist(lapply(shared, function(trait) {
      pos_val <- as.numeric(h[i, get(paste0(trait, "_pos"))])
      neg_val <- as.numeric(h[i, get(paste0(trait, "_neg"))])
      signed <- pos_val - neg_val

      # Clean trait name: optionally strip GWAS source (e.g., _Verma2024)
      clean_name <- trait
      if (strip_source) {
        clean_name <- sub("_[A-Za-z]+\\d{4}$", "", trait)
      }

      data.table(cluster = clusters[i], trait = clean_name,
                 loading = signed, abs_loading = abs(signed))
    }))
  }))

  # For publication arm: same trait from multiple GWAS → keep the one with max abs loading per cluster
  if (strip_source) {
    result <- result[, .SD[which.max(abs_loading)], by = .(cluster, trait)]
  }

  result
}

# --- Helper: create barplot ---
make_barplot <- function(dt, top_n, arm_prefix, nrow_facet, ncol_facet, title,
                         cluster_names = NULL) {
  # Label clusters
  dt[, cluster_label := paste0(arm_prefix, " ", cluster)]
  k_levels <- paste0(arm_prefix, " K", seq_len(length(unique(dt$cluster))))

  # Rename clusters if mapping provided
  if (!is.null(cluster_names)) {
    dt[, cluster_label := ifelse(cluster_label %in% names(cluster_names),
                                 cluster_names[cluster_label], cluster_label)]
    k_levels <- cluster_names[k_levels]
  }
  dt[, cluster_label := factor(cluster_label, levels = k_levels)]

  # Top N per cluster
  top <- dt[, {
    .SD[order(-abs_loading)][seq_len(min(.N, top_n))]
  }, by = cluster_label]

  # Direction for fill
  top[, direction := fifelse(loading > 0, "Increased", "Decreased")]

  # Create unique trait labels per facet so each facet can have its own ordering
  top[, trait_facet := paste(trait, cluster_label, sep = "__")]
  top[, trait_ordered := reorder(trait_facet, loading), by = cluster_label]

  ggplot(top, aes(x = trait_ordered, y = loading, fill = direction)) +
    geom_col(width = 0.7) +
    scale_x_discrete(labels = function(x) sub("__.*$", "", x)) +
    scale_fill_manual(values = c("Increased" = "#E74C3C", "Decreased" = "#3498DB"),
                      name = "Cluster Weight") +
    coord_flip() +
    facet_wrap(~cluster_label, nrow = nrow_facet, ncol = ncol_facet,
               scales = "free") +
    labs(x = NULL, y = "H-Matrix Loading (pos - neg)", title = title) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 8),
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      panel.spacing = unit(0.8, "lines"),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# --- TAILS arm (K=8) ---
cat("Processing TAILS arm...\n")
tails <- process_h_matrix("results/a1_analysis_splitukb/EUR/H_matrix_EUR.tsv",
                           strip_source = FALSE)

tails_names <- c(
  "TAILS K1" = "Glucose+",
  "TAILS K2" = "Erythrocyte",
  "TAILS K3" = "Metabolic Syndrome",
  "TAILS K4" = "Body Composition-",
  "TAILS K5" = "Triglycerides ApoA",
  "TAILS K6" = "Body Composition+",
  "TAILS K7" = "Reticulocyte",
  "TAILS K8" = "Triglycerides-"
)

p_tails <- make_barplot(tails, top_n = 20, arm_prefix = "TAILS",
                         nrow_facet = 2, ncol_facet = 4,
                         title = "TAILS Arm: Top 20 Trait Loadings per Cluster (K=8)",
                         cluster_names = tails_names)

ggsave(file.path(out_dir, "barplot_tails_trait_loadings.png"),
       p_tails, width = 18, height = 10, dpi = 300)
cat("  Saved: barplot_tails_trait_loadings.png\n")

# --- Publication arm (K=9) ---
cat("Processing Publication arm...\n")
pub <- process_h_matrix("results/a1_analysis/EUR/H_matrix_EUR.tsv",
                         strip_source = TRUE)

dk_names <- c(
  "DK K1" = "Glucose-",
  "DK K2" = "Lpa",
  "DK K3" = "Metabolic Syndrome",
  "DK K4" = "Body Composition-",
  "DK K5" = "CRP ApoB",
  "DK K6" = "Liver",
  "DK K7" = "Erythrocyte",
  "DK K8" = "Glucose+",
  "DK K9" = "Triglycerides-"
)

p_pub <- make_barplot(pub, top_n = 20, arm_prefix = "DK",
                       nrow_facet = 3, ncol_facet = 3,
                       title = "Publication Arm: Top 20 Trait Loadings per Cluster (K=9)",
                       cluster_names = dk_names)

ggsave(file.path(out_dir, "barplot_dk_trait_loadings.png"),
       p_pub, width = 18, height = 14, dpi = 300)
cat("  Saved: barplot_dk_trait_loadings.png\n")

# --- Core Published arm (K=8) ---
cat("Processing Core Published arm...\n")
core_pub <- process_h_matrix("results/a1_analysis_core_published/EUR/H_matrix_EUR.tsv",
                              strip_source = TRUE)

core_pub_names <- c(
  "Core K1" = "Height+",
  "Core K2" = "Height-",
  "Core K3" = "Body Composition-",
  "Core K4" = "Erythrocyte",
  "Core K5" = "Glucose-",
  "Core K6" = "Metabolic Syndrome",
  "Core K7" = "Glucose+",
  "Core K8" = "Triglycerides-"
)

p_core_pub <- make_barplot(core_pub, top_n = 20, arm_prefix = "Core",
                            nrow_facet = 2, ncol_facet = 4,
                            title = "Core Published Arm: Top 20 Trait Loadings per Cluster (K=8)",
                            cluster_names = core_pub_names)

ggsave(file.path(out_dir, "barplot_core_published_trait_loadings.png"),
       p_core_pub, width = 18, height = 10, dpi = 300)
cat("  Saved: barplot_core_published_trait_loadings.png\n")

# --- Core Split-UKB arm (K=6) ---
cat("Processing Core Split-UKB arm...\n")
core_split <- process_h_matrix("results/a1_analysis_core_splitukb/EUR/H_matrix_EUR.tsv",
                                strip_source = FALSE)

core_split_names <- c(
  "CoreSplit K1" = "Glucose+",
  "CoreSplit K2" = "Body Composition-",
  "CoreSplit K3" = "Reticulocytes",
  "CoreSplit K4" = "Metabolic Syndrome",
  "CoreSplit K5" = "Triglycerides ApoA",
  "CoreSplit K6" = "Triglycerides-"
)

p_core_split <- make_barplot(core_split, top_n = 20, arm_prefix = "CoreSplit",
                              nrow_facet = 2, ncol_facet = 3,
                              title = "Core Split-UKB Arm: Top 20 Trait Loadings per Cluster (K=6)",
                              cluster_names = core_split_names)

ggsave(file.path(out_dir, "barplot_core_splitukb_trait_loadings.png"),
       p_core_split, width = 14, height = 10, dpi = 300)
cat("  Saved: barplot_core_splitukb_trait_loadings.png\n")

cat("Done.\n")
