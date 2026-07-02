#!/usr/bin/env Rscript
#
# A2: POPS-Prioritized Gene Sets for Each bNMF Cluster
#
# For each cluster variant, finds all protein-coding genes within a window
# (default ±500kb), then selects the gene with the highest max-POPS-score
# across all available traits. This produces a cluster-specific gene set
# that uses functional prioritization rather than simple proximity.
#
# Outputs diagnostic tables and a gene-set file consumed by a2_gtex_heatmaps.R.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"

# -- Config -------------------------------------------------------------------

out_dir <- file.path(base_dir, "results/a2_analysis/pops_trait_selection")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gtf_path  <- file.path(base_dir, "reference/gencode.v19.annotation.gtf.gz")
w_path    <- file.path(base_dir, "results/a1_analysis/META/W_matrix_META.tsv")
pops_path <- file.path(base_dir, "data/PoPS_FullResults.txt")

pops_cohort  <- "UKB"
window_kb    <- 500
window_bp    <- window_kb * 1000
top_pct      <- 0.10   # include genes in top 10% of max-POPS scores

cluster_labels <- c(
  K1 = "Lpa",                       K2 = "Insulin Sensitivity",
  K3 = "Metabolic Syndrome",        K4 = "Lipid Protective",
  K5 = "Atherogenic-Inflammatory",  K6 = "Adiponectin",
  K7 = "SHBG",                      K8 = "Platelet",
  K9 = "ALP"
)

stopifnot(
  "PoPS_FullResults.txt not found" = file.exists(pops_path),
  "W matrix not found"             = file.exists(w_path),
  "GTF not found"                  = file.exists(gtf_path)
)

message("=== A2: POPS-Prioritized Gene Selection ===\n")
message(sprintf("  Cohort: %s | Window: +/-%d kb | Top: %g%%\n",
                pops_cohort, window_kb, top_pct * 100))

# -- Step 1: Parse GTF for protein-coding genes ------------------------------

message("Step 1: Parsing GTF for protein-coding genes...")

gene_df <- read_tsv(
  gtf_path, comment = "#",
  col_names = c("chr", "source", "feature", "start", "end",
                "score", "strand", "frame", "attributes"),
  col_types = cols(chr = col_character(), start = col_double(),
                   end = col_double(), .default = col_character())
) %>%
  filter(feature == "gene",
         str_detect(attributes, 'gene_type "protein_coding"')) %>%
  mutate(
    gene_name = str_replace(attributes, '.*gene_name "([^"]+)".*', "\\1"),
    gene_id   = str_replace(attributes, '.*gene_id "([^"]+)".*', "\\1"),
    gene_id   = str_remove(gene_id, "\\.\\d+$"),
    tss       = ifelse(str_detect(attributes, 'strand "\\+"'), start, end)
  ) %>%
  mutate(chr = str_remove(chr, "^chr")) %>%
  filter(chr %in% as.character(1:22)) %>%
  select(chr, start, end, tss, gene_name, gene_id)

gene_dt <- as.data.table(gene_df)
id_to_symbol <- setNames(gene_df$gene_name, gene_df$gene_id)

message(sprintf("  %d protein-coding genes on autosomes", nrow(gene_df)))

# -- Step 2: Assign variants to clusters ------------------------------------

message("Step 2: Loading W matrix and assigning variants to clusters...")

w_mat <- read_tsv(w_path, show_col_types = FALSE)
k_cols <- grep("^K\\d+$", names(w_mat), value = TRUE)

variant_clusters <- w_mat %>%
  mutate(
    cluster = k_cols[max.col(across(all_of(k_cols)))],
    chr     = str_extract(VAR_ID, "^[^_]+"),
    pos     = as.numeric(str_extract(VAR_ID, "(?<=_)\\d+"))
  ) %>%
  select(VAR_ID, chr, pos, cluster)

var_dt <- as.data.table(variant_clusters)

message(sprintf("  %d variants across %d clusters",
                nrow(var_dt), n_distinct(var_dt$cluster)))

# -- Step 3: Load POPS scores and compute max across traits per gene ---------

message(sprintf("\nStep 3: Loading POPS scores (cohort = %s)...", pops_cohort))

pops_raw <- fread(pops_path)
pops_raw <- pops_raw[cohort == pops_cohort]
pops_raw[, ensgid := str_remove(ensgid, "\\.\\d+$")]
pops_raw <- pops_raw[ensgid %in% gene_df$gene_id]

traits <- unique(pops_raw$trait)
message(sprintf("  %d traits, %d unique genes", length(traits), uniqueN(pops_raw$ensgid)))

# For each gene: max POPS score across all traits and the trait that achieved it
pops_max <- pops_raw[,
  .(max_pops_score = max(pops_score),
    best_trait     = trait[which.max(pops_score)]),
  by = ensgid
]

pops_threshold <- quantile(pops_max$max_pops_score, 1 - top_pct)
n_above <- sum(pops_max$max_pops_score >= pops_threshold)
message(sprintf("  Computed max-across-traits POPS score for %d genes", nrow(pops_max)))
message(sprintf("  Top %g%% threshold: %.3f (%d genes above)",
                top_pct * 100, pops_threshold, n_above))

# -- Step 4: For each variant, find all top-POPS genes within window ---------

message(sprintf("\nStep 4: Finding top-%g%% POPS genes within +/-%d kb of each variant...",
                top_pct * 100, window_kb))

# For each variant, include ALL genes in the window whose max-POPS score
# is in the genome-wide top 10%
pops_prioritized <- var_dt[, {
  g <- gene_dt[chr == .BY$chr]
  if (nrow(g) == 0) {
    list(gene_id = NA_character_, gene_name = NA_character_,
         distance = NA_real_, max_pops_score = NA_real_,
         best_trait = NA_character_, n_window_genes = 0L)
  } else {
    dist_to_gene <- pmin(
      pmax(g$start - pos, 0),
      pmax(pos - g$end, 0)
    )
    in_window <- dist_to_gene <= window_bp

    if (sum(in_window) == 0) {
      list(gene_id = NA_character_, gene_name = NA_character_,
           distance = NA_real_, max_pops_score = NA_real_,
           best_trait = NA_character_, n_window_genes = 0L)
    } else {
      candidates <- g[in_window]
      candidates[, distance := dist_to_gene[in_window]]
      n_window <- nrow(candidates)

      candidates <- merge(candidates, pops_max, by.x = "gene_id", by.y = "ensgid",
                          all.x = TRUE)
      # Keep only genes above the top-10% threshold
      candidates <- candidates[!is.na(max_pops_score) & max_pops_score >= pops_threshold]

      if (nrow(candidates) == 0) {
        list(gene_id = NA_character_, gene_name = NA_character_,
             distance = NA_real_, max_pops_score = NA_real_,
             best_trait = NA_character_, n_window_genes = n_window)
      } else {
        list(gene_id        = candidates$gene_id,
             gene_name      = candidates$gene_name,
             distance       = candidates$distance,
             max_pops_score = candidates$max_pops_score,
             best_trait     = candidates$best_trait,
             n_window_genes = n_window)
      }
    }
  }
}, by = .(VAR_ID, chr, pos, cluster)]

n_var_with_genes <- pops_prioritized[!is.na(gene_id), uniqueN(VAR_ID)]
n_total_hits     <- sum(!is.na(pops_prioritized$gene_id))
message(sprintf("  %d/%d variants had top-%g%% POPS genes in window",
                n_var_with_genes, nrow(var_dt), top_pct * 100))
message(sprintf("  %d total variant-gene pairs (multiple genes per variant allowed)",
                n_total_hits))

# -- Step 5: Build cluster-specific gene sets --------------------------------

message("\nStep 5: Building cluster-specific POPS gene sets...\n")

# Deduplicate: one row per cluster × gene, keep closest variant mapping
gene_set_detail <- pops_prioritized[!is.na(gene_id)]
gene_set_detail <- gene_set_detail[order(cluster, gene_id, distance)]
gene_set_dedup  <- gene_set_detail[, .SD[1], by = .(cluster, gene_id)]
setnames(gene_set_dedup, "gene_id", "ensgid")
setnames(gene_set_dedup, "gene_name", "gene_symbol")

for (cl in names(cluster_labels)) {
  cl_genes <- gene_set_dedup[cluster == cl]
  n_genes  <- nrow(cl_genes)
  top_traits <- cl_genes[order(-max_pops_score), head(unique(best_trait), 5)]

  message(sprintf("  %s (%s): %d POPS genes | top traits: %s",
                  cl, cluster_labels[cl], n_genes,
                  paste(top_traits, collapse = ", ")))
}

# -- Step 6: Write outputs ---------------------------------------------------

message("\nStep 6: Writing output files...")

# Per-variant detail: shows every variant and all its top-POPS genes
out_detail <- file.path(out_dir, "pops_variant_gene_mapping.csv")
fwrite(pops_prioritized, out_detail)
message(sprintf("  %s (%d rows)", basename(out_detail), nrow(pops_prioritized)))

# Deduplicated gene sets for heatmap script consumption
out_genesets <- file.path(out_dir, "pops_cluster_gene_sets.tsv")
fwrite(gene_set_dedup[, .(cluster, ensgid, gene_symbol, best_trait, max_pops_score)],
       out_genesets, sep = "\t")
message(sprintf("  %s (%d rows)", basename(out_genesets), nrow(gene_set_dedup)))

# Summary per cluster
cluster_summary <- gene_set_dedup[,
  .(n_pops_genes       = .N,
    mean_max_pops_score = mean(max_pops_score),
    top_traits          = paste(head(unique(best_trait[order(-max_pops_score)]), 5),
                                collapse = "; ")),
  by = cluster
]
cluster_summary[, cluster_label := cluster_labels[cluster]]
setcolorder(cluster_summary, c("cluster", "cluster_label"))

out_summary <- file.path(out_dir, "pops_cluster_summary.csv")
fwrite(cluster_summary, out_summary)
message(sprintf("  %s (%d rows)", basename(out_summary), nrow(cluster_summary)))

message("\nDone. Review outputs before running a2_gtex_heatmaps.R")
