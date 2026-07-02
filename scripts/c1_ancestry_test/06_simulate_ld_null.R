#!/usr/bin/env Rscript
# 06_simulate_ld_null.R
# LD simulation: demonstrate that shared causal genes with ancestry-specific
# LD produce low SNP Jaccard but high gene/pathway Jaccard.
#
# Addresses reviewer concern: "Are SNP differences really just LD, or are
# you overinterpreting?"
#
# Strategy:
#   1. Select N causal genes present in both 1KG EUR and AFR reference panels
#   2. For each causal gene, pick DIFFERENT tag SNPs in EUR vs AFR (simulating
#      different LD structures pointing to the same causal locus)
#   3. Generate synthetic GWAS sumstats with signal at tag SNPs + LD-correlated SNPs
#   4. Run the clumping + gene mapping + pathway pipeline
#   5. Show: SNP Jaccard LOW, Gene Jaccard HIGH, Pathway metrics HIGH
#
# Usage:
#   Rscript scripts/c1_ancestry_test/06_simulate_ld_null.R
#   Rscript scripts/c1_ancestry_test/06_simulate_ld_null.R --config config/c1_config.yaml

library(data.table)
library(yaml)
library(ggplot2)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/c1_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}

# --- Load config ---
cfg <- yaml::read_yaml(config_path)

# --- Determine project root ---
project_root <- getwd()
if (!file.exists(file.path(project_root, "scripts", "c1_ancestry_test", "jaccard_utils.R"))) {
  script_dir <- dirname(sys.frame(1)$ofile %||% ".")
  project_root <- normalizePath(file.path(script_dir, "..", ".."))
}

# --- Source utilities ---
source(file.path(project_root, "scripts", "c1_ancestry_test", "jaccard_utils.R"))
source(file.path(project_root, "scripts", "a1_analysis", "figure_utils.R"))

# --- Configuration ---
sim_dir <- file.path(project_root, cfg$results_dir, "simulation")
dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)

eur_ref_prefix <- file.path(project_root, cfg$ld_clump$ref_panel$EUR)
afr_ref_prefix <- file.path(project_root, cfg$ld_clump$ref_panel$AFR)
gtf_path       <- cfg$gtf_path
gene_flank_kb  <- if (is.null(cfg$jaccard$gene_flank_kb)) 50 else cfg$jaccard$gene_flank_kb
gene_flank_bp  <- gene_flank_kb * 1000
top_n          <- cfg$jaccard$top_n
pathway_gmt    <- cfg$jaccard$pathway_gmt
pathway_fdr    <- if (is.null(cfg$jaccard$pathway_fdr)) 0.05 else cfg$jaccard$pathway_fdr

n_causal       <- 75      # causal genes to simulate
n_clusters     <- 3       # split causal genes into this many simulated clusters
sim_seed       <- 42
sample_size    <- 200000  # simulated GWAS sample size

cat(sprintf("\n%s\n=== C1.1 LD Simulation ===\n%s\n",
            strrep("=", 50), strrep("=", 50)))
cat(sprintf("  Causal genes: %d (split into %d clusters)\n", n_causal, n_clusters))
cat(sprintf("  EUR ref: %s\n", eur_ref_prefix))
cat(sprintf("  AFR ref: %s\n", afr_ref_prefix))


# =====================================================================
# STEP 1: Load reference panel BIM files
# =====================================================================

cat("\n--- Step 1: Loading reference panel BIM files ---\n")

load_all_bim <- function(prefix) {
  rbindlist(lapply(1:22, function(chr) {
    bim_file <- paste0(prefix, ".", chr, ".bim")
    if (!file.exists(bim_file)) return(data.table())
    bim <- fread(bim_file, header = FALSE,
                 col.names = c("CHR", "RSID", "CM", "POS", "A1", "A2"))
    # Build VAR_ID with sorted alleles
    bim[, allele_min := pmin(A1, A2)]
    bim[, allele_max := pmax(A1, A2)]
    bim[, VAR_ID := paste(CHR, POS, allele_min, allele_max, sep = "_")]
    bim[, .(CHR, POS, RSID, A1, A2, VAR_ID)]
  }))
}

eur_bim <- load_all_bim(eur_ref_prefix)
afr_bim <- load_all_bim(afr_ref_prefix)
cat(sprintf("  EUR reference: %d SNPs\n", nrow(eur_bim)))
cat(sprintf("  AFR reference: %d SNPs\n", nrow(afr_bim)))


# =====================================================================
# STEP 2: Select causal genes with different tag SNPs per ancestry
# =====================================================================

cat("\n--- Step 2: Selecting causal genes ---\n")

gene_dt <- parse_gtf_genes_dt(gtf_path)
cat(sprintf("  Protein-coding genes: %d\n", nrow(gene_dt)))

set.seed(sim_seed)

# For each gene, find SNPs within gene body +/- flank in both panels
eligible_genes <- rbindlist(lapply(seq_len(nrow(gene_dt)), function(gi) {
  g <- gene_dt[gi]
  region_start <- max(0L, as.integer(g$start - gene_flank_bp))
  region_end   <- as.integer(g$end + gene_flank_bp)

  eur_in_region <- eur_bim[CHR == g$chr & POS >= region_start & POS <= region_end]
  afr_in_region <- afr_bim[CHR == g$chr & POS >= region_start & POS <= region_end]

  if (nrow(eur_in_region) >= 5 && nrow(afr_in_region) >= 5) {
    data.table(
      gene_name = g$gene_name, gene_id = g$gene_id,
      chr = g$chr, start = g$start, end = g$end,
      n_eur_snps = nrow(eur_in_region),
      n_afr_snps = nrow(afr_in_region)
    )
  }
}))

cat(sprintf("  Genes with >=5 SNPs in both panels: %d\n", nrow(eligible_genes)))

if (nrow(eligible_genes) < n_causal) {
  cat(sprintf("  WARNING: Only %d eligible genes, reducing n_causal.\n", nrow(eligible_genes)))
  n_causal <- nrow(eligible_genes)
  n_clusters <- min(n_clusters, n_causal)
}

# Randomly select causal genes
causal_genes <- eligible_genes[sample(.N, n_causal)]

# For each causal gene, pick one EUR tag SNP and a DIFFERENT AFR tag SNP
causal_dt <- rbindlist(lapply(seq_len(nrow(causal_genes)), function(i) {
  g <- causal_genes[i]
  region_start <- max(0L, as.integer(g$start - gene_flank_bp))
  region_end   <- as.integer(g$end + gene_flank_bp)

  eur_snps <- eur_bim[CHR == g$chr & POS >= region_start & POS <= region_end]
  afr_snps <- afr_bim[CHR == g$chr & POS >= region_start & POS <= region_end]

  # Pick a random EUR SNP
  eur_pick <- eur_snps[sample(.N, 1)]

  # Pick an AFR SNP at a DIFFERENT position (to simulate LD-driven tag difference)
  afr_candidates <- afr_snps[POS != eur_pick$POS]
  if (nrow(afr_candidates) == 0) {
    afr_candidates <- afr_snps  # fallback: allow same position
  }
  afr_pick <- afr_candidates[sample(.N, 1)]

  data.table(
    gene_name  = g$gene_name, gene_id = g$gene_id,
    chr        = g$chr, start = g$start, end = g$end,
    eur_var_id = eur_pick$VAR_ID, eur_pos = eur_pick$POS,
    afr_var_id = afr_pick$VAR_ID, afr_pos = afr_pick$POS,
    same_pos   = (eur_pick$POS == afr_pick$POS)
  )
}))

# Assign to clusters
causal_dt[, cluster := paste0("K", rep(1:n_clusters, length.out = .N))]

cat(sprintf("  Selected %d causal genes across %d clusters\n", nrow(causal_dt), n_clusters))
cat(sprintf("  EUR-AFR tag SNPs at same position: %d/%d\n",
            sum(causal_dt$same_pos), nrow(causal_dt)))

fwrite(causal_dt, file.path(sim_dir, "causal_genes.tsv"), sep = "\t")


# =====================================================================
# STEP 3: Generate synthetic GWAS summary statistics
# =====================================================================

cat("\n--- Step 3: Generating synthetic GWAS ---\n")

generate_synthetic_gwas <- function(bim_dt, tag_var_ids, n_traits = 7,
                                     sample_size = 200000, seed = 42) {
  set.seed(seed)

  # All SNPs from reference panel
  all_var_ids <- unique(bim_dt$VAR_ID)

  # For each trait: causal SNPs get significant signal, rest get null
  trait_names <- paste0("trait_", 1:n_traits)
  gwas_list <- list()

  for (t_i in seq_along(trait_names)) {
    trait <- trait_names[t_i]

    # Null GWAS: uniform p-values
    dt <- data.table(
      VAR_ID        = all_var_ids,
      Effect_Allele = "A",
      P_VALUE       = runif(length(all_var_ids)),
      BETA          = rnorm(length(all_var_ids), 0, 0.001),
      SE            = rep(0.01, length(all_var_ids)),
      N             = sample_size,
      MAF           = runif(length(all_var_ids), 0.05, 0.5)
    )

    # Causal SNPs: significant p-values and realistic betas
    is_causal <- dt$VAR_ID %in% tag_var_ids
    n_causal_found <- sum(is_causal)

    if (n_causal_found > 0) {
      dt[is_causal, P_VALUE := 10^(-runif(.N, 8, 20))]
      dt[is_causal, BETA := rnorm(.N, 0, 0.05)]
      dt[is_causal, SE := abs(BETA) / sqrt(-log10(P_VALUE) * 2)]

      # LD-based signal spreading: nearby SNPs get attenuated signal
      causal_positions <- bim_dt[VAR_ID %in% tag_var_ids, .(CHR, POS, VAR_ID)]
      for (ci in seq_len(nrow(causal_positions))) {
        cp <- causal_positions[ci]
        causal_p <- dt[VAR_ID == cp$VAR_ID, P_VALUE]

        # Find nearby non-causal SNPs (within 250kb)
        nearby <- bim_dt[CHR == cp$CHR &
                          abs(POS - cp$POS) > 0 &
                          abs(POS - cp$POS) <= 250000 &
                          !(VAR_ID %in% tag_var_ids)]

        if (nrow(nearby) > 0) {
          distances <- abs(nearby$POS - cp$POS)
          # Attenuated signal proportional to distance
          spread_log10p <- -log10(causal_p) * exp(-distances / 100000)
          spread_p <- 10^(-spread_log10p)

          # Only update if spread signal is more significant than null
          for (ni in seq_len(nrow(nearby))) {
            idx <- which(dt$VAR_ID == nearby$VAR_ID[ni])
            if (length(idx) == 1 && spread_p[ni] < dt$P_VALUE[idx]) {
              dt$P_VALUE[idx] <- spread_p[ni]
              dt$BETA[idx] <- dt$BETA[idx] + rnorm(1, 0, 0.02) * exp(-distances[ni] / 100000)
            }
          }
        }
      }
    }

    cat(sprintf("    %s: %d causal SNPs found in panel, %d significant (P<5e-8)\n",
                trait, n_causal_found, sum(dt$P_VALUE < 5e-8)))

    gwas_list[[trait]] <- dt
  }

  gwas_list
}

eur_gwas <- generate_synthetic_gwas(eur_bim, causal_dt$eur_var_id,
                                     n_traits = 7, sample_size = sample_size,
                                     seed = sim_seed)
afr_gwas <- generate_synthetic_gwas(afr_bim, causal_dt$afr_var_id,
                                     n_traits = 7, sample_size = sample_size,
                                     seed = sim_seed + 1)

# Write synthetic GWAS to disk
eur_gwas_dir <- file.path(sim_dir, "eur_gwas")
afr_gwas_dir <- file.path(sim_dir, "afr_gwas")
dir.create(eur_gwas_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(afr_gwas_dir, recursive = TRUE, showWarnings = FALSE)

eur_gwas_files <- list()
afr_gwas_files <- list()
for (trait in names(eur_gwas)) {
  eur_file <- file.path(eur_gwas_dir, sprintf("sim_%s_EUR.txt.gz", trait))
  afr_file <- file.path(afr_gwas_dir, sprintf("sim_%s_AFR.txt.gz", trait))
  fwrite(eur_gwas[[trait]], eur_file, sep = "\t")
  fwrite(afr_gwas[[trait]], afr_file, sep = "\t")
  eur_gwas_files[[trait]] <- eur_file
  afr_gwas_files[[trait]] <- afr_file
}


# =====================================================================
# STEP 4: Build synthetic W matrices from GWAS signal
# =====================================================================

cat("\n--- Step 4: Building synthetic W matrices ---\n")

# For each ancestry, build a W matrix where cluster membership is
# defined by the causal gene assignments, and weights are -log10(p)
build_sim_W <- function(gwas_list, tag_var_ids_by_cluster, all_var_ids) {
  # Use the first trait's p-values as weights
  trait1 <- gwas_list[[1]]

  cluster_names <- names(tag_var_ids_by_cluster)
  W <- data.table(VAR_ID = all_var_ids)

  for (k in cluster_names) {
    cluster_tags <- tag_var_ids_by_cluster[[k]]
    weights <- rep(0, length(all_var_ids))

    # Tag SNPs get their -log10(p) weight
    for (tag in cluster_tags) {
      idx <- which(all_var_ids == tag)
      if (length(idx) == 1) {
        p_val <- trait1[VAR_ID == tag, P_VALUE]
        if (length(p_val) == 1 && !is.na(p_val) && p_val > 0) {
          weights[idx] <- -log10(p_val)
        }
      }
    }

    W[, (k) := weights]
  }

  W
}

# Group tag SNPs by cluster
eur_tags_by_cluster <- split(causal_dt$eur_var_id, causal_dt$cluster)
afr_tags_by_cluster <- split(causal_dt$afr_var_id, causal_dt$cluster)

# Use only SNPs that passed a basic p-value threshold in at least one trait
eur_sig_snps <- unique(unlist(lapply(eur_gwas, function(dt) dt[P_VALUE < 1e-5, VAR_ID])))
afr_sig_snps <- unique(unlist(lapply(afr_gwas, function(dt) dt[P_VALUE < 1e-5, VAR_ID])))

cat(sprintf("  EUR significant SNPs (P<1e-5): %d\n", length(eur_sig_snps)))
cat(sprintf("  AFR significant SNPs (P<1e-5): %d\n", length(afr_sig_snps)))

eur_W <- build_sim_W(eur_gwas, eur_tags_by_cluster, eur_sig_snps)
afr_W <- build_sim_W(afr_gwas, afr_tags_by_cluster, afr_sig_snps)

cat(sprintf("  EUR W matrix: %d variants x %d clusters\n",
            nrow(eur_W), ncol(eur_W) - 1))
cat(sprintf("  AFR W matrix: %d variants x %d clusters\n",
            nrow(afr_W), ncol(afr_W) - 1))


# =====================================================================
# STEP 5: Compute Jaccard at all levels
# =====================================================================

cat("\n--- Step 5: Computing Jaccard similarities ---\n")

# Build SNP-to-gene mapping for all simulation variants
all_sim_var_ids <- unique(c(eur_W$VAR_ID, afr_W$VAR_ID))
snp_to_gene_map <- build_snp_to_gene_map(all_sim_var_ids, gene_dt, flank_bp = gene_flank_bp)
cat(sprintf("  SNPs mapped to genes: %d / %d\n",
            length(unique(snp_to_gene_map$VAR_ID)), length(all_sim_var_ids)))

# Load pathway data if available
pathway_map <- NULL
background_gene_ids <- NULL
has_pathway <- FALSE

if (!is.null(pathway_gmt) && file.exists(pathway_gmt)) {
  pathway_map <- parse_gmt(pathway_gmt)
  background_gene_ids <- unique(snp_to_gene_map$gene_id)
  has_pathway <- TRUE
  cat(sprintf("  Pathways loaded: %d\n", length(unique(pathway_map$pathway))))
}

# Match EUR vs AFR clusters
cluster_names <- setdiff(colnames(eur_W), "VAR_ID")
sim_results <- rbindlist(lapply(cluster_names, function(k) {
  eur_snps <- get_cluster_snp_set(eur_W, k, top_n = top_n)
  afr_snps <- get_cluster_snp_set(afr_W, k, top_n = top_n)

  # SNP-level Jaccard
  snp_j <- jaccard_similarity(eur_snps, afr_snps)

  # Gene-level Jaccard
  gene_j <- gene_level_jaccard(eur_snps, afr_snps, snp_to_gene_map)

  # Shared genes count
  eur_genes <- unique(snp_to_gene_map[VAR_ID %in% eur_snps, gene_name])
  afr_genes <- unique(snp_to_gene_map[VAR_ID %in% afr_snps, gene_name])
  shared_genes <- intersect(eur_genes, afr_genes)

  # Expected shared genes (genes in causal set for this cluster)
  expected_genes <- causal_dt[cluster == k, gene_name]

  result <- data.table(
    cluster              = k,
    snp_jaccard          = snp_j,
    gene_jaccard         = gene_j,
    n_eur_genes          = length(eur_genes),
    n_afr_genes          = length(afr_genes),
    n_shared_genes       = length(shared_genes),
    n_expected_causal    = length(expected_genes),
    causal_gene_recovery = length(intersect(shared_genes, expected_genes)) / length(expected_genes)
  )

  # Pathway-level metrics
  if (has_pathway) {
    pw_j <- pathway_level_jaccard(eur_snps, afr_snps, snp_to_gene_map,
                                    pathway_map, background_gene_ids, pathway_fdr)
    ranked <- ranked_pathway_similarity(eur_snps, afr_snps, snp_to_gene_map,
                                          pathway_map, background_gene_ids)
    result[, pathway_jaccard := pw_j]
    result[, pathway_spearman := ranked$pathway_spearman]
    result[, pathway_top25_overlap := ranked$pathway_top25_overlap]
  }

  cat(sprintf("  %s: SNP J=%.3f, Gene J=%.3f, Shared genes=%d/%d%s\n",
              k, snp_j, gene_j, length(shared_genes), length(expected_genes),
              if (has_pathway) sprintf(", Pathway J=%.3f, Spearman=%.3f",
                pw_j, ranked$pathway_spearman %||% NA) else ""))

  result
}))


# =====================================================================
# STEP 6: Save results and generate figure
# =====================================================================

cat("\n--- Step 6: Saving results and generating figure ---\n")

fwrite(sim_results, file.path(sim_dir, "simulation_summary.tsv"), sep = "\t")
cat(sprintf("  Written: simulation_summary.tsv\n"))

# Build the figure: grouped bar chart comparing levels
metric_cols <- c("snp_jaccard", "gene_jaccard")
metric_labels <- c("SNP", "Gene")
if (has_pathway) {
  metric_cols <- c(metric_cols, "pathway_jaccard", "pathway_spearman")
  metric_labels <- c(metric_labels, "Pathway (Jaccard)", "Pathway (Spearman)")
}

sim_long <- melt(sim_results, id.vars = "cluster",
                  measure.vars = metric_cols,
                  variable.name = "level", value.name = "similarity")
sim_long[, level := factor(level, levels = metric_cols, labels = metric_labels)]

bar_colors <- c("SNP" = "steelblue4", "Gene" = "seagreen4",
                "Pathway (Jaccard)" = "mediumpurple4", "Pathway (Spearman)" = "darkorange3")

p_sim <- ggplot(sim_long, aes(x = cluster, y = similarity, fill = level)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = bar_colors, name = "Comparison Level") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  labs(
    title = "LD Simulation: Shared Causal Genes, Different Tag SNPs",
    subtitle = sprintf("Expected: SNP Jaccard LOW (different tags), Gene/Pathway HIGH (same loci)\n%d causal genes per cluster, top-%d SNPs compared",
                        n_causal %/% n_clusters, top_n),
    x = "Simulated Cluster",
    y = "Similarity"
  ) +
  theme_big_text(base_size = 14) +
  coord_cartesian(ylim = c(min(0, min(sim_long$similarity, na.rm = TRUE) - 0.1),
                            max(1, max(sim_long$similarity, na.rm = TRUE) + 0.1)))

sim_fig_file <- file.path(sim_dir, "simulation_figure")
ggsave(paste0(sim_fig_file, ".png"), p_sim, width = 10, height = 7, dpi = 300)
cat(sprintf("  Written: simulation_figure.png\n"))

# Also make a summary text panel
cat("\n--- Simulation Summary ---\n")
cat(sprintf("  Mean SNP Jaccard:  %.4f\n", mean(sim_results$snp_jaccard)))
cat(sprintf("  Mean Gene Jaccard: %.4f\n", mean(sim_results$gene_jaccard)))
if (has_pathway) {
  cat(sprintf("  Mean Pathway Jaccard:  %.4f\n", mean(sim_results$pathway_jaccard, na.rm = TRUE)))
  cat(sprintf("  Mean Pathway Spearman: %.4f\n", mean(sim_results$pathway_spearman, na.rm = TRUE)))
}
cat(sprintf("  Gene/SNP ratio:    %.2f\n",
            mean(sim_results$gene_jaccard) / max(mean(sim_results$snp_jaccard), 1e-10)))
cat(sprintf("  Causal gene recovery: %.1f%%\n",
            100 * mean(sim_results$causal_gene_recovery)))

cat("\n=== LD Simulation complete ===\n")
cat(sprintf("Results in: %s\n", sim_dir))
cat("Files:\n")
cat("  causal_genes.tsv        - Selected causal genes + EUR/AFR tag SNPs\n")
cat("  simulation_summary.tsv  - Per-cluster Jaccard at all levels\n")
cat("  simulation_figure.png   - Grouped bar chart (SNP vs Gene vs Pathway)\n")
