#!/usr/bin/env Rscript
#
# A2: CATLAS Weighted Specificity Score (WSS) Summary Heatmap
#
# Computes WSS enrichment across pancreatic and esophageal cell types from
# CATLAS cCRE data, combining both tissues into a single summary heatmap.
#
# WSS(k, c) = sum( W[i,k] * specificity[gene_i, c] )
# where specificity = cCRE_count[gene, cell_type] / sum(cCRE_count[gene, all])
#
# Usage:
#   Rscript scripts/a2_analysis/a2_catlas_wss_summary.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"

# -- Config -------------------------------------------------------------------

out_dir <- file.path(base_dir, "results/a2_analysis/catlas/catlas_summary")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gtf_path <- file.path(base_dir, "reference/gencode.v19.annotation.gtf.gz")
w_path   <- file.path(base_dir, "results/a1_analysis/META/W_matrix_META.tsv")

n_perm <- 10000
window <- 50000
set.seed(42)  # reproducible permutation p-values

cluster_labels <- c(
  K1  = "Lpa",
  K2  = "Adiponectin",
  K3  = "Platelet",
  K4  = "SHBG",
  K5  = "Blood Pressure-Stature",
  K6  = "Metabolic",
  K7  = "Triglycerides-HDL",
  K8  = "ALP-LDL",
  K9  = "Obesity",
  K10 = "Glycemic"
)

# Tissue configs: directory, nice labels
tissue_configs <- list(
  Pancreas = list(
    ccre_dir = file.path(base_dir, "data/catlas/pancreas/hg19"),
    labels = c(
      Pancreatic_Acinar_Cell         = "Acinar Cell",
      Pancreatic_Alpha_Cell_1        = "Alpha Cell 1",
      Pancreatic_Alpha_Cell_2        = "Alpha Cell 2",
      Pancreatic_Beta_Cell_1         = "Beta Cell 1",
      Pancreatic_Beta_Cell_2         = "Beta Cell 2",
      `Pancreatic_Delta,Gamma_cell`  = "Delta/Gamma Cell",
      Ductal_Cell_Pancreatic         = "Ductal Cell",
      Fetal_Pancreatic_Acinar_Cell_1 = "Fetal Acinar Cell 1",
      Fetal_Pancreatic_Acinar_Cell_2 = "Fetal Acinar Cell 2",
      Fetal_Pancreatic_Ductal_Cell   = "Fetal Ductal Cell",
      Fetal_Pancreatic_Islet_Cell    = "Fetal Islet Cell"
    )
  ),
  Esophagus = list(
    ccre_dir = file.path(base_dir, "data/catlas/esophagus/hg19"),
    labels = c(
      Esophageal_Epithelial_Cell            = "Epithelial Cell",
      Smooth_Muscle_Esophageal_Mucosal      = "Smooth Muscle (Mucosal)",
      Smooth_Muscle_Esophageal_Muscularis_1 = "Smooth Muscle (Muscularis 1)",
      Smooth_Muscle_Esophageal_Muscularis_2 = "Smooth Muscle (Muscularis 2)",
      Smooth_Muscle_Esophageal_Muscularis_3 = "Smooth Muscle (Muscularis 3)"
    )
  )
)

message("=== A2: CATLAS WSS Summary Heatmap (Pancreas + Esophagus) ===\n")

# -- Step 1: Parse GTF --------------------------------------------------------

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
    chr = str_remove(chr, "^chr")
  ) %>%
  filter(chr %in% as.character(1:22)) %>%
  select(chr, start, end, gene_name)

gene_dt <- as.data.table(gene_df)
message(sprintf("  %d protein-coding genes on autosomes", nrow(gene_dt)))

# -- Step 2: Load W matrix and map variants to genes --------------------------

message("Step 2: Loading W matrix and mapping variants to genes...")

w_mat <- read_tsv(w_path, show_col_types = FALSE)
k_cols <- grep("^K\\d+$", names(w_mat), value = TRUE)

var_dt <- as.data.table(w_mat)[, .(
  VAR_ID = VAR_ID,
  chr = str_extract(VAR_ID, "^[^_]+"),
  pos = as.numeric(str_extract(VAR_ID, "(?<=_)\\d+"))
)]

nearest_genes <- var_dt[, {
  g <- gene_dt[chr == .BY$chr]
  if (nrow(g) == 0) {
    list(gene_name = NA_character_)
  } else {
    dist_to_gene <- pmax(g$start - pos, 0) + pmax(pos - g$end, 0)
    list(gene_name = g$gene_name[which.min(dist_to_gene)])
  }
}, by = .(VAR_ID, chr, pos)]

var_gene <- nearest_genes[!is.na(gene_name), .(VAR_ID, gene_name)]
var_gene <- var_gene[!duplicated(VAR_ID)]

w_long <- as.data.table(w_mat)[, c("VAR_ID", k_cols), with = FALSE]
w_gene <- merge(w_long, var_gene, by = "VAR_ID")

message(sprintf("  %d variants with gene mappings and W weights", nrow(w_gene)))

# -- Step 3: Build gene windows for cCRE overlap ------------------------------

gene_windows <- copy(gene_dt)
gene_windows[, win_start := pmax(start - window, 0L)]
gene_windows[, win_end   := end + window]
setkey(gene_windows, chr, win_start, win_end)

# -- Step 4: WSS per tissue ---------------------------------------------------

message("Step 3: Computing WSS for each tissue...\n")

wss_results <- list()

for (tissue_name in names(tissue_configs)) {
  tc <- tissue_configs[[tissue_name]]
  message(sprintf("  --- WSS for tissue: %s ---", tissue_name))

  bed_files <- list.files(tc$ccre_dir, pattern = "\\.bed$", full.names = TRUE)
  bed_files <- bed_files[!grepl("_unmapped\\.bed$", bed_files)]
  cell_type_ids <- tools::file_path_sans_ext(basename(bed_files))

  # Build gene × cell-type cCRE count matrix
  ccre_mat <- matrix(0L, nrow = nrow(gene_dt), ncol = length(cell_type_ids),
                     dimnames = list(gene_dt$gene_name, cell_type_ids))

  for (i in seq_along(bed_files)) {
    bed <- fread(bed_files[i], header = FALSE, select = 1:3,
                 col.names = c("chr", "start", "end"))
    bed[, chr := sub("^chr", "", chr)]
    bed <- bed[chr %in% as.character(1:22)]
    setkey(bed, chr, start, end)

    overlaps <- foverlaps(bed, gene_windows, by.x = c("chr", "start", "end"),
                          by.y = c("chr", "win_start", "win_end"),
                          nomatch = NULL)
    if (nrow(overlaps) > 0) {
      gc <- overlaps[, .N, by = gene_name]
      idx <- match(gc$gene_name, gene_dt$gene_name)
      ccre_mat[idx[!is.na(idx)], i] <- gc$N[!is.na(idx)]
    }
  }

  # Specificity matrix
  row_totals <- rowSums(ccre_mat)
  spec_full <- ccre_mat / row_totals
  spec_full[!is.finite(spec_full)] <- 0

  # Subset to variants whose nearest gene has cCRE data
  var_in_tissue <- w_gene[gene_name %in% gene_dt$gene_name]
  if (nrow(var_in_tissue) < 2) {
    message(sprintf("    Skipping %s: too few variants", tissue_name))
    next
  }

  gene_idx <- match(var_in_tissue$gene_name, gene_dt$gene_name)
  spec_mat <- spec_full[gene_idx, , drop = FALSE]
  w_weights <- as.matrix(var_in_tissue[, ..k_cols])

  # Observed WSS
  obs_wss <- crossprod(w_weights, spec_mat)

  # Permutation null
  n_variants <- nrow(spec_mat)
  n_all <- nrow(spec_full)
  exceed_count <- matrix(0L, nrow = length(k_cols), ncol = length(cell_type_ids))
  for (p in seq_len(n_perm)) {
    perm_idx <- sample.int(n_all, n_variants)
    perm_spec <- spec_full[perm_idx, , drop = FALSE]
    perm_wss <- crossprod(w_weights, perm_spec)
    exceed_count <- exceed_count + (perm_wss >= obs_wss)
  }
  p_mat <- exceed_count / n_perm

  for (ki in seq_along(k_cols)) {
    for (ci in seq_along(cell_type_ids)) {
      wss_results[[length(wss_results) + 1]] <- data.table(
        tissue = tissue_name,
        cell_type = cell_type_ids[ci],
        cluster = k_cols[ki],
        wss_observed = obs_wss[ki, ci],
        p_value = p_mat[ki, ci]
      )
    }
  }

  message(sprintf("    %d variants x %d cell types x %d clusters",
                  n_variants, length(cell_type_ids), length(k_cols)))
}

# -- Step 5: FDR and save CSV -------------------------------------------------

message("\nStep 4: Applying FDR correction and saving results...")

wss_dt <- rbindlist(wss_results)
wss_dt[, log10_p := -log10(pmax(p_value, 1 / (n_perm + 1)))]
wss_dt[, fdr := p.adjust(p_value, method = "BH"), by = cluster]

for (cl in k_cols) {
  n_cl  <- sum(wss_dt$cluster == cl)
  n_sig <- sum(wss_dt$cluster == cl & wss_dt$fdr < 0.01)
  message(sprintf("  %s (%s): %d tests, %d significant at FDR < 0.01",
                  cl, cluster_labels[cl], n_cl, n_sig))
}

wss_csv <- file.path(out_dir, "catlas_wss_results.csv")
fwrite(wss_dt, wss_csv)
message(sprintf("  Saved %s (%d rows)", basename(wss_csv), nrow(wss_dt)))

# -- Step 6: Summary heatmap --------------------------------------------------

message("\nStep 5: Generating summary heatmap...")

sig_wss <- wss_dt[fdr < 0.01]
message(sprintf("  Total FDR < 0.01: %d", nrow(sig_wss)))

if (nrow(sig_wss) == 0) {
  message("  No significant results — skipping summary heatmap.")
  # Fall back to FDR < 0.05
  sig_wss <- wss_dt[fdr < 0.05]
  message(sprintf("  Trying FDR < 0.05: %d significant", nrow(sig_wss)))
  if (nrow(sig_wss) == 0) {
    message("  Still no significant results. Showing top results by p-value instead.")
    sig_wss <- wss_dt[p_value < 0.05]
    message(sprintf("  Nominal p < 0.05: %d results", nrow(sig_wss)))
  }
}

if (nrow(sig_wss) > 0) {
  # Build nice display labels
  all_labels <- do.call(c, lapply(tissue_configs, `[[`, "labels"))

  sig_wss[, display_label := {
    nice <- all_labels[cell_type]
    ifelse(is.na(nice),
           paste0(tissue, " — ", str_replace_all(cell_type, "_", " ")),
           paste0(tissue, " — ", nice))
  }]

  # Columns in the desired META display order (Glycemic -> Lpa); any clusters
  # not in the list fall back to natural K order at the end.
  desired_k_order <- c("K10", "K9", "K4", "K2", "K7", "K8", "K6", "K3", "K5", "K1")
  present_clusters <- unique(sig_wss$cluster)
  sig_clusters <- c(desired_k_order[desired_k_order %in% present_clusters],
                    sort(setdiff(present_clusters, desired_k_order)))
  sig_labels   <- sort(unique(sig_wss$display_label))

  mat <- matrix(NA_real_, nrow = length(sig_labels), ncol = length(sig_clusters),
                dimnames = list(sig_labels, cluster_labels[sig_clusters]))

  for (i in seq_len(nrow(sig_wss))) {
    rl <- sig_wss$display_label[i]
    cl <- cluster_labels[sig_wss$cluster[i]]
    mat[rl, cl] <- sig_wss$log10_p[i]
  }

  # Tissue annotation
  row_tissues <- str_extract(sig_labels, "^[^—]+") |> str_trim()
  tissue_colors <- c("Pancreas" = "#66C2A5", "Esophagus" = "#FC8D62")

  row_ha <- rowAnnotation(
    Tissue = row_tissues,
    col = list(Tissue = tissue_colors),
    show_legend = TRUE,
    annotation_legend_param = list(Tissue = list(title = "Tissue"))
  )

  # Purple color scale
  max_val <- ceiling(max(mat, na.rm = TRUE))
  max_val <- max(max_val, 1)
  col_fun <- colorRamp2(
    seq(0, max_val, length.out = 10),
    colorRampPalette(c("#F2F0F7", "#CBC9E2", "#9E9AC8", "#756BB1", "#54278F"))(10)
  )

  mat[is.na(mat)] <- 0

  n_rows  <- length(sig_labels)
  panel_h <- max(6, n_rows * 0.25 + 2)

  hm <- Heatmap(
    mat,
    col = col_fun,
    na_col = "white",
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_order = sig_labels,
    column_order = cluster_labels[sig_clusters],
    row_names_side = "left",
    row_names_max_width = unit(12, "cm"),
    row_names_gp = gpar(fontsize = 10),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 10),
    rect_gp = gpar(col = "white", lwd = 1),
    heatmap_legend_param = list(title = "-log10(p)", direction = "horizontal"),
    left_annotation = row_ha,
    name = "summary"
  )

  fdr_threshold <- if (any(wss_dt$fdr < 0.01)) "0.01" else if (any(wss_dt$fdr < 0.05)) "0.05" else "nominal p < 0.05"

  summary_path <- file.path(out_dir, "catlas_summary_heatmap.png")
  png(summary_path, width = 12, height = panel_h, units = "in", res = 600)
  draw(hm,
       column_title = sprintf("CATLAS: Significant Cell Types (WSS, per-cluster FDR < %s)", fdr_threshold),
       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
       heatmap_legend_side = "bottom")
  dev.off()
  message(sprintf("  Saved %s", basename(summary_path)))
}

message("\nDone.")
