#!/usr/bin/env Rscript
# Quick regeneration of the WSS summary heatmap from pre-computed results.
# Reads tabula_sapiens_wss_results.csv and rebuilds the FDR-filtered heatmap.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"

out_dir <- file.path(base_dir, "results/a2_analysis/tabula_sapiens")
ts_data_dir <- file.path(base_dir,
  "data/tabula_sapiens/tabula_sapiens_cell_types/tabula_sapiens_processed_data")
all_meta_path <- file.path(ts_data_dir, "all_metadata.csv")

# Canonical META labels (source of truth: scripts/a1_analysis/06_ancestry_trait_barplots.R);
# this WSS uses the same W_matrix_META.tsv as the main tabula/catlas scripts, so the
# K1..K10 mapping is identical. (Previous values here were a stale, non-canonical set.)
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

merge_rules <- list(
  Blood = list(
    list(cols = c("NK_cell", "nk_cell"), keep = "NK_cell", label = "NK Cells")
  ),
  Large_Intestine = list(
    list(cols = c("cd8pos__T_cell", "cd8pos_t_cell"),
         keep = "cd8pos__T_cell", label = "CD8-positive alpha-beta T Cells"),
    list(cols = c("intestinal_transient_amplifying_cell", "transit_amplifying_cell_of_lg_int"),
         keep = "intestinal_transient_amplifying_cell", label = "Intestinal Transient Amplifying Cells")
  ),
  Small_Intestine = list(
    list(cols = c("Cd8pos__T_cell", "cd8pos_t_cell"),
         keep = "Cd8pos__T_cell", label = "CD8 alpha-beta T Cells"),
    list(cols = c("intestinal_transient_amplifying_cell", "transit_amplifying_cell_of_small_int"),
         keep = "intestinal_transient_amplifying_cell", label = "Intestinal Transient Amplifying Cells"),
    list(cols = c("paneth_cell", "paneth_cell_of_epithelium_of_small_int"),
         keep = "paneth_cell", label = "Paneth Cells")
  ),
  Lymph_Node = list(
    list(cols = c("CD141positive_myeloid_dendritic_cell", "cd141positive_myeloid_dendritic_cell"),
         keep = "CD141positive_myeloid_dendritic_cell", label = "CD141 Myeloid Dendritic Cells"),
    list(cols = c("CD1cpositive_myeloid_dendritic_cell", "cd1cpositive_myeloid_dendritic_cell"),
         keep = "CD1cpositive_myeloid_dendritic_cell", label = "CD1c Myeloid Dendritic Cells"),
    list(cols = c("regulatory_T_cell", "regulatory_t_cell"),
         keep = "regulatory_T_cell", label = "Regulatory T Cells")
  ),
  Thymus = list(
    list(cols = c("Macrophages", "macrophage"),
         keep = "Macrophages", label = "Macrophages"),
    list(cols = c("Monocytes", "monocyte"),
         keep = "Monocytes", label = "Monocytes"),
    list(cols = c("Natural_killer_cells", "nk_cell"),
         keep = "Natural_killer_cells", label = "NK Cells"),
    list(cols = c("Plasma_cells", "plasma_cell"),
         keep = "Plasma_cells", label = "Plasma Cells")
  )
)

label_overrides <- c(
  memory_b_cell = "Memory B Cells"
)

all_meta_dt <- fread(all_meta_path)

tissue_files <- list.files(ts_data_dir, pattern = "_expression\\.csv$", full.names = FALSE)
tissues <- str_remove(tissue_files, "_expression\\.csv$")

format_celltype <- function(x) str_replace_all(x, "_", " ")
format_tissue <- function(x) str_replace_all(x, "_", " ")

build_nice_label_map <- function(tissue_name, all_meta_dt, merge_rules,
                                 label_overrides, ts_data_dir) {
  expr_path <- file.path(ts_data_dir, sprintf("%s_expression.csv", tissue_name))
  header <- names(fread(expr_path, nrows = 0))
  cell_types_raw <- setdiff(header, "Gene_name")

  placeholder <- grepl("^not.?exist", cell_types_raw, ignore.case = TRUE)
  cell_types <- cell_types_raw[!placeholder]

  tissue_meta <- all_meta_dt[Tissue == tissue_name]
  if (nrow(tissue_meta) == length(cell_types)) {
    nice <- setNames(tissue_meta$Cell_Type, cell_types)
  } else {
    nice <- setNames(format_celltype(cell_types), cell_types)
  }

  for (col_name in names(label_overrides)) {
    if (col_name %in% names(nice)) {
      nice[col_name] <- label_overrides[col_name]
    }
  }

  if (tissue_name %in% names(merge_rules)) {
    for (rule in merge_rules[[tissue_name]]) {
      if (rule$keep %in% cell_types) {
        nice[rule$keep] <- rule$label
      }
      others <- setdiff(rule$cols, rule$keep)
      cell_types <- setdiff(cell_types, others)
    }
    nice <- nice[cell_types]
  }

  nice
}

# Load WSS results and recompute -log10(p) with permutation floor
n_perm <- 10000
wss_dt <- fread(file.path(out_dir, "tabula_sapiens_wss_results.csv"))
wss_dt[, log10_p := -log10(pmax(p_value, 1 / (n_perm + 1)))]

sig_wss <- wss_dt[fdr < 0.01]
message(sprintf("Total FDR < 0.01: %d", nrow(sig_wss)))

if (nrow(sig_wss) == 0) {
  message("No significant results — skipping summary heatmap.")
  quit(status = 0)
}

# Build nice display labels
nice_lookup <- list()
for (tis in tissues) {
  lmap <- build_nice_label_map(tis, all_meta_dt, merge_rules,
                               label_overrides, ts_data_dir)
  for (ct in names(lmap)) {
    nice_lookup[[paste(tis, ct, sep = "|")]] <- lmap[ct]
  }
}

sig_wss[, display_label := {
  key <- paste(tissue, cell_type, sep = "|")
  nice <- nice_lookup[key]
  ifelse(
    vapply(nice, is.null, logical(1)),
    paste0(format_tissue(tissue), " — ", format_celltype(cell_type)),
    paste0(format_tissue(tissue), " — ", unlist(nice))
  )
}]

# Columns in the desired META display order (Glycemic -> Lpa); any clusters not
# in the list fall back to natural K order at the end.
desired_k_order <- c("K10", "K9", "K4", "K2", "K7", "K8", "K6", "K3", "K5", "K1")
present_clusters <- unique(sig_wss$cluster)
sig_clusters <- c(desired_k_order[desired_k_order %in% present_clusters],
                  sort(setdiff(present_clusters, desired_k_order)))
sig_labels <- sort(unique(sig_wss$display_label))

mat <- matrix(NA_real_, nrow = length(sig_labels), ncol = length(sig_clusters),
              dimnames = list(sig_labels, cluster_labels[sig_clusters]))

for (i in seq_len(nrow(sig_wss))) {
  rl <- sig_wss$display_label[i]
  cl <- cluster_labels[sig_wss$cluster[i]]
  mat[rl, cl] <- sig_wss$log10_p[i]
}

row_tissues <- str_extract(sig_labels, "^[^—]+") |> str_trim()
unique_tissues <- sort(unique(row_tissues))
tissue_colors <- setNames(
  hcl.colors(length(unique_tissues), palette = "Set 2"),
  unique_tissues
)

row_ha <- rowAnnotation(
  Tissue = row_tissues,
  col = list(Tissue = tissue_colors),
  show_legend = TRUE,
  annotation_legend_param = list(Tissue = list(title = "Tissue"))
)

max_val <- ceiling(max(mat, na.rm = TRUE))
max_val <- max(max_val, 1)
summary_col_fun <- colorRamp2(
  seq(0, max_val, length.out = 10),
  colorRampPalette(c("#FFF5EB", "#FDD0A2", "#FD8D3C", "#D94801", "#7F2704"))(10)
)

mat[is.na(mat)] <- 0

n_rows <- length(sig_labels)
panel_h <- max(6, n_rows * 0.25 + 2)

hm <- Heatmap(
  mat,
  col = summary_col_fun,
  na_col = "white",
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_order = sig_labels,
  column_order = cluster_labels[sig_clusters],
  row_names_side = "left",
  row_names_max_width = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(title = "-log10(p)", direction = "horizontal"),
  left_annotation = row_ha,
  name = "summary"
)

summary_path <- file.path(out_dir, "tabula_sapiens_summary_heatmap.png")
png(summary_path, width = 12, height = panel_h, units = "in", res = 600)
draw(hm,
     column_title = "Tabula Sapiens: FDR-Significant Cell Types (WSS, per-cluster FDR < 0.01)",
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     heatmap_legend_side = "bottom")
dev.off()
message(sprintf("Saved %s", basename(summary_path)))
