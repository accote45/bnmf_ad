#!/usr/bin/env Rscript
# Tabulate number of SNVs per META bNMF cluster
# Assigns each variant to its max-weight cluster from the W matrix

library(tidyverse)

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
w_file   <- file.path(base_dir, "results/a1_analysis/META/W_matrix_META.tsv")
out_file <- file.path(base_dir, "results/a1_analysis/META/snv_counts_per_cluster_META.csv")

# Load W matrix and assign each variant to its max-weight cluster
w_mat <- read_tsv(w_file, show_col_types = FALSE)

k_cols <- str_subset(colnames(w_mat), "^K\\d+$")

snv_counts <- w_mat %>%
  mutate(cluster = k_cols[max.col(pick(all_of(k_cols)))]) %>%
  count(cluster, name = "n_snvs") %>%
  arrange(as.integer(str_extract(cluster, "\\d+")))

print(snv_counts)

write_csv(snv_counts, out_file)
cat("\nWritten to:", out_file, "\n")
