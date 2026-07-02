#!/usr/bin/env Rscript
# Diagnostic: distribution of -log10(p) values in FDR < 0.01 WSS results

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
out_dir <- file.path(base_dir, "results/a2_analysis/tabula_sapiens")

wss_dt <- fread(file.path(out_dir, "tabula_sapiens_wss_results.csv"))
sig <- wss_dt[fdr < 0.01]

# Show the raw p-value and the -log10(p) as computed in the heatmap
sig[, log10p_raw := -log10(pmax(p_value, 1e-300))]

message(sprintf("FDR < 0.01: %d results", nrow(sig)))
message(sprintf("  p = 0: %d", sum(sig$p_value == 0)))
message(sprintf("  p > 0: %d", sum(sig$p_value > 0)))
message(sprintf("  -log10(p) range: %.1f to %.1f", min(sig$log10p_raw), max(sig$log10p_raw)))

p <- ggplot(sig, aes(x = log10p_raw)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(
    title = "Distribution of -log10(p) for FDR < 0.01 WSS results",
    subtitle = sprintf("n = %d results; %d with p = 0 (-log10 = 300)",
                        nrow(sig), sum(sig$p_value == 0)),
    x = "-log10(p-value)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)

ggsave(file.path(out_dir, "pval_diagnostic.png"), p,
       width = 8, height = 5, dpi = 300)
message("Saved pval_diagnostic.png")
