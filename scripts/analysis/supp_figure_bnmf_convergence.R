# supp_figure_bnmf_convergence.R
# Supplementary figure: scatterplot of iterations vs K_converged.
# X-axis: iterations, Y-axis: K_converged, color by ancestry, fill by is_optimal_k.

library(data.table)
library(ggplot2)

source("scripts/a1_analysis/figure_utils.R")

# 1. Load & combine run_summary.csv for each ancestry
ancestries <- c("META", "EUR", "AFR")
dt <- rbindlist(lapply(ancestries, function(anc) {
  d <- fread(file.path("results", anc, sprintf("run_summary_%s.csv", anc)))
  d[, ancestry := anc]
}))

# 2. Prepare types for aesthetics
dt[, ancestry      := factor(ancestry, levels = ancestries)]
dt[, is_optimal_k  := as.logical(is_optimal_k)]

# 3. Plot: iterations vs K_converged, color=ancestry, fill=is_optimal_k
p <- ggplot(dt, aes(x = iterations, y = K_converged, color = ancestry,
                     fill = is_optimal_k, shape = is_optimal_k)) +
  geom_jitter(size = 3, stroke = 1, width = 0, height = 0.1) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Yes", "FALSE" = "No"),
                     name   = "Optimal K") +
  scale_fill_manual(values = c("TRUE" = "black", "FALSE" = NA),
                    labels = c("TRUE" = "Yes", "FALSE" = "No"),
                    name   = "Optimal K", guide = "none") +
  scale_y_continuous(breaks = seq(0, max(dt$K_converged))) +
  labs(x = "Iterations", y = expression(K[converged]), color = "Ancestry") +
  theme_big_text() +
  theme(legend.position = "right")

# 4. Save
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/supp_figure_bnmf_convergence.png", p,
       width = 8, height = 6, dpi = 300)

cat("Saved: results/figures/supp_figure_bnmf_convergence.png\n")
