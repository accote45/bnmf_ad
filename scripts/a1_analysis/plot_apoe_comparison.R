#!/usr/bin/env Rscript
# plot_apoe_comparison.R
# Before/after heatmap for the APOE-exclusion decision: the with-APOE bNMF (K=3)
# has an amyloid (CSF Abeta42) cluster driven by APOE; removing APOE (K=4)
# dissolves it and reveals immune / lipid / metabolic / neurodegeneration subtypes.
#
# One stacked faceted heatmap (net weight = pos - neg per cluster x trait), so the
# CSF Abeta42 column can be compared top (with APOE) vs bottom (without APOE).
#
# Usage:
#   Rscript scripts/a1_analysis/plot_apoe_comparison.R \
#     [--withapoe <H_matrix.tsv>] [--noapoe <H_matrix.tsv>] [--out <png>]

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- which(args == f); if (length(i)) args[i + 1] else d }
withapoe <- getarg("--withapoe", "results/ad_analysis/META_withAPOE/H_matrix_META.tsv")
noapoe   <- getarg("--noapoe",   "results/ad_analysis/META/H_matrix_META.tsv")
out      <- getarg("--out",      "figures/apoe_comparison_heatmap.png")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

clean_trait <- function(x) gsub("_", " ", sub("_[A-Za-z]+\\d{4}$", "", x))

# Cluster labels per run (keyed to each run's K-numbering; verified from the
# summaries: the with-APOE K3 is CSF-Abeta42-dominated -> "Amyloid (APOE)").
lab_withapoe <- c(K1 = "Neurodegeneration", K2 = "Immune", K3 = "Amyloid (APOE)")
lab_noapoe   <- c(K1 = "Neurodegeneration", K2 = "Immune", K3 = "Metabolic", K4 = "Lipid")

net_long <- function(h_path, analysis, labels) {
  h <- read_tsv(h_path, show_col_types = FALSE)
  trait_cols <- setdiff(names(h), "Cluster")
  traits <- unique(sub("_(pos|neg)$", "", trait_cols))
  df <- bind_rows(lapply(traits, function(tr) {
    pos <- paste0(tr, "_pos"); neg <- paste0(tr, "_neg")
    v <- (if (pos %in% names(h)) h[[pos]] else 0) - (if (neg %in% names(h)) h[[neg]] else 0)
    tibble(Cluster = h$Cluster, trait = tr, net = v)
  }))
  df %>% mutate(analysis = analysis,
                cluster_label = ifelse(Cluster %in% names(labels), labels[Cluster], Cluster))
}

d1 <- net_long(withapoe, "With APOE  (K=3)", lab_withapoe)
d2 <- net_long(noapoe,   "Without APOE  (K=4)", lab_noapoe)

# quick verification to console: which cluster carries CSF Abeta42?
csf <- bind_rows(d1, d2) %>% filter(grepl("CSF_Abeta42", trait))
cat("CSF Abeta42 net weight by cluster:\n"); print(csf %>% select(analysis, cluster_label, net))

dat <- bind_rows(d1, d2) %>%
  mutate(analysis = factor(analysis, levels = c("With APOE  (K=3)", "Without APOE  (K=4)")),
         cluster_label = factor(cluster_label,
           levels = c("Amyloid (APOE)", "Lipid", "Metabolic", "Immune", "Neurodegeneration")))

# Trait order: grouped, with CSF Abeta42 highlighted at the end of the neuro block.
trait_order <- c("HDL","Triglycerides","ApoA","ApoB","LDL","TotalCholesterol",
                 "BMI","Hba1c","RandomGlucose","T2D",
                 "WBC","LymphocyteCount","MonocyteCount","NeutrophilCount","Eosinophilcount","CRP",
                 "SBP","DBP",
                 "ALS","Parkinsons","CSF_Abeta42",
                 "Bipolar","MDD","Neuroticism","Schizophrenia")
present <- intersect(trait_order, unique(dat$trait))
present <- c(present, setdiff(unique(dat$trait), present))
dat$trait <- factor(dat$trait, levels = present)

# Cap the scale so the extreme amyloid cell (~-6) saturates but the subtype
# signals (~+/-3) stay readable; the true value is off-scale by design.
cap <- 4

p <- ggplot(dat, aes(trait, cluster_label, fill = net)) +
  geom_tile(color = "white", linewidth = 0.5) +
  # outline the CSF Abeta42 cells (the APOE-driven amyloid signal) in both panels
  geom_tile(data = function(d) subset(d, trait == "CSF_Abeta42"),
            color = "#111111", linewidth = 1.1) +
  scale_fill_gradient2(low = "steelblue4", mid = "white", high = "firebrick3",
                       midpoint = 0, limits = c(-cap, cap), oob = scales::squish,
                       name = "Net weight\n(pos - neg)") +
  scale_x_discrete(labels = clean_trait) +
  facet_wrap(~ analysis, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = NULL,
       title = "Excluding APOE dissolves the amyloid cluster",
       subtitle = "Net cluster-trait weights. The CSF Abeta42 column is the APOE-driven amyloid signal.") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 13),
        plot.title = element_text(face = "bold"))

ggsave(out, p, width = 12, height = 7, dpi = 300, bg = "white")
cat("Wrote", out, "\n")
