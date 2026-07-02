library(tidyverse)

# ── Parameters ──────────────────────────────────────────────────────────────────
base_dir    <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
pheno_file  <- file.path(base_dir, "results/a0_analysis/prs_ct/phenotypes_combined.txt")
output_file <- file.path(base_dir, "phenotypes/phenotype_prevalence.csv")

# ── Load data ───────────────────────────────────────────────────────────────────
pheno <- read_tsv(pheno_file, show_col_types = FALSE)

# ── Compute prevalence ──────────────────────────────────────────────────────────
prevalence <- tibble(
  phenotype = c("T2D", "CAD", "T2D_CAD"),
  n_cases = c(
    sum(pheno$T2D == 1),
    sum(pheno$CAD == 1),
    sum(pheno$T2D == 1 & pheno$CAD == 1)
  ),
  n_total = nrow(pheno)
) %>%
  mutate(
    n_controls    = n_total - n_cases,
    prevalence_pct = round(n_cases / n_total * 100, 2)
  ) %>%
  select(phenotype, n_cases, n_controls, n_total, prevalence_pct)

print(prevalence)

# ── Write output ────────────────────────────────────────────────────────────────
write_csv(prevalence, output_file)
cat("\nWritten to:", output_file, "\n")
