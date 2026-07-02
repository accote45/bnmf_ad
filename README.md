# multiancestry_polygenic

Private repository for investigating multiancestry polygenic mechanisms using Bayesian Non-negative Matrix Factorization (bNMF).

## Overview

This project decomposes GWAS effect sizes across multiple cardiometabolic traits into latent genetic clusters using bNMF, then evaluates cross-ancestry portability of these clusters. The workflow covers GWAS harmonization, bNMF factorization, PRS computation, and cross-ancestry Jaccard similarity testing.

## Repository Structure

```
├── Snakefile_a0                    # A0 pipeline: LDSC genetic correlations
├── Snakefile_a1                    # A1 pipeline: bNMF per ancestry + visualization
├── config/
│   ├── a0_config.yaml              # A0 config (LDSC traits, filters)
│   ├── a1_config.yaml              # A1 config (5-ancestry GWAS paths, bNMF params, QC)
│   ├── b1_config.yaml              # B1 config (genotypes, phenotypes, A1 results)
│   ├── b2_config.yaml              # B2 config (UKB DB, date fields, survival params)
│   ├── b3_1_config.yaml            # B3.1 config (biomarker traits, medication codes, corrections)
│   ├── b3_2_config.yaml            # B3.2 config (decile biomarker/risk plots)
│   └── config.yaml                 # Legacy pipeline config
├── profiles/lsf/                   # Snakemake LSF profile for HPC submission
├── scripts/
│   ├── a0_analysis/                # A0 LDSC genetic correlations + C+T PRS + cross-trait regression
│   ├── a1_analysis/                # A1 bNMF analysis (shared with C1)
│   ├── a1_3_analysis/              # A1.3 cross-ancestry GWAS signal consistency (EUR vs EAS)
│   ├── a2_analysis/                # A2 GTEx tissue expression enrichment
│   ├── b1_analysis/                # B1 cluster PRS association (8 scripts + supplementary tables; aou/ = All of Us replication kit)
│   ├── b2_analysis/                # B2 survival + B2.2 PheWAS (6 scripts)
│   ├── b3_1_analysis/              # B3.1 cluster PRS vs continuous biomarkers (4 scripts)
│   ├── b3_2_analysis/              # B3.2 decile biomarker/risk plots + risk stratification + 3D hyperplanes + cluster x GW-PRS interaction (4 scripts)
│   ├── gwas_processing/            # GWAS download, build checking, harmonization
│   ├── pipeline/                   # Legacy orchestrator + 1KG download
│   ├── analysis/                   # Figure 2/3 scripts (PRS, expression)
│   ├── compute_prs.R               # PRS computation (cluster-specific + total GRS)
│   └── wrangle_aou_phenotypes.R    # AoU phenotype extraction via BigQuery
├── sumstats/                       # Raw + harmonized GWAS summary statistics
├── reference/                      # 1KG LD reference panels (EUR, AFR, EAS)
└── results/                        # Analysis outputs per ancestry
```

## A0 Pipeline — LDSC Genetic Correlations

Computes pairwise genetic correlations between a reference trait (T2D, Suzuki 2024) and comparison traits using LD score regression, then converts observed-scale SNP heritabilities to the liability scale (Lee et al. 2011). Preprocessing includes MHC exclusion, HapMap3 SNP filtering, allele alignment, and strand-ambiguous removal.

```bash
snakemake -s Snakefile_a0 --profile profiles/lsf    # HPC with LSF
snakemake -s Snakefile_a0 -n                        # Dry run
```

### Output

`results/a0_analysis/rg_results.csv` — Tidy CSV: `trait1, trait2, rg, se, z, p, h2_obs_*, h2_int_*, gcov_int*, h2_liab_*`
`results/a0_analysis/liability_h2.csv` — Per-trait observed and liability-scale h² (with P, K, conversion factor)

### C+T PRS and cross-trait regression

Computes genome-wide C+T PRS (PRSice-2) for 6 traits in UKB, then runs cross-trait logistic regressions (T2D ~ each CVD PRS, and each CVD ~ T2D PRS). Uses UKB-free base GWAS (Mahajan-noukb / Tcheandjieu), independent of the LDSC pipeline's GWAS. `02_run_prsice.sh` takes an optional `KEEP_FILE` to restrict the target to a sample subset (currently the 478,636 ancestry-assigned cohort).

```bash
bsub -n 4 -R "rusage[mem=12000] span[hosts=1]" -W 4:00 -q premium -P acc_paul_oreilly \
  bash scripts/a0_analysis/02_run_prsice.sh LSF_ARRAY    # PRSice array (6 traits)
Rscript scripts/a0_analysis/03_cross_trait_regression.R   # cross-trait regressions
```

Output: `results/a0_analysis/prs_ct/prsice_output/{trait}.all_score`, `results/a0_analysis/prs_ct/cross_trait_regression_results.csv`

### T2D–CVD locus overlap

Clumps each GWAS to independent lead SNPs (plink, r²>0.05, 500 kb, P<5e-8) and quantifies how T2D loci overlap each CVD trait's loci.

```bash
bash scripts/a0_analysis/05a_clump_gwas_loci.sh prep
bash scripts/a0_analysis/05a_clump_gwas_loci.sh clump <trait>   # per trait (or LSF_ARRAY)
Rscript scripts/a0_analysis/05b_find_overlapping_loci.R         # pairwise overlap
Rscript scripts/a0_analysis/05c_write_loci_overlap_xlsx.R       # summary table
```

Output: `results/a0_analysis/loci_overlap/{overlap_pairwise,overlap_summary_cvd_side,overlap_all5}.csv`

## A1 Pipeline — bNMF Analysis

Runs bNMF factorization per ancestry on the published arm GWAS trait set, then generates multi-panel visualizations and cross-ancestry comparison. The published arm was validated as the optimal input via B1 comparison analyses (shared-trait and expanded arm sensitivity tests).

```bash
# Dry run
snakemake -s Snakefile_a1 -n

# HPC with LSF (5 ancestry jobs in parallel)
snakemake -s Snakefile_a1 --profile profiles/lsf --jobs 5

# Sensitivity clustering (k-means/hierarchical/half-max/DynamicTreeCut + bNMF-vs-hclust logistic), LSF
bash scripts/a1_analysis/submit_sensitivity.sh META
```

### Outputs

- `results/a1_analysis/{ancestry}/` — H/W matrices, QC reports, heatmaps, summaries
- `results/a1_analysis/{ancestry}/sensitivity/` — K-means/hierarchical/half-max/Dynamic Tree Cut assignments, profiles, heatmaps, dendrogram, silhouette plot, summary, logistic regression comparison results
- `results/a1_analysis/figures/` — `figure_main.png` (3 panels: heatmap, pos/neg bars, gene bars), per-ancestry gene bar plots (`panel_C_gene_bars_{ancestry}.png`), per-ancestry trait loading barplots (`barplot_{ancestry}_trait_loadings.png`), gene metric barplots (`barplot_{ancestry}_gene_{mean|max|specificity}.png`), gene lollipop and max-loading-vs-specificity scatter (`barplot_{ancestry}_gene_{lollipop|scatter}.png`), combined per-cluster trait+gene figures (`barplot_{ancestry}_combined_trait_gene.png` and the shared-trait split variant `barplot_{ancestry}_combined_trait_gene_split.png`), `heatmap_bnmf_vs_hclust_{ancestry}.png`, `cad_t2d_spectrum_{ancestry}.png`, `cross_ancestry_scatter.png`, `supp_bnmf_convergence.png`
- `results/a1_analysis/{ancestry}/cad_t2d_spectrum_scores_{ancestry}.tsv` — Per-cluster spectrum positions (z²-based h2_spectrum), donut proportions, GWS variant counts

## A1 Sensitivity Arms — GWAS Source & Trait Set Comparisons

Tests whether bNMF cluster structure is robust to GWAS source and trait selection. Five arms total, with pairwise comparisons:

```bash
# Run core + expanded arms and all comparisons
snakemake -s Snakefile_a1_shared --profile profiles/lsf
```

### Outputs

- `results/a1_analysis_core_published/EUR/` — Core published arm results
- `results/a1_analysis_core_splitukb/EUR/` — Core split-UKB arm results
- `results/a1_analysis_expanded/EUR/` — Expanded union arm results
- `results/a1_comparison_*/` — Pairwise comparison reports, cluster matching, heatmaps
- `results/b1_analysis/sensitivity_splitukb/` — B1 association results for split-UKB clusters

## B1 Analysis — Cluster PRS Association Testing

Tests whether A1 bNMF clusters are predictive of T2D/CAD susceptibility via individual-level PRS in UK Biobank across ancestries. Currently configured with META clusters (K=10) and META GWAS BETAs. Synthetic outcome = T2D AND CAD (intersection). All models adjusted for age, age², sex, genotyping batch, and PC1-10.

### Outputs (`results/b1_analysis/`)

- `forest_eur_validation.png` — Cluster PRS (per SD) vs genome-wide PRS, colored by T2D/CAD/T2D and CAD
- `forest_noneur_validation.png` — AFR/EAS/SAS faceted forest plot (per SD)
- `forest_combined_validation.png` — All four ancestries (EUR/AFR/EAS/SAS) faceted in one figure
- `forest_eur_top10pct.png` — Top 10% cluster PRS vs genome-wide PRS (top-decile model)
- `forest_noneur_top10pct.png` — Cross-ancestry top 10% forest plot
- `heatmap_cluster_prs_corr.png` — Cluster PRS Pearson correlation heatmap (upper triangle + diagonal; lower triangle blank)
- `association_results_all.csv` — Per-SD regression results (individual + joint × 3 outcomes × K clusters × 5 groups)
- `quantile_prs_associations.csv` — Quantile associations (top 20/10/5/1% × 10 clusters + 2 GW × 3 outcomes × 4 groups)
- `cross_disease_pgs_associations.csv` — Cross-disease results (T2D among CAD cases + CAD among T2D cases × 10 clusters × continuous + top 10% × 4 validation groups)
- `decile_event_rate.png` — 5×2 grid of T2D/CAD/T2D+CAD event rates by PGS decile per cluster
- `gw_threshold_selection.csv` — best-fit genome-wide PRS p-value threshold per trait, selected by Nagelkerke R² on eur_train (T2D `Pt_0.05`, CAD `Pt_1e-05`)
- `joint_vs_gw_comparison.csv` — R² and AUC for joint cluster PRS, GW T2D PRS, GW CAD PRS, and combined models (3 outcomes × 4 models)
- `joint_vs_gw_comparison.png` — Grouped bar chart comparing R² and AUC across models

Config: `config/b1_config.yaml` — includes `cluster_labels` block with META labels (Lpa, Metabolic Syndrome, Glucose 1, Lipid, ALP, Adiponectin-SHBG, Metabolic Syndrome-Inflammatory, Platelet, Glucose 2, Lipid-Liver)

### All of Us replication kit (`scripts/b1_analysis/aou/`)

Port of B1 for external validation in the All of Us Researcher Workbench. `make_aou_weights.R`
(repo root) builds one self-contained, GRCh38-lifted weights file to upload; `extract_acaf_weights.sh`
pulls the weight loci out of the AoU srWGS ACAF-threshold callset (FUSE-free, via `gsutil` byte-range
reads); `compute_prs_aou.sh` scores the extracted genotypes; `b1_aou_validation.ipynb` runs the
associations interactively in R (T2D/CAD via OMOP, forest plots across AoU ancestries).

```bash
Rscript scripts/b1_analysis/make_aou_weights.R       # -> results/b1_analysis/aou/b1_cluster_weights.tsv (upload this)
bash    scripts/b1_analysis/aou/extract_acaf_weights.sh # in AoU: -> acaf_extracted/ (small per-chr plink)
bash    scripts/b1_analysis/aou/compute_prs_aou.sh     # in AoU: GENO_DIR=acaf_extracted -> cluster_prs_all.tsv
```

Output: `b1_cluster_weights.tsv` (885 variants × cluster weights), then `cluster_prs_all.tsv` and
`association_results_*.csv` + forest PNGs inside the AoU workspace.

## B1.2 Analysis — Pathway PRS vs bNMF Cluster PRS Comparison

Sensitivity analysis comparing data-driven bNMF cluster PRS (from B1) to experimentally-defined pathway PRS from PRSet. Meta-analyzes T2D + CAD META GWAS via METAL IVW, runs PRSet with MSigDB c2 pathways (~5,878), then compares the top 6 pathway PRS (by Synthetic OR/SD on validation) to the 9 META bNMF cluster PRS.

Config: `config/b1_2_config.yaml`

## B2 Analysis — Survival Analysis (Time-to-Comorbidity)

Survival analysis testing whether bNMF cluster PRS predicts time from first T2D/CAD diagnosis to development of the other condition (comorbidity).

Scripts 01-03 are shared between B2 and B2.1 via `--config` flag.

### Outputs (`results/b2_analysis/`)

- `cox_results.csv` — 18 Cox model results (3 directions × 6 clusters)
- `cumhaz_cad_to_t2d.png`, `cumhaz_t2d_to_cad.png`, `cumhaz_combined.png` — Cumulative hazard curves
- `cluster_prevalence_{train,validation}.tsv` — Comorbidity prevalence per cluster
- `cohort_tabulation.tsv` — CAD-first/T2D-first/dual counts

Config: `config/b2_config.yaml`

## B2.1 Analysis — Simplified Survival Analysis (Weighted + Hard-Assignment PRS)

Simplified survival analysis pipeline that runs Cox proportional hazards models and adjusted cumulative hazard plots for 10 META cluster PRS at 3 percentile thresholds (top 33%, 20%, 10%). Outcome: T2D/CAD composite — time from first diagnosis to developing the other condition.

Supports two PRS modes via `--hard-assignment` flag:
- **Default (soft assignment)**: Each variant contributes to all cluster PRS weighted by W matrix
- **Hard assignment**: Each variant assigned to only its max-W cluster, zeroing out all others

```bash
# Default (soft-assignment PRS)
bash scripts/b2_analysis/submit_b2_1_pipeline.sh

# Hard-assignment PRS
bash scripts/b2_analysis/submit_b2_1_pipeline.sh --hard-assignment
```

### Outputs

| Mode | Results dir | Key outputs |
|------|-------------|-------------|
| Default | `results/b2_1_analysis/` | `cox_results.csv`, `cumhaz_top{33,20,10}.png`, `cumhaz_top{33,10}_{t2d_to_cad,cad_to_t2d}.png` (direction-split), `diagnosis_interval_summary.csv`, `age_at_diagnosis_summary.csv` |
| Hard assignment | `results/b2_1_analysis_hard_assignment/` | Same file names |

Config: `config/b2_1_config.yaml` — 10 clusters (K1-K10 with canonical META labels), 3 thresholds (0.67, 0.80, 0.90)

## B2.2 Analysis — PheWAS (Cluster PRS vs ICD10 Phenotypes)

Tests whether META cluster PRS are associated with a broad range of ICD10-coded disease phenotypes via logistic regression. Analyzes 473 ICD10 codes across 10 disease chapters (E: Endocrine through N: Genitourinary) against all 10 META cluster PRS in EUR validation samples.

Supports `--hard-assignment` flag (same pattern as B1/B2.1) to run on hard-assignment PRS with output to `results/b2_2_analysis/prs_hard_assignment/`.

```bash
# Default (soft-assignment PRS)
bash scripts/b2_analysis/submit_b2_2_pipeline.sh

# Hard-assignment PRS
bash scripts/b2_analysis/submit_b2_2_pipeline.sh --hard-assignment
```

### Outputs

| Mode | Results dir | Key outputs |
|------|-------------|-------------|
| Default | `results/b2_2_analysis/` | `phewas_results.csv`, `manhattan_K{1-10}.png`, `manhattan_combined.png`, `manhattan_2x2_selected.png` (Lpa, Glycemic, Blood Pressure-Stature, Adiponectin) |
| Hard assignment | `results/b2_2_analysis/prs_hard_assignment/` | Same file names |

Config: `config/b2_2_config.yaml` — 10 clusters (K1-K10), 10 ICD10 chapters (E-N), min 50 cases

## B3.1 Analysis — Cluster PRS vs Continuous Biomarker Traits

Tests associations between the 10 META bNMF cluster PRS and 14 continuous biomarker/anthropometric traits in UKB EUR validation via linear regression, with sex-stratified analyses. Extends B1 (binary outcomes) to continuous biological signatures. Medication corrections applied to avoid confounding by treatment effects.

### Outputs (`results/b3_1_analysis/`)

- `biomarker_phenotypes.tsv` — Extracted + corrected phenotype values (includes PCE_ASCVD, PREVENT_ASCVD)
- `medication_flags.tsv` — Per-individual medication status (statin, antihypertensive, diabetes med, has_diabetes, current_smoker)
- `phenotype_summary.tsv` — Descriptive statistics per trait
- `biomarker_association_all.csv` — Combined individual + joint regression results with `sex_group` column (All/Male/Female)
- `forest_K*_*.png` — 10 per-cluster forest plots (sex-stratified, category-faceted)
- `forest_combined.png` — Combined grid: trait categories (rows) × clusters (columns)
- `heatmap_cluster_biomarker.png` — Summary heatmap (10 clusters × traits, All only)
- `forest_combined_selected.png` — Selected 4-cluster forest plot (Lpa, Glycemic, Blood Pressure-Stature, Adiponectin); pass `--selected-clusters` to choose
- `forest_individual_vs_joint.png` — FDR-significant associations: individual vs joint model comparison

Config: `config/b3_1_config.yaml`

## B3.2 Analysis — Decile Biomarker/Risk Plots & Risk Stratification

Visualizes the dose-response relationship between cluster PGS deciles and biomarker/clinical risk values, and demonstrates that cluster PGS provides additional risk stratification within standard clinical risk categories.

### Outputs (`results/b3_2_analysis/prs_hard_assignment/`)

- `decile_biomarker_K*_*.png` — Per-cluster 2×2 biomarker/risk plots (10 files)
- `decile_biomarker_GW_T2D.png`, `decile_biomarker_GW_CAD.png` — Genome-wide PRS sanity checks
- `decile_biomarker_Combined_Cluster_PGS.png` — Combined cluster PGS (T2D+CAD joint model) decile plot
- `risk_stratification_PCE.png` — PCE risk stratification by cluster PGS (4 panels)
- `risk_stratification_PREVENT.png` — PREVENT risk stratification by cluster PGS (4 panels)
- `risk_stratification_PCE_combined.png` — PCE risk stratification by combined cluster PGS (single panel)
- `risk_stratification_PREVENT_combined.png` — PREVENT risk stratification by combined cluster PGS (single panel)
- `hyperplane_3d_K*_*.png` — Per-cluster 3D GAM surface plots (PGS percentile ~ PCE + HbA1c)
- `forest_top_bottom_decile.png` — Top/bottom 20% mean ± SD forest plot (PCE + HbA1c, vermillion/cerulean)
- `hyperplane_3d_combined_2x2.png` — Combined 2×2 figure (Lipid, Glucose 1, MetSyn-Inflammatory, Lipid-Liver)
- `decile_line_cluster_gw_prs.png` — Cluster PGS decile vs biomarker stratified by GW PRS: Glucose 2 × GW T2D PRS → HbA1c (cerulean), Lipid × GW CAD PRS → PCE risk (vermillion)
- `decile_line_combined_score.png` — GW PGS alone vs equal-weight combined (GW PGS + cluster cPGS) score; 2×2 with left = HbA1c (T2D + Glycemic), right = PCE (CAD + Adiponectin); top row = score decile, bottom row = cumulative top 10%→1%. Points annotated with the relative % difference in individuals above the clinical threshold (HbA1c 5.7%, PCE 7.5%). Pass `--combined`
- `decile_line_cad_cluster_grid.png` — 5×2 grid of mean PCE risk by score decile, GW CAD PGS alone vs GW CAD PGS + each META cluster cPGS (one panel per cluster K1–K10), shared y-axis. Pass `--cad-cluster-grid`
- `decile_line_t2d_cluster_grid.png` — 5×2 grid of mean HbA1c by score decile, GW T2D PGS alone vs GW T2D PGS + each META cluster cPGS (one panel per cluster K1–K10), shared y-axis. Pass `--t2d-cluster-grid`
- `decile_line_hba1c_pce.png` — All 10 cluster PGS deciles vs mean HbA1c and PCE risk (one line per cluster)

Config: `config/b3_2_config.yaml`

## A1.2 Analysis — Published Cluster Comparison

Compares the 10 META bNMF clusters against clusters from three published studies using Pearson correlation of trait-weight profiles. Comparator studies: Suzuki 2024 (Nature, 8 hard clusters), Smith 2024 (Nature Medicine, 12 soft clusters), Pascat 2026 (Nature Communications, 5 hard clusters).

### Execution

```bash
bash scripts/a1_2_analysis/submit_a1_2_cluster_comparison.sh
```

### Outputs (`results/a1_2_analysis/`)

- `heatmap_pearson_{suzuki2024,smith2024,pascat2026}.png` — Individual correlation heatmaps
- `combined_heatmap_published_comparison.png` — Multi-panel combined figure
- `correlation_summary.csv` — All pairwise r and p-values
- `best_match_summary.csv` — Best match (max |r|) per META cluster per study

## A1.3 Analysis — Cross-ancestry GWAS Signal Consistency (EUR vs EAS)

Compares ancestry-specific CAD and T2D GWAS arms (CAD: Aragam2022 EUR vs Sakaue2021 EAS; T2D: Suzuki2024 EUR vs EAS) to assess signal consistency. Per trait: defines independent lead SNPs per arm (ancestry-matched LD clumping), quantifies shared vs ancestry-unique loci (±500kb), regresses EUR beta on EAS beta at EAS lead SNPs (R² + 95% CI, with allele-orientation QC), and compares allele frequencies. EUR = reference; trans-ancestry META_ files excluded.

### Execution

```bash
bsub -P acc_paul_oreilly -q premium -n 4 -R "rusage[mem=32000] span[hosts=1]" -W 4:00 \
  "module load R/4.2.0 && module load plink/1.90b6.21 && Rscript scripts/a1_3_analysis/01_cross_ancestry_consistency.R --trait both"
```

### Outputs (`results/a1_3_analysis/`)

- `cross_ancestry_summary.csv` — Lead counts, locus sharing, beta slope/R²/CI, allele-orientation flag
- `lead_snps_{CAD,T2D}_{EUR,EAS}.csv`, `locus_sharing_{CAD,T2D}.csv`
- `barplot_locus_sharing_{CAD,T2D}.png` — EUR-only / shared / EAS-only loci
- `scatter_beta_{CAD,T2D}.png`, `scatter_AF_T2D.png` — effect-size and allele-frequency concordance

## A2 Analysis — GTEx Tissue Expression Enrichment

Tests whether genes near cluster variants are preferentially expressed in specific tissues compared to background protein-coding genes. Uses one-sided Welch t-test (cluster > control) across 46 GTEx tissues. Each cluster heatmap has 2 columns: Absolute Expression (log2(TPM+1), labeled "Abs.") and Relative Expression (Skene specificity, labeled "Rel."). The Absolute test restricts the background to genes expressed (>0) in each tissue/cell type to avoid trivial saturation against silent genes; Relative uses the full background. GTEx uses a green color palette; Tabula Sapiens uses an orange color palette.

### Execution

```bash
# Submit to LSF (both scripts are independent and can run in parallel)
bsub -P acc_paul_oreilly -q premium -R "rusage[mem=16000]" -W 1:00 \
  Rscript scripts/a2_analysis/a2_gtex_heatmaps.R
bsub -P acc_paul_oreilly -q premium -R "rusage[mem=16000]" -W 1:00 \
  Rscript scripts/a2_analysis/a2_tabula_sapiens_heatmaps.R
```

### Outputs (`results/a2_analysis/`)

- `gtex_ttest_results_nearest.csv` — 920 rows (10 clusters × 46 tissues × 2 expression types)
- `gene_mapping_summary.csv` — Per-cluster gene counts, mean/SD mapping distance
- `gtex_heatmap_CH_K{1-10}.png` — Individual per-cluster heatmaps (2-column ComplexHeatmap, green palette)
- `gtex_heatmap_combined_CH_part1.png` — K1–K6 combined grid (2×3, panels a–f)
- `gtex_heatmap_combined_CH_part2.png` — K7–K10 combined grid (2×2, panels g–j)

#### Single-gene-set mode (sanity check)

`a2_gtex_heatmaps.R` accepts an optional external gene list to run the same Abs/Rel
heatmap on one gene set (e.g. reproducing a published trait's nearest-gene heatmap).
`sanity_check_t2d_clump_nearest.R` regenerates a trait's nearest-to-hit genes by
clumping a GWAS (reuses `prep_bnmf.R`'s `ld_clump`).

```bash
Rscript scripts/a2_analysis/a2_gtex_heatmaps.R \
  --gene_list_file <ENSG.csv> --label T2D --out_dir results/a2_analysis/sanity_check
```

Output (`results/a2_analysis/sanity_check/`): `gtex_heatmap_CH_<label>.png`,
`gtex_ttest_results_nearest.csv`, plus clumping/concordance CSVs.

#### Standalone legends

`a2_make_legends.R` renders the horizontal "-logP" color-bar legends (GTEx green,
Tabula Sapiens orange, catlas purple) for manual figure assembly, since the heatmaps
suppress their own legend.

```bash
Rscript scripts/a2_analysis/a2_make_legends.R
```

Output (`results/a2_analysis/legends/`): `{gtex,tabula_sapiens,catlas}_legend.png` + combined.

### Tabula Sapiens Cell-Type Expression

Tests cell-type-level enrichment across 14 tissues (Blood, Fat, Heart, Large Intestine, Liver, Lung, Lymph Node, Mammary, Prostate, Skin, Small Intestine, Spleen, Thymus, Vasculature) using single-cell RNA-seq from the Tabula Sapiens atlas.

#### Outputs (`results/a2_analysis/tabula_sapiens/`)

- `tabula_sapiens_{tissue}_K{1-10}.png` — Per-tissue × per-cluster 2-column heatmaps (orange palette)
- `tabula_sapiens_combined_{tissue}_part{1,2}.png` — Combined grid figures for 6 focus tissues (12 files total)
- `tabula_sapiens_permtest_results.csv` — Per-tissue permutation test results (nearest genes only)
- `tabula_sapiens_wss_results.csv` — WSS enrichment results with per-cluster FDR (3,070 rows)
- `tabula_sapiens_summary_heatmap.png` — FDR-significant cell types across all 10 clusters (orange palette)

### CATLAS cCRE Enrichment (Pancreas + Esophagus)

Tests whether genes near cluster variants are enriched for candidate cis-regulatory elements (cCREs) from CATLAS single-cell ATAC-seq data. Uses gene-level cCRE counts within ±50kb windows as a regulatory activity proxy, with 10k permutation tests. Purple color palette. CATLAS BED files (hg38) are lifted to hg19 via liftOver.

```bash
# Download + liftOver
bash scripts/a2_analysis/download_catlas_pancreas.sh
bash scripts/a2_analysis/download_catlas_esophagus.sh

# Per-tissue heatmaps
bsub -P acc_paul_oreilly -q premium -R "rusage[mem=16000]" -W 1:00 \
  "module load R/4.2.0 && Rscript scripts/a2_analysis/a2_catlas_pancreas_heatmaps.R"
bsub -P acc_paul_oreilly -q premium -R "rusage[mem=16000]" -W 1:00 \
  "module load R/4.2.0 && Rscript scripts/a2_analysis/a2_catlas_esophagus_heatmaps.R"

# Combined WSS summary
bsub -P acc_paul_oreilly -q premium -R "rusage[mem=16000]" -W 1:00 \
  "module load R/4.2.0 && Rscript scripts/a2_analysis/a2_catlas_wss_summary.R"
```

#### Outputs

- `results/a2_analysis/catlas/catlas_pancreas/` — Per-cluster heatmaps (11 cell types), combined grids, CSV
- `results/a2_analysis/catlas/catlas_esophagus/` — Per-cluster heatmaps (5 cell types), combined grids, CSV
- `results/a2_analysis/catlas/catlas_summary/catlas_wss_results.csv` — WSS enrichment with per-cluster BH FDR (160 rows)
- `results/a2_analysis/catlas/catlas_summary/catlas_summary_heatmap.png` — FDR-significant cell types (purple palette, teal/orange tissue strip)

#### CATLAS Firth Enrichment (supplementary, SNP-level)

Additive SNP-level alternative to the permutation analysis above. Per cluster, lead SNPs (Y=1) are compared against matched null SNPs (Y=0; 1000G EUR variants within ±50kb of a lead, r²<0.05 to all cluster leads) via Firth logistic regression adjusting for genic context (CDS/5'UTR/3'UTR, GENCODE v19), testing direct overlap with each cell type's ATAC peaks. Bonferroni per cluster. CATLAS-only (requires interval data; not applicable to GTEx/Tabula Sapiens expression). Run under the rocky9 R 4.2.0 recipe.

```bash
Rscript scripts/a2_analysis/a2_catlas_firth_enrichment.R
```

Output:
- `results/a2_analysis/catlas/catlas_firth/catlas_firth_enrichment_results.csv` + combined heatmap
- `results/a2_analysis/catlas/catlas_{pancreas,esophagus}/firth/` — per-cluster heatmaps + combined grid