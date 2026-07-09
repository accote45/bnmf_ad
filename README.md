# bnmf_ad

Discovering **genetic subtypes of Alzheimer's disease (AD)** by decomposing GWAS
effect sizes across an AD-relevant trait panel with Bayesian Non-negative Matrix
Factorization (bNMF).

This is a retarget of a cardiometabolic (CAD/T2D) bNMF pipeline to AD. See
[AD_PIPELINE_SETUP.md](AD_PIPELINE_SETUP.md) for the full setup, data provenance,
and Minerva bring-up notes.

## Concept

The pipeline builds a **variant × trait Z-score matrix** and factorizes it with bNMF:

- **Rows = variant universe.** The genome-wide-significant, MAF-filtered,
  LD-clumped lead SNPs from the *reference* AD GWAS (127 PGC3 full-cohort loci).
- **Columns = trait panel.** A set of AD-relevant quantitative-trait GWAS; each
  row variant's Z-score is pulled from each trait GWAS.
- **bNMF → K clusters** = latent genetic subtypes, each defined by a trait
  signature (W = variants × K, H = K × traits).

The **primary analysis excludes the APOE region** (chr19 gene ±1 Mb) so subtypes
reflect non-APOE genetic structure; otherwise APOE's outsized effect dominates the
factorization and anchors an amyloid cluster. With APOE excluded the v1 run
resolves **4 subtypes**: Neurodegeneration (FTD/ALS), Immune, Metabolic, Lipid.

## Scope

**v1 = trans-ancestry / META.** The 127 AD loci and the trait GWAS are both
multi-ancestry, so the panel uses the `.META.` version of every trait; the LD
reference is EUR 1000G (the standard pragmatic approximation for a trans-ancestry
set). Run label = `META`. A single ancestry means **Snakemake is not required** —
the R scripts are run directly.

## Repository Structure

```
├── config/
│   └── ad_config.yaml              # The AD run config: ref/trait GWAS paths, bNMF params, APOE exclusion, cluster labels
├── scripts/
│   ├── gwas_processing/            # GWAS harmonization → *.processed.txt.gz (VAR_ID, aligned alleles, Z/BETA/SE, GRCh37)
│   │   ├── harmonize_sumstats.py   #   genome-wide harmonizer (+ check_build.py, column_map.py)
│   │   ├── convert_ad_loci.R       #   127-loci AD reference → harmonized reference file (variant universe)
│   │   ├── harmonize_ad_traits.sh  #   harmonize the 7 new AD-specific trait GWAS (+ prep_ad_trait_gwas.py)
│   │   └── …                       #   (legacy per-trait harmonize batch scripts from the cardiometabolic build)
│   ├── a1_analysis/                # bNMF engine + visualization
│   │   ├── 01_run_bnmf.R           #   QC + clump + Z-matrix + proxy + allele-align + bNMF
│   │   ├── bnmf_algorithm.R        #   bNMF math (unchanged engine)
│   │   ├── prep_bnmf.R             #   QC / clumping / Z-matrix / allele-alignment helpers
│   │   ├── 02_visualize_results.R  #   W/H heatmaps + trait-signature figures
│   │   ├── plot_apoe_comparison.R  #   with-APOE vs no-APOE before/after heatmap
│   │   └── make_concept_schematic.R#   presentation schematic of the pipeline
│   └── setup/
│       └── link_reuse_sumstats.sh  # symlink the 19 reused cardiometabolic .META. trait GWAS into sumstats/harmonized/
├── sumstats/                       # Raw + harmonized GWAS summary statistics (git-ignored; lives on Minerva)
├── reference/                      # 1000G EUR LD reference panel (git-ignored; symlinked on Minerva)
└── results/ad_analysis/            # bNMF outputs + figures
```

> Note: `scripts/` also retains code from the original cardiometabolic pipeline
> (e.g. `a0/a2/b1/b2/b3` analyses, extra harmonize batch scripts). These are **not
> part of the AD analysis** and are kept only as reference; the AD run uses the
> paths documented here.

## Environment

Runs on **Minerva (Mount Sinai LSF)**. Code is edited locally, pushed to GitHub,
pulled on Minerva, and run there.

- **R:** `module load R/4.2.0` + `module load plink/1.90b6.21`.
- **Harmonization:** a `gwaslab` conda env (python 3.10, gwaslab 4.1.9). Build
  recipe and the one-line import-bug patch are documented in the header of
  `scripts/gwas_processing/harmonize_ad_traits.sh`. The batch script self-activates
  the env.
- Large data (AD GWAS, trait GWAS, 1000G panels) lives under `sumstats/` and
  `reference/` — both git-ignored, never committed.

## Pipeline

Config: [config/ad_config.yaml](config/ad_config.yaml) — reference/trait GWAS
paths, bNMF parameters (`K=10` initial, `nreps=50`, `p_threshold=5e-8`,
`maf_threshold=0.001`), APOE region exclusion, and per-run cluster labels.

### 1. Reference GWAS (variant universe / matrix rows)

Convert the 127-loci AD "top loci" table into the harmonized reference format:

```bash
Rscript scripts/gwas_processing/convert_ad_loci.R
#   → sumstats/harmonized/AD_PGC3_top127.META.GRCh37.processed.txt.gz
```

### 2. Trait panel (matrix columns)

25 traits total: 18 reused cardiometabolic `.META.` GWAS (lipid/APOE axis, immune,
glycemic, vascular) plus 7 new AD-specific GWAS (ALS, Parkinson's, CSF Aβ42,
Bipolar, MDD, Neuroticism, Schizophrenia).

```bash
# Symlink the 19 reused (already-harmonized) trait GWAS
bash scripts/setup/link_reuse_sumstats.sh

# Harmonize the 7 new AD-specific trait GWAS
bsub < scripts/gwas_processing/harmonize_ad_traits.sh
#   → sumstats/harmonized/*.processed.txt.gz
```

### 3. Run bNMF + visualize

```bash
module load R/4.2.0 && module load plink/1.90b6.21

Rscript scripts/a1_analysis/01_run_bnmf.R \
  --config config/ad_config.yaml --ancestry META

Rscript scripts/a1_analysis/02_visualize_results.R \
  --config config/ad_config.yaml
```

`01_run_bnmf.R` accepts an optional `--select_k` override to force a specific K
instead of the modal rank (K sensitivity analysis).

### Outputs (`results/ad_analysis/`)

- `META/` — H/W matrices (`H_matrix_META.tsv`, `W_matrix_META.tsv`), QC reports,
  cluster summaries
- `META_withAPOE/` — the with-APOE comparison run (K=3), used by the APOE
  before/after figure
- `figures/` — bNMF heatmaps, trait-signature barplots, `apoe_comparison_heatmap.png`
  (with- vs without-APOE), `concept_schematic.png`

## Trait panel

| Axis | Traits (`.META.` unless noted) |
|------|--------------------------------|
| Lipid / APOE | HDL, LDL, TotalCholesterol, Triglycerides (Graham 2021); ApoA, ApoB (Karczewski 2025) |
| Immune / inflammatory | CRP, Lymphocyte, Monocyte, Neutrophil, WBC, Eosinophil (Karczewski 2025) |
| Metabolic / glycemic | T2D (Suzuki 2024), RandomGlucose (Lagou 2023), HbA1c (Verma 2024), BMI (Karczewski 2025) |
| Vascular | SBP, DBP (Karczewski 2025) |
| AD-specific | ALS (van Rheenen 2021), Parkinson's (Kim 2023), CSF Aβ42 (Timsina 2026, EUR), Bipolar (PGC 2025), MDD (PGC 2025), Neuroticism (Nagel 2018, EUR), Schizophrenia (Trubetskoy 2022, EUR) |

See [config/ad_config.yaml](config/ad_config.yaml) for exact file paths and the
notes on dropped traits (FastingInsulin, FastingGlucose, Chen HbA1c).
