# AD bNMF Pipeline — Setup & Plan

Adapting the cardiometabolic bNMF pipeline (from `latlio/multiancestry_polygenic`)
to discover **genetic subtypes of Alzheimer's disease (AD)**.

Repo: https://github.com/accote45/bnmf_ad

## Concept

The pipeline builds a **variant × trait Z-score matrix** and factorizes it with bNMF:

- **Rows = variant universe.** Genome-wide-significant, MAF-filtered, LD-clumped SNPs
  from the *reference* GWAS. In the original: CAD + T2D. **Here: the AD GWAS.**
- **Columns = trait panel.** A set of quantitative-trait GWAS; each row variant's
  Z-score is pulled from each trait GWAS.
- **bNMF → K clusters** = latent genetic subtypes, each defined by a trait signature
  (W = variants × K, H = K × traits).

## Environment

- **Runs on Minerva (Mount Sinai LSF).** `module load R/4.2.0`, plink, and 1000G
  reference panels are expected there.
- **Code is edited on the Mac**, pushed to GitHub, pulled on Minerva, run there.
- **Scope (v1): EUR only** — one variant universe, one EUR LD panel.

## Workflow

```
edit on Mac  ->  git push  ->  git pull on Minerva  ->  run Rscript / bsub
```

Data (AD GWAS, trait GWAS, 1000G panels) lives on Minerva under
`sumstats/`, `reference/` — both git-ignored, never committed.

## Pipeline entry points (what we keep)

| Step | Script | Role |
|------|--------|------|
| Harmonize GWAS | `scripts/gwas_processing/harmonize_sumstats.py` (+ `check_build.py`, `column_map.py`) | Raw sumstats -> `*.processed.txt.gz` (VAR_ID, aligned alleles, Z/BETA/SE, GRCh37) |
| Run bNMF | `scripts/a1_analysis/01_run_bnmf.R` | QC + clump + Z-matrix + proxy + allele-align + bNMF |
| bNMF math | `scripts/a1_analysis/bnmf_algorithm.R` | Unchanged engine |
| Prep helpers | `scripts/a1_analysis/prep_bnmf.R` | QC / clumping / Z-matrix / allele alignment |
| Visualize | `scripts/a1_analysis/02_visualize_results.R` | W/H heatmaps, trait-signature figures |

Single ancestry => **Snakemake not required**; run the Rscript directly.

Dropped (CAD/T2D-specific): `09_cad_t2d_spectrum.R`, CAD/T2D cluster-label blocks,
the multi-ancestry comparison scripts.

## Key config changes (new `config/ad_config.yaml`)

- `ancestries: [EUR]`
- `ref_gwas.EUR: { AD: sumstats/harmonized/<AD_gwas>.processed.txt.gz }`
- `trait_gwas.EUR: { ...AD-relevant panel... }`
- `bnmf.p_threshold.EUR: 5e-8` (loosen if AD GWAS is underpowered)
- `allele_alignment.conflict_rule: first_occurrence`
  (single reference => orient each variant to the AD risk allele; no two-disease
  conflict to resolve — the CAD/T2D `strongest` logic isn't needed)
- EUR 1000G panel + HapMap3 paths as in the original `a1_config.yaml`

## Trait panel (draft — refine)

**Reuse from existing cardiometabolic harmonized set (if readable on Minerva):**
- Lipids / ApoE axis: LDL, HDL, TotalCholesterol, Triglycerides, ApoB, ApoA
- Immune / inflammatory: CRP, Lymphocyte, Monocyte, Neutrophil, WBC counts
- Metabolic: BMI, FastingGlucose, Hba1c (+ T2D as a trait)
- Vascular: SBP, DBP

**New AD-specific GWAS to source + harmonize:**
- Fluid biomarkers: CSF/plasma Abeta42, p-tau, total-tau, pTau181, GFAP, NfL
- Neuroimaging: hippocampal volume, total brain volume, WMH, cortical thickness
- Cognitive / reserve: general cognitive function, educational attainment
- Related: PD / ALS / FTD-LBD, parental lifespan / longevity, depression, sleep

## Reuse source (confirmed)

- Harmonized cardiometabolic sumstats: `SRC = /sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic/sumstats/harmonized/`
- **Read access CONFIRMED** (group `paul_oreilly`; 569 files). Format:
  `VAR_ID  RSID  Effect_Allele  P_VALUE  BETA  SE  N  MAF  EAF`, VAR_ID = `chr_pos_A1_A2`
  with alleles alphabetically sorted.
- Files are ~1 GB each -> **symlink, never copy**. `scripts/setup/link_reuse_sumstats.sh`
  links the 19 v1 reuse files into `sumstats/harmonized/`.

## AD GWAS reference (confirmed)

- Source = the **127 significant loci** from the full-cohort AD analysis (not a subset).
- Build: **hg19/GRCh37** (no liftover). Lead SNP coded `chr:pos:effect:other`.
- Rich file: rsID, alleles (ref/other/effect), beta, SE, OR, z_value, n_case/n_control,
  gene annotations (Predicted effector gene, ClosestGene, FLAMES).
- Processing: small dedicated converter (127 rows) -> harmonized reference file
  (`VAR_ID` alphabetical, `Effect_Allele`, `BETA`/Z, `P`). NOT the genome-wide harmonizer.

## Open items / TODO

- [x] Confirm read access to `lioul01`'s harmonized sumstats — DONE.
- [x] AD GWAS form/build confirmed (127 loci, hg19).
- [ ] Decide data home on Minerva (personal `cotea02`/`accote45` space vs shared lab).
- [ ] Clone `bnmf_ad` on Minerva; run `scripts/setup/link_reuse_sumstats.sh`.
- [ ] **Source the new AD-specific trait GWAS** (in progress — takes time).
- [ ] Draft the 127-loci converter (ready to write; needs the actual file to test).
- [ ] Finalize v1 trait list once new GWAS are in hand.

## Phases

0. Repo/workflow setup  ← (this doc; .gitignore cleanup done)
1. Harmonize AD GWAS (reference / variant universe)
2. Assemble trait panel (reuse + new)
3. Write `config/ad_config.yaml`
4. Run `01_run_bnmf.R` + `02_visualize_results.R` (EUR)
5. Iterate: tune K, p-threshold, correlation filter; interpret AD subtypes
