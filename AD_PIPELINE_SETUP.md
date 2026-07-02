# AD bNMF Pipeline — Setup & Plan

Adapting the cardiometabolic bNMF pipeline (from `latlio/multiancestry_polygenic`)
to discover **genetic subtypes of Alzheimer's disease (AD)**.

Repo: https://github.com/accote45/bnmf_ad

## ⏩ RESUME HERE (status as of 2026-07-02)

- **Minerva project root:** `/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/results/bnmf_ad`
  (raw AD trait GWAS at `.../pathway_prs_ad/data/gwas_for_bnmf/`).
- **Done:** repo cloned on Minerva; `reference/1kg_eur` + 19 reuse `.META.` traits symlinked;
  AD reference converted (127 loci, `AD_PGC3_top127.META.GRCh37`); gwaslab env built
  (`/hpc/users/cotea02/.conda/envs/gwaslab`, gwaslab 4.1.9 + pd.NA import patch);
  6 of 7 new trait GWAS harmonized (ALS, Parkinsons, CSF Aβ42, Bipolar, MDD, Schizophrenia).
- **In progress:** neuroticism harmonization (last of 7) — rerun after the `dtype=str` prep fix.
  Check: `ls -l sumstats/harmonized/*Neuroticism*.processed.txt.gz`.
- **NEXT STEP once neuroticism lands →** launch the bNMF:
  ```
  module load R/4.2.0 && module load plink/1.90b6.21
  Rscript scripts/a1_analysis/01_run_bnmf.R --config config/ad_config.yaml --ancestry META
  ```
- **Machine switch:** all code is in GitHub, all data + jobs on Minerva. On a new computer:
  `git clone https://github.com/accote45/bnmf_ad.git`, then SSH to Minerva as usual.

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
- **gwaslab conda env** (for harmonization): `gwaslab` env, python 3.10, gwaslab 4.1.9.
  Built with `ml anaconda3/2025.06` + classic-solver `conda create` + pip. Needed a
  one-line fix for a gwaslab 4.1.9 import bug (`Union[str, pd.NA]` invalid annotation) —
  see the recipe in `scripts/gwas_processing/harmonize_ad_traits.sh` header. The batch
  script activates this env itself.
- **Code is edited on the Mac**, pushed to GitHub, pulled on Minerva, run there.
- **Scope (v1): trans-ancestry / META** — the 127 AD loci and the trait GWAS are
  both multi-ancestry, so the panel uses the `.META.` version of every trait.
  LD reference = EUR 1000G (standard approximation for a trans-ancestry set, as in
  the original pipeline's META arm). Run label = `META`.

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

**Reuse from existing harmonized set — 19 traits, `.META.` versions** (confirmed readable):
- Lipids / ApoE axis: HDL, LDL, TotalCholesterol, Triglycerides (Graham2021); ApoA, ApoB (Karczewski2025)
- Immune / inflammatory: CRP, Lymphocyte, Monocyte, Neutrophil, WBC, Eosinophil (Karczewski2025)
- Metabolic: T2D (Suzuki2024), FastingGlucose (Downie2022), RandomGlucose (Lagou2023), Hba1c (Chen2021), BMI (Karczewski2025)
- Vascular: SBP, DBP (Karczewski2025)
- Note: FastingInsulin dropped (no META version); RandomGlucose used instead. New AD
  traits must also be `.META.` for column consistency.

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

## 1000G EUR reference panel (confirmed)

- Location: `/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic/reference/1kg_eur/`
- **Read access CONFIRMED**; per-chromosome plink1 bfiles `1000G.EUR.QC.{1..22}.{bed,bim,fam}` (GRCh37).
- Wire in via directory symlink (no copy):
  ```
  REF=/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic/reference/1kg_eur
  mkdir -p reference && ln -sfn "$REF" reference/1kg_eur
  ```
- Config `ld_clump.ref_panel.EUR: reference/1kg_eur/1000G.EUR.QC` resolves through it.

## AD GWAS reference (confirmed)

- File on Minerva: `/sc/arion/projects/paul_oreilly/data/GWASs/NonBiobanks/raw_data/ad/PGC3_Unpublished/uffleman_top_loci_ad_gwas.xlsx`
- Convert with `scripts/gwas_processing/convert_ad_loci.R` (default `--input` points at the file above)
  -> `sumstats/harmonized/AD_PGC3_top127.META.GRCh37.processed.txt.gz`.

- Source = the **127 significant loci** from the full-cohort AD analysis (not a subset).
- Build: **hg19/GRCh37** (no liftover). Lead SNP coded `chr:pos:effect:other`.
- Rich file: rsID, alleles (ref/other/effect), beta, SE, OR, z_value, n_case/n_control,
  gene annotations (Predicted effector gene, ClosestGene, FLAMES).
- Processing: small dedicated converter (127 rows) -> harmonized reference file
  (`VAR_ID` alphabetical, `Effect_Allele`, `BETA`/Z, `P`). NOT the genome-wide harmonizer.

## Open items / TODO

- [x] Confirm read access to `lioul01`'s harmonized sumstats — DONE.
- [x] AD GWAS form/build confirmed (127 loci, hg19).
- [x] 1000G EUR reference panel access confirmed.
- [x] Draft the 127-loci converter (`convert_ad_loci.R`) + EUR config skeleton (`ad_config.yaml`).
- [ ] Decide data home on Minerva (personal `cotea02` space vs shared lab).
- [ ] Clone `bnmf_ad` on Minerva; symlink `reference/1kg_eur`; run `link_reuse_sumstats.sh`.
- [ ] Run `convert_ad_loci.R`; paste console output to confirm column parsing.
- [x] Source + review the 7 new AD-specific trait GWAS.
- [x] Harmonization tooling for them (`harmonize_ad_traits.sh` + `prep_ad_trait_gwas.py`
      + column_map patches); appended to `config/ad_config.yaml` (26 traits total).
- [ ] Run `harmonize_ad_traits.sh` on Minerva; verify the 7 QC JSONs + first-run output.
- [ ] Symlink reuse (.META.) files; run converter; run bNMF.

## New AD trait GWAS (raw at .../lab/cotea02/pathway_prs_ad/data/gwas_for_bnmf/)

| Config key | Raw file | Build | Notes |
|-----------|----------|-------|-------|
| ALS_vanRheenen2021 | GCST90027163_buildGRCh37.tsv.gz | 37 | cross-anc; turnkey |
| Parkinsons_Kim2023 | GCST90275127.tsv | 37 | multi-anc MR-MEGA; turnkey |
| CSF_Abeta42_Timsina2026 | GCST90726396.tsv | **38** | EUR; auto liftover 38->37 |
| Bipolar_PGC2025 | daner_pgc4_bd_multi_...HRCfrq | 37 | daner, space-delim, OR->beta |
| MDD_PGC2025 | pgc-mdd2025_...tsv | 37 | PGC-VCF: prep strips `##` |
| Neuroticism_Nagel2018 | sumstats_neuroticism_ctg_format.txt.gz | 37 | EUR; prep derives BETA/SE from Z |
| Schizophrenia_Trubetskoy2022 | PGC3_SCZ_wave3.primary...vcf.tsv | 37 | EUR; PGC-VCF: prep strips `##` |

Harmonize: `bsub < scripts/gwas_processing/harmonize_ad_traits.sh` (needs a python env
with gwaslab/pandas/numpy/scipy). Outputs -> `sumstats/harmonized/`.

## Minerva bring-up (once new GWAS are ready)

```
git clone https://github.com/accote45/bnmf_ad.git && cd bnmf_ad
REF=/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic/reference/1kg_eur
mkdir -p reference && ln -sfn "$REF" reference/1kg_eur
bash scripts/setup/link_reuse_sumstats.sh
module load R/4.2.0 && module load plink/1.90b6.21
Rscript scripts/gwas_processing/convert_ad_loci.R
Rscript scripts/a1_analysis/01_run_bnmf.R --config config/ad_config.yaml --ancestry META
```

## Phases

0. Repo/workflow setup  ← (this doc; .gitignore cleanup done)
1. Harmonize AD GWAS (reference / variant universe)
2. Assemble trait panel (reuse + new)
3. Write `config/ad_config.yaml`
4. Run `01_run_bnmf.R` + `02_visualize_results.R` (EUR)
5. Iterate: tune K, p-threshold, correlation filter; interpret AD subtypes
