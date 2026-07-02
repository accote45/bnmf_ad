#!/bin/bash
# ------------------------------------------------------------------
# submit_harmonize_splitukb.sh
# ------------------------------------------------------------------
# Submits 74 LSF jobs (1 per trait) to convert PLINK GWAS files
# to harmonized format for the split-UKB bNMF comparison.
#
# Queue: premium, Walltime: 4:00, Memory: 16GB
# ------------------------------------------------------------------

set -euo pipefail

PROJECT_DIR="/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
SOURCE_DIR="/sc/arion/projects/psychgen/projects/prs/extreme_traits/July_EUR_noCancerCase_residFirst/gwas/inverse/all"
OUTPUT_DIR="${PROJECT_DIR}/sumstats/harmonized_splitukb"
SCRIPT="${PROJECT_DIR}/scripts/a1_analysis/convert_plink_to_harmonized.R"
LOG_DIR="${PROJECT_DIR}/logs/harmonize_splitukb"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

# Trait map: field_id trait_name
# Derived from config/splitukb_trait_map.yaml
declare -a TRAITS=(
  "3761 AgeHayFever"
  "30600 Albumin"
  "20157 AlphanumTrail"
  "3143 AnkleSpacing"
  "30630 ApolipoproteinA"
  "21001 BMI"
  "23105 BMR"
  "20022 BirthWeight"
  "23099 BodyFatPct"
  "23106 BodyImpedance"
  "30680 Calcium"
  "30690 Cholesterol"
  "30700 Creatinine"
  "30210 EosinophilPct"
  "30150 Eosinophils"
  "20150 FEVBest"
  "20256 FEVZscore"
  "23129 FFTrunkMass"
  "20151 FVCBest"
  "20016 FluidIntellect"
  "2217 GlassesAge"
  "30740 Glucose"
  "30750 GlycatedHaem"
  "46 GripStrength"
  "3536 HRTAge"
  "30030 HaematocritPct"
  "30020 HaemoglobinConc"
  "3148 HeelBMD"
  "3144 HeelBUA"
  "50 Height"
  "30770 IGF1"
  "30280 ImmRTFrac"
  "23115 LegFatPct"
  "23107 LegImpedance"
  "30180 LymphocytePct"
  "30120 Lymphocytes"
  "20023 MatchIDTime"
  "30060 MeanCPHeme"
  "30040 MeanCPVol"
  "30260 MeanRTVol"
  "2714 MenarcheAge"
  "3581 MenopauseAge"
  "30190 MonocytePct"
  "30130 Monocytes"
  "20127 Neuroticism"
  "30200 NeutrophilPct"
  "30140 Neutrophils"
  "20156 NumericTrail"
  "1050 OutdoorTime"
  "4199 PWNotch"
  "30810 Phosphate"
  "30090 PlateletCrit"
  "20154 PredFEVPct"
  "4194 PulseRate"
  "30070 RBCDistWidth"
  "30010 RBCs"
  "30240 ReticulocytePct"
  "30250 Reticulocytes"
  "20159 SDMatches"
  "30830 SHBG"
  "20015 SittingHeight"
  "30530 SodiumInUrine"
  "30270 SpheredCellV"
  "1070 TVWatching"
  "30850 Testosterone"
  "4288 TimeToAnswer"
  "30860 TotalProtein"
  "30870 Triglycerides"
  "30670 Urea"
  "30890 VitaminD"
  "30000 WBCs"
  "48 WaistCircum"
  "21002 Weight"
)

echo "Submitting ${#TRAITS[@]} harmonization jobs..."

for entry in "${TRAITS[@]}"; do
  field_id=$(echo "${entry}" | awk '{print $1}')
  trait_name=$(echo "${entry}" | awk '{print $2}')

  bsub -P acc_paul_oreilly \
       -q premium \
       -W 4:00 \
       -n 1 \
       -R "rusage[mem=16000]" \
       -J "harm_${field_id}" \
       -o "${LOG_DIR}/harm_${field_id}.stdout" \
       -e "${LOG_DIR}/harm_${field_id}.stderr" \
       "module load R/4.1.0 && Rscript ${SCRIPT} \
         --source-dir ${SOURCE_DIR} \
         --output-dir ${OUTPUT_DIR} \
         --field-id ${field_id} \
         --trait-name ${trait_name}"

  echo "  Submitted: ${field_id} (${trait_name})"
done

echo "All ${#TRAITS[@]} jobs submitted."
