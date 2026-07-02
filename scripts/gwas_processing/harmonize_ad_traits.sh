#!/bin/bash
#BSUB -J harmonize_ad_traits
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 8
#BSUB -R "rusage[mem=8000] span[hosts=1]"
#BSUB -W 8:00
#BSUB -o logs/lsf/harmonize_ad_traits_%J.out
#BSUB -e logs/lsf/harmonize_ad_traits_%J.err

# Harmonize the 7 AD-specific trait GWAS into the bNMF pipeline format.
# Chains: (optional preprocess) -> check_build.py -> harmonize_sumstats.py.
# Output -> <bnmf_ad>/sumstats/harmonized/*.processed.txt.gz
#
# REQUIRES the gwaslab conda env (activated below). To (re)build it:
#   ml anaconda3/2025.06
#   conda create -y -n gwaslab python=3.10 --solver=classic
#   source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate gwaslab
#   pip install gwaslab pandas numpy scipy
#   # gwaslab 4.1.9 import fix (invalid pd.NA type annotation):
#   sed -i 's/Union\[str, pd\.NA\]/Union[str, None]/g' \
#     "$(python -c 'import gwaslab,os;print(os.path.dirname(gwaslab.__file__))')/hm/hm_harmonize_sumstats.py"

set -uo pipefail
export OPENBLAS_NUM_THREADS=1

# --- Activate gwaslab env (self-contained for the LSF job) ---------------
ANACONDA_MODULE="${ANACONDA_MODULE:-anaconda3/2025.06}"
CONDA_ENV="${CONDA_ENV:-gwaslab}"
module load "${ANACONDA_MODULE}" 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"
python -c "import gwaslab" 2>/dev/null || { echo "ERROR: gwaslab not importable in env '${CONDA_ENV}'"; exit 1; }

# --- Paths ---------------------------------------------------------------
# PROJECT_ROOT = the bnmf_ad clone (auto-detected from this script's location).
PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${PROJECT_ROOT}"

SCRIPTS="${PROJECT_ROOT}/scripts/gwas_processing"
# Raw AD trait GWAS (override with RAW_DIR=... if you move them).
RAW_DIR="${RAW_DIR:-/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/gwas_for_bnmf}"
HARMONIZED="${PROJECT_ROOT}/sumstats/harmonized"
STAGE_DIR="${PROJECT_ROOT}/sumstats/staged"
BUILD_DIR="${PROJECT_ROOT}/sumstats/build_checks"
QC_DIR="${PROJECT_ROOT}/sumstats/qc"
PREFERRED_BUILD="19"
N_CORES=8

mkdir -p "${HARMONIZED}" "${STAGE_DIR}" "${BUILD_DIR}" "${QC_DIR}" logs/lsf

# --- File table: RAW | OUTPUT_BASE | KIND | SRC_BUILD -------------------
#   KIND:      standard | neuroticism | pgc_vcf  (preprocessing needed?)
#   SRC_BUILD: 19 or 38  (used ONLY to override a failed build detection;
#              harmonize still lifts to PREFERRED_BUILD=19)
FILES=(
  "GCST90027163_buildGRCh37.tsv.gz|vanRheenen_NatGenet_2021.ALS.META.GRCh37|standard|19"
  "GCST90275127.tsv|Kim_NatGenet_2023.Parkinsons.META.GRCh37|standard|19"
  "GCST90726396.tsv|Timsina_NatComms_2026.CSF_Abeta42.EUR.GRCh37|standard|38"
  "daner_pgc4_bd_multi_no23andMe_neff75_dentrem_HRCfrq|PGC_BD_2025.Bipolar.META.GRCh37|standard|19"
  "pgc-mdd2025_no23andMe_div_v3-49-46-01.tsv|PGC_MDD_2025.MDD.META.GRCh37|pgc_vcf|19"
  "sumstats_neuroticism_ctg_format.txt.gz|Nagel_NatGenet_2018.Neuroticism.EUR.GRCh37|neuroticism|19"
  "PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv|Trubetskoy_Nature_2022.Schizophrenia.EUR.GRCh37|pgc_vcf|19"
)

TOTAL=${#FILES[@]}; COUNT=0; FAILED=0
for ENTRY in "${FILES[@]}"; do
  IFS='|' read -r RAW OUTBASE KIND SRCBUILD <<< "${ENTRY}"
  COUNT=$((COUNT + 1))

  RAW_FILE="${RAW_DIR}/${RAW}"
  OUT_FILE="${HARMONIZED}/${OUTBASE}.processed.txt.gz"
  BUILD_JSON="${BUILD_DIR}/${OUTBASE}.build.json"
  QC_JSON="${QC_DIR}/${OUTBASE}.qc.json"

  if [ -f "${OUT_FILE}" ]; then echo "[${COUNT}/${TOTAL}] SKIP (exists): ${OUTBASE}"; continue; fi
  if [ ! -f "${RAW_FILE}" ]; then echo "[${COUNT}/${TOTAL}] SKIP (not found): ${RAW}"; FAILED=$((FAILED+1)); continue; fi

  echo ""; echo "[${COUNT}/${TOTAL}] ${OUTBASE}  (kind=${KIND}, src_build=${SRCBUILD})"

  # 1. Preprocess awkward formats into a staged clean file
  case "${KIND}" in
    neuroticism)
      INFILE="${STAGE_DIR}/${OUTBASE}.staged.txt.gz"
      python "${SCRIPTS}/prep_ad_trait_gwas.py" --input "${RAW_FILE}" --output "${INFILE}" --kind neuroticism || { echo "  FAILED (prep)"; FAILED=$((FAILED+1)); continue; } ;;
    pgc_vcf)
      INFILE="${STAGE_DIR}/${OUTBASE}.staged.tsv"
      python "${SCRIPTS}/prep_ad_trait_gwas.py" --input "${RAW_FILE}" --output "${INFILE}" --kind pgc_vcf || { echo "  FAILED (prep)"; FAILED=$((FAILED+1)); continue; } ;;
    *)
      INFILE="${RAW_FILE}" ;;
  esac

  # 2. Build check
  if [ ! -f "${BUILD_JSON}" ]; then
    echo "  Build check..."
    python "${SCRIPTS}/check_build.py" --input-file "${INFILE}" --output-file "${BUILD_JSON}"
  fi
  # Override only a FAILED detection, using the known source build
  if grep -q '"actual_build": "Unknown"' "${BUILD_JSON}" 2>/dev/null; then
    echo "  Build detection failed -> overriding to GRCh source build ${SRCBUILD}"
    python -c "
import json
d=json.load(open('${BUILD_JSON}'));
d['actual_build']='${SRCBUILD}'; d['override_reason']='detection failed; source build from file table'
json.dump(d, open('${BUILD_JSON}','w'), indent=2)"
  fi

  # 3. Harmonize (no --ref-gwas: trait alignment happens inside the bNMF prep step)
  echo "  Harmonizing -> ${OUT_FILE}"
  if python "${SCRIPTS}/harmonize_sumstats.py" \
       --input-file "${INFILE}" \
       --build-info "${BUILD_JSON}" \
       --output-file "${OUT_FILE}" \
       --qc-json "${QC_JSON}" \
       --preferred-build "${PREFERRED_BUILD}" \
       --n-cores "${N_CORES}"; then
    echo "  Done: ${OUTBASE}"
  else
    echo "  FAILED: ${OUTBASE}"; FAILED=$((FAILED+1))
  fi
done

echo ""; echo "=== AD trait harmonization complete ==="
echo "  Total: ${TOTAL}, Failed: ${FAILED}"
echo "  Harmonized files in ${HARMONIZED}:"
ls -1 "${HARMONIZED}"/*.META.GRCh37.processed.txt.gz "${HARMONIZED}"/*.EUR.GRCh37.processed.txt.gz 2>/dev/null | sed 's/^/    /'
