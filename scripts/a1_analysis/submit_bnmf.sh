#!/bin/bash
#BSUB -J bnmf_ad_META
#BSUB -P acc_paul_oreilly
#BSUB -q express
#BSUB -n 1
#BSUB -R "rusage[mem=64000] span[hosts=1]"
#BSUB -W 12:00
#BSUB -o logs/lsf/bnmf_ad_META_%J.out
#BSUB -e logs/lsf/bnmf_ad_META_%J.err

# Launch the AD bNMF (single ancestry arm) on Minerva LSF.
# Wraps: Rscript scripts/a1_analysis/01_run_bnmf.R --config <cfg> --ancestry <ANC>
# Resources mirror config/ad_config.yaml `hpc:` (express, 64 GB, 1 core, 12h).
#
# Submit from the bnmf_ad repo root:
#     bsub < scripts/a1_analysis/submit_bnmf.sh
# Override defaults via env vars, e.g.:
#     ANCESTRY=META CONFIG=config/ad_config.yaml bsub < scripts/a1_analysis/submit_bnmf.sh

set -uo pipefail
export OPENBLAS_NUM_THREADS=1

# --- Paths ---------------------------------------------------------------
# `bsub < script` runs from stdin, so $0 can't locate the script. LSF preserves
# the submission directory as CWD, so default PROJECT_ROOT to it.
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
if [ ! -f "${PROJECT_ROOT}/scripts/a1_analysis/01_run_bnmf.R" ]; then
  echo "ERROR: PROJECT_ROOT='${PROJECT_ROOT}' is not the bnmf_ad repo root."
  echo "       Submit from the repo root, or pass PROJECT_ROOT=/path/to/bnmf_ad."
  exit 1
fi
cd "${PROJECT_ROOT}"

CONFIG="${CONFIG:-config/ad_config.yaml}"
ANCESTRY="${ANCESTRY:-META}"
SELECT_K="${SELECT_K:-}"          # optional: force a specific K (e.g. SELECT_K=4)
[ -f "${CONFIG}" ] || { echo "ERROR: config not found: ${CONFIG}"; exit 1; }
EXTRA_ARGS=""
[ -n "${SELECT_K}" ] && EXTRA_ARGS="--select-k ${SELECT_K}"

mkdir -p logs/lsf

# --- Environment ---------------------------------------------------------
# `module` is not initialized in a non-login LSF batch shell; source it first.
source /etc/profile.d/modules.sh 2>/dev/null || true
module load R/4.2.0
module load plink/1.90b6.21

echo "Host:        $(hostname)"
echo "PROJECT_ROOT: ${PROJECT_ROOT}"
echo "Config:      ${CONFIG}"
echo "Ancestry:    ${ANCESTRY}"
echo "Select K:    ${SELECT_K:-<modal>}"
echo "R:           $(which Rscript)"
echo "plink:       $(which plink)"
echo "Started:     $(date)"

# --- Run -----------------------------------------------------------------
Rscript scripts/a1_analysis/01_run_bnmf.R --config "${CONFIG}" --ancestry "${ANCESTRY}" ${EXTRA_ARGS}
status=$?

echo "Finished:    $(date)  (exit ${status})"
exit ${status}
