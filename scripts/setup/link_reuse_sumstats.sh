#!/bin/bash
# link_reuse_sumstats.sh
# Symlink the v1 "reuse" cardiometabolic trait GWAS (already harmonized by
# lioul01) into this project's sumstats/harmonized/ so the AD config's
# repo-relative paths resolve through them. Files are ~1GB each — symlink,
# never copy. New AD-specific harmonized files coexist in the same dir.
#
# Run on Minerva from the bnmf_ad project root:
#   bash scripts/setup/link_reuse_sumstats.sh
#
# Read access to SRC is required (confirmed: group paul_oreilly can read).

set -uo pipefail

SRC="/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic/sumstats/harmonized"
DEST="sumstats/harmonized"

mkdir -p "$DEST"

# One canonical EUR file per v1 reuse trait (names taken from config/a1_config.yaml).
FILES=(
  # --- Lipid / APOE axis ---
  "Graham_Nature_2021.HDL.EUR.GRCh37.processed.txt.gz"
  "Graham_Nature_2021.LDL.EUR.GRCh37.processed.txt.gz"
  "Graham_Nature_2021.TotalCholesterol.EUR.GRCh37.processed.txt.gz"
  "Graham_Nature_2021.Triglycerides.EUR.GRCh37.processed.txt.gz"
  "Karczewski_NatureGenetics_2025.ApoA.EUR.GRCh37.processed.txt.gz"
  "Karczewski_NatureGenetics_2025.ApoB.EUR.GRCh37.processed.txt.gz"
  # --- Immune / inflammatory ---
  "Said_NatureComms_2022.CRP.EUR.GRCh37.processed.txt.gz"
  "Chen_Cell_2020.LymphocyteCount.EUR.GRCh37.processed.txt.gz"
  "Chen_Cell_2020.MonocyteCount.EUR.GRCh37.processed.txt.gz"
  "Chen_Cell_2020.NeutrophilCount.EUR.GRCh37.processed.txt.gz"
  "Chen_Cell_2020.WBC.EUR.GRCh37.processed.txt.gz"
  "Karczewski_NatGenet_2025.Eosinophilcount.EUR.GRCh37.processed.txt.gz"
  # --- Metabolic / glycemic ---
  "Suzuki_Nature_2024.t2d.EUR.GRCh37.processed.txt.gz"
  "Chen_NatGenet_2021.FastingGlucose.EUR.GRCh37.processed.txt.gz"
  "Chen_NatGenet_2021.FastingInsulin.EUR.GRCh37.processed.txt.gz"
  "Chen_NatureGenetics_2021.Hba1c.EUR.GRCh37.processed.txt.gz"
  "Karczewski_NatGenet_2025.BMI.EUR.GRCh37.processed.txt.gz"
  # --- Vascular ---
  "Keaton_NatureGenetics_2024.SBP.EUR.GRCh37.processed.txt.gz"
  "Evangelou_NatGenet_2018.DBP.EUR.GRCh37.processed.txt.gz"
)

linked=0; missing=0
missing_files=()
for f in "${FILES[@]}"; do
  if [[ -r "$SRC/$f" ]]; then
    ln -sfn "$SRC/$f" "$DEST/$f"
    echo "  linked: $f"
    ((linked++))
  else
    echo "  MISSING/unreadable: $f"
    missing_files+=("$f")
    ((missing++))
  fi
done

echo ""
echo "Linked $linked / ${#FILES[@]} files into $DEST/ (missing: $missing)"
if (( missing > 0 )); then
  echo "Missing files (check exact name in $SRC, or plan re-harmonization):"
  printf '  %s\n' "${missing_files[@]}"
  echo "Tip: ls \"$SRC\" | grep -iE 'HDL|CRP|BMI'  # find the real filename"
fi
