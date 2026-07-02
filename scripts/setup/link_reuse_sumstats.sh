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

# One canonical .META. file per v1 reuse trait (trans-ancestry framing;
# names taken from config/a1_config.yaml's META block).
FILES=(
  # --- Lipid / APOE axis ---
  "Graham_Nature_2021.HDL.META.GRCh37.processed.txt.gz"
  "Graham_Nature_2021.LDL.META.GRCh37.processed.txt.gz"
  "Graham_Nature_2021.TotalCholesterol.META.GRCh37.processed.txt.gz"
  "Graham_Nature_2021.Triglycerides.META.GRCh37.processed.txt.gz"
  "Karczewski_NatureGenetics_2025.ApoA.META.GRCh37.processed.txt.gz"
  "Karczewski_NatureGenetics_2025.ApoB.META.GRCh37.processed.txt.gz"
  # --- Immune / inflammatory ---
  "Karczewski_NatGenet_2025.CRP.META.GRCh37.processed.txt.gz"
  "Karczewski_NatGenet_2025.LymphocyteCount.META.GRCh37.processed.txt.gz"
  "Karczewski_NatGenet_2025.MonocyteCount.META.GRCh37.processed.txt.gz"
  "Karczewski_NatGenet_2025.NeutrophilCount.META.GRCh37.processed.txt.gz"
  "Karczewski_NatGenet_2025.WBC.META.GRCh37.processed.txt.gz"
  "Karczewski_NatGenet_2025.Eosinophilcount.META.GRCh37.processed.txt.gz"
  # --- Metabolic / glycemic ---
  "Suzuki_Nature_2024.t2d.META.GRCh37.processed.txt.gz"
  "Downie_Diabetologia_2022.FastingGlucose.META.GRCh37.processed.txt.gz"
  "Lagou_NatGenet_2023.RandomGlucose.META.GRCh37.processed.txt.gz"
  "Chen_NatureGenetics_2021.Hba1c.META.GRCh37.processed.txt.gz"
  "Karczewski_NatGenet_2025.BMI.META.GRCh37.processed.txt.gz"
  # --- Vascular ---
  "Karczewski_NatGenet_2025.SBP.META.GRCh37.processed.txt.gz"
  "Karczewski_NatGenet_2025.DBP.META.GRCh37.processed.txt.gz"
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
