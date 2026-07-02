#!/bin/bash
set -euo pipefail

SUMSTATS="/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic/sumstats"
mkdir -p "${SUMSTATS}"

# Array of "URL|TARGET_FILENAME" pairs
FILES=(
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90692001-GCST90693000/GCST90692778/harmonised/GCST90692778.h.tsv.gz|Karczewski_NatGenet_2025.Creatinine.META.GRCh37.txt.gz"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90691001-GCST90692000/GCST90691572/harmonised/GCST90691572.h.tsv.gz|Karczewski_NatGenet_2025.Creatinine.EUR.GRCh37.txt.gz"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90692001-GCST90693000/GCST90692116/harmonised/GCST90692116.h.tsv.gz|Karczewski_NatGenet_2025.Creatinine.SAS.GRCh37.txt.gz"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90301001-GCST90302000/GCST90301952/harmonised/GCST90301952.h.tsv.gz|Karjalainen_Nature_2024.Creatinine.META.GRCh37.tsv.gz"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90278001-GCST90279000/GCST90278624/harmonised/GCST90278624.h.tsv.gz|Chen_CellGenomics_2023.Creatinine.EAS.GRCh37.tsv.gz"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90018001-GCST90019000/GCST90018759/harmonised/34594039-GCST90018759-EFO_0004518-Build37.f.tsv.gz|Sakaue_NatGenet_2021.Creatinine.EAS.GRCh37.txt.gz"
)

SUCCESS=0
SKIPPED=0
FAILED=0

for ENTRY in "${FILES[@]}"; do
    IFS='|' read -r URL FNAME <<< "${ENTRY}"
    OUTPATH="${SUMSTATS}/${FNAME}"
    if [ -f "${OUTPATH}" ]; then
        echo "SKIP (exists): ${FNAME}"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi
    echo "Downloading: ${FNAME}"
    if wget -q -O "${OUTPATH}" "${URL}"; then
        echo "  OK: ${FNAME}"
        SUCCESS=$((SUCCESS + 1))
    else
        echo "  FAILED: ${FNAME}"
        rm -f "${OUTPATH}"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "Done. Success: ${SUCCESS}, Skipped: ${SKIPPED}, Failed: ${FAILED}"
