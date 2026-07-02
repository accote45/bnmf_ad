#!/usr/bin/env bash
# compute_prs_aou.sh
# Terminal PRS computation for the All of Us replication of B1.
# Mirrors 01_compute_cluster_prs.R exactly, but scores the AoU ACAF-threshold
# PLINK genotypes (GRCh38) instead of UKB. Matches the GRCh37-derived cluster
# weights to AoU variants by lifted CHR38:POS38 + alleles (primary) with rsID as
# a fallback, then runs `plink2 --score` per PLINK fileset and merges to one
# per-person cluster-PGS table.
#
# Run from the AoU terminal:
#   bash compute_prs_aou.sh
#
# Inputs : $WEIGHTS  (b1_cluster_weights.tsv, uploaded to the workspace)
#          $GENO_DIR (AoU ACAF-threshold PLINK .bed/.bim/.fam)
# Output : $OUTDIR/cluster_prs_all.tsv  (FID, IID=person_id, PRS_K1..PRS_K10)
set -euo pipefail

# ============================ PARAMETERS ====================================
# Edit these to match your workspace, then run. All echoed below before work.
GENO_DIR="${GENO_DIR:-workspaces/duplicateofstatinresponse/genotype/plink}"  # ACAF PLINK dir
WEIGHTS="${WEIGHTS:-b1_cluster_weights.tsv}"          # uploaded weights file
OUTDIR="${OUTDIR:-cluster_prs_out}"                   # local output dir
PLINK2="${PLINK2:-plink2}"                            # plink2 on PATH in AoU
KEEP="${KEEP:-}"                                      # optional --keep file (FID IID); empty = all
N_CLUSTERS="${N_CLUSTERS:-10}"                        # K1..K10
UPLOAD="${UPLOAD:-1}"                                 # 1 = gsutil cp result to $WORKSPACE_BUCKET

mkdir -p "$OUTDIR"
SCORE_COL_NUMS="3-$((2 + N_CLUSTERS))"               # score file: col1=ID col2=A1 col3..=Wbeta

echo "=== compute_prs_aou: parameters ==="
echo "  GENO_DIR : $GENO_DIR"
echo "  WEIGHTS  : $WEIGHTS"
echo "  OUTDIR   : $OUTDIR"
echo "  PLINK2   : $PLINK2"
echo "  KEEP     : ${KEEP:-<none, score all samples>}"
echo "  N_CLUSTERS / score-col-nums : $N_CLUSTERS / $SCORE_COL_NUMS"
echo "  UPLOAD to bucket : $UPLOAD"
echo

[ -f "$WEIGHTS" ] || { echo "ERROR: weights file not found: $WEIGHTS" >&2; exit 1; }

# Discover all PLINK filesets in GENO_DIR (handles single genome-wide or per-chr)
mapfile -t BEDS < <(ls "$GENO_DIR"/*.bed 2>/dev/null | sort)
[ "${#BEDS[@]}" -gt 0 ] || { echo "ERROR: no .bed files under $GENO_DIR" >&2; exit 1; }
echo "Found ${#BEDS[@]} PLINK fileset(s):"
printf '  %s\n' "${BEDS[@]}"
echo

# ===================== STEP 1: build per-fileset score files =================
# An embedded R step matches the weights to each fileset's .bim and writes a
# plink2 score file (ID = the fileset's own variant ID; A1 = effect allele,
# strand-aligned). Primary join: CHR38:POS38 + allele-set match. Fallback: rsID.
echo "--- Step 1: matching weights to genotype .bim files ---"
MANIFEST="$OUTDIR/score_manifest.txt"   # lines: <bed_prefix>\t<score_file>
: > "$MANIFEST"

Rscript - "$WEIGHTS" "$OUTDIR" "$MANIFEST" "${BEDS[@]}" <<'RSCRIPT'
suppressMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
weights_path <- a[1]; outdir <- a[2]; manifest <- a[3]
beds <- a[-(1:3)]

w <- fread(weights_path)
wbeta_cols <- grep("^Wbeta_K", colnames(w), value = TRUE)
w[, CHR38 := as.character(CHR38)]
w[, CHR38 := sub("^chr", "", CHR38)]

compl <- c(A = "T", T = "A", C = "G", G = "C")
is_palindrome <- function(a1, a2) (a1 == compl[a2])  # A/T or C/G — strand-ambiguous

total_matched <- 0L
for (bed in beds) {
  prefix  <- sub("\\.bed$", "", bed)
  bimfile <- paste0(prefix, ".bim")
  if (!file.exists(bimfile)) { cat(sprintf("  SKIP %s: no .bim\n", prefix)); next }

  bim <- fread(bimfile, header = FALSE,
               col.names = c("bCHR", "bID", "bCM", "bPOS", "bA1", "bA2"))
  bim[, bCHR := sub("^chr", "", as.character(bCHR))]

  ## --- primary: position match, then allele-set agreement ---
  m <- merge(w, bim, by.x = c("CHR38", "POS38"), by.y = c("bCHR", "bPOS"),
             allow.cartesian = TRUE)
  if (nrow(m) > 0) {
    m[, exact := (REF == bA1 & ALT == bA2) | (REF == bA2 & ALT == bA1)]
    # strand flip only for non-palindromic variants where exact failed
    m[, flip := !exact & !is_palindrome(REF, ALT) &
        ((compl[REF] == bA1 & compl[ALT] == bA2) |
         (compl[REF] == bA2 & compl[ALT] == bA1))]
    m <- m[exact | flip]
    # effect allele on the genotype strand
    m[, A1 := fifelse(flip, compl[Effect_Allele], Effect_Allele)]
    setorder(m, VAR_ID37, -exact)             # prefer exact over flip
    m <- m[!duplicated(VAR_ID37)]
  }
  matched_ids <- if (nrow(m) > 0) m$VAR_ID37 else character(0)

  ## --- fallback: rsID match for weights not matched by position ---
  w_un <- w[!VAR_ID37 %in% matched_ids & RSID != "." & !is.na(RSID)]
  if (nrow(w_un) > 0) {
    r <- merge(w_un, bim[, .(bID, bA1, bA2)], by.x = "RSID", by.y = "bID")
    if (nrow(r) > 0) {
      r[, exact := (REF == bA1 & ALT == bA2) | (REF == bA2 & ALT == bA1)]
      r <- r[exact]
      r[, A1 := Effect_Allele]
      r[, bID := RSID]
      r <- r[!duplicated(VAR_ID37)]
      m <- if (nrow(m) > 0) rbindlist(list(m, r), use.names = TRUE, fill = TRUE) else r
    }
  }

  if (is.null(m) || nrow(m) == 0) { cat(sprintf("  %s: 0 matched\n", basename(prefix))); next }

  # score file: ID (genotype variant id), A1 (effect allele), Wbeta_K1..K10
  score <- m[, c("bID", "A1", wbeta_cols), with = FALSE]
  setnames(score, c("bID", "A1"), c("ID", "A1"))
  score <- score[!duplicated(ID)]            # plink also gets ignore-dup-ids
  sf <- file.path(outdir, paste0("score_", basename(prefix), ".tsv"))
  fwrite(score, sf, sep = "\t")
  cat(sprintf("%s\t%s\n", prefix, sf), file = manifest, append = TRUE)

  total_matched <- total_matched + nrow(score)
  cat(sprintf("  %s: %d variants matched\n", basename(prefix), nrow(score)))
}
cat(sprintf("Total matched across filesets: %d / %d weights\n",
            total_matched, nrow(w)))
RSCRIPT

[ -s "$MANIFEST" ] || { echo "ERROR: no variants matched any fileset" >&2; exit 1; }
echo

# ===================== STEP 2: plink2 --score per fileset ====================
echo "--- Step 2: plink2 --score ---"
KEEP_ARG=""
[ -n "$KEEP" ] && KEEP_ARG="--keep $KEEP"

SSCORES=()
while IFS=$'\t' read -r prefix sf; do
  [ -n "$prefix" ] || continue
  out="$OUTDIR/$(basename "$prefix")"
  echo "  scoring $(basename "$prefix") ..."
  OMP_NUM_THREADS=1 "$PLINK2" \
    --bfile "$prefix" \
    $KEEP_ARG \
    --score "$sf" 1 2 header cols=+scoresums ignore-dup-ids \
    --score-col-nums "$SCORE_COL_NUMS" \
    --threads 1 \
    --out "$out" 2>&1 | grep -iE "error|warning|variants? loaded|valid predictors|samples" || true
  [ -f "$out.sscore" ] && SSCORES+=("$out.sscore") \
    || echo "    WARNING: $out.sscore not produced"
done < "$MANIFEST"

[ "${#SSCORES[@]}" -gt 0 ] || { echo "ERROR: no .sscore files produced" >&2; exit 1; }
echo

# ===================== STEP 3: merge sscores (sum across filesets) ===========
# Per script 01: plink2 names the K score-sum columns SCORE1_SUM..SCORE{N}_SUM in
# cluster order; sum them across filesets per person and rename to PRS_K1..K{N}.
echo "--- Step 3: merge per-fileset scores ---"
OUT_TSV="$OUTDIR/cluster_prs_all.tsv"

Rscript - "$N_CLUSTERS" "$OUT_TSV" "${SSCORES[@]}" <<'RSCRIPT'
suppressMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
n_clusters <- as.integer(a[1]); out_tsv <- a[2]; files <- a[-(1:2)]

merged <- NULL
for (f in files) {
  dt <- fread(f)
  setnames(dt, "#FID", "FID", skip_absent = TRUE)
  sum_cols <- grep("^SCORE.*_SUM$", colnames(dt), value = TRUE)   # file order = K1..KN
  chr_dt <- dt[, c("FID", "IID", sum_cols), with = FALSE]
  if (is.null(merged)) {
    merged <- chr_dt
  } else {
    merged <- merge(merged, chr_dt, by = c("FID", "IID"), all = TRUE,
                    suffixes = c("", ".new"))
    for (sc in sum_cols) {
      nc <- paste0(sc, ".new")
      if (nc %in% colnames(merged)) {
        merged[[sc]] <- fifelse(is.na(merged[[sc]]), 0, merged[[sc]]) +
                        fifelse(is.na(merged[[nc]]), 0, merged[[nc]])
        merged[[nc]] <- NULL
      }
    }
  }
}
sum_cols <- grep("^SCORE.*_SUM$", colnames(merged), value = TRUE)
setnames(merged, sum_cols, paste0("PRS_K", seq_len(n_clusters)))
fwrite(merged, out_tsv, sep = "\t")
cat(sprintf("Wrote %s: %d individuals x %d cluster PGS\n",
            out_tsv, nrow(merged), n_clusters))
RSCRIPT
echo

# ===================== STEP 4: upload to workspace bucket ====================
if [ "$UPLOAD" = "1" ] && [ -n "${WORKSPACE_BUCKET:-}" ]; then
  echo "--- Uploading to $WORKSPACE_BUCKET/b1_aou/ ---"
  gsutil cp "$OUT_TSV" "$WORKSPACE_BUCKET/b1_aou/cluster_prs_all.tsv"
fi

echo "=== Done. PRS file: $OUT_TSV ==="
