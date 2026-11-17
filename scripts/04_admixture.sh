#!/bin/bash
# =============================================================================
# Script 04: ADMIXTURE — Ancestry Inference
# Author: Niharika Joshi
# Description: Runs ADMIXTURE at K=2 to K=7 with cross-validation to
#              identify the optimal number of ancestry components.
#              Results are used for visualisation in R.
# =============================================================================

set -euo pipefail

# --- Configuration ---
QC_DIR="data/qc"
ADM_DIR="results/admixture"
CHR=22
K_MIN=2
K_MAX=7
THREADS=4   # adjust to your CPU count

mkdir -p "$ADM_DIR"

# ADMIXTURE must be run from the directory containing the BED file
cd "$ADM_DIR"
BED_FILE="../../$QC_DIR/chr${CHR}_pruned.bed"

echo "============================================"
echo "  ADMIXTURE Analysis (K=$K_MIN to K=$K_MAX)"
echo "============================================"

# --- Step 1: Run ADMIXTURE for each K with cross-validation ---
echo "[1/2] Running ADMIXTURE with cross-validation..."
echo "      This may take 10–30 minutes depending on SNP count."
echo ""

CV_ERRORS=()

for K in $(seq $K_MIN $K_MAX); do
    echo "  Running K=$K..."
    admixture --cv=10 -j${THREADS} "$BED_FILE" $K | tee "admixture_K${K}.log"
    
    # Extract cross-validation error from log
    CV_ERR=$(grep "CV error" "admixture_K${K}.log" | awk '{print $4}')
    CV_ERRORS+=("K=$K: $CV_ERR")
    echo "    CV error (K=$K): $CV_ERR"
    echo ""
done

echo "  ✓ ADMIXTURE runs complete"

# --- Step 2: Print cross-validation summary ---
echo ""
echo "[2/2] Cross-validation error summary:"
echo "      (Lower CV error = better model fit)"
echo ""
for err in "${CV_ERRORS[@]}"; do
    echo "  $err"
done

echo ""
echo "============================================"
echo "  Done! ADMIXTURE results saved to: $ADM_DIR"
echo "  Files: chr${CHR}_pruned.K.Q (ancestry fractions)"
echo "         chr${CHR}_pruned.K.P (allele frequencies)"
echo ""
echo "  Next: Rscript scripts/05_roh_analysis.R"
echo "  Visualise ADMIXTURE results in R using:"
echo "    read.table('results/admixture/chr22_pruned.5.Q')"
echo "============================================"
