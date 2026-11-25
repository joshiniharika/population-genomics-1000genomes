#!/bin/bash
# =============================================================================
# Script 04: Ancestry inference with fastmixture
# Author: Niharika Joshi
# Description: Runs fastmixture (ADMIXTURE-style ancestry estimation) at
#              K=2..7 with several seeds and 5-fold cross-validation.
#              ADMIXTURE 1.3.0 was tried first but segfaulted on this system,
#              so fastmixture was used instead.
# =============================================================================

set -euo pipefail

QC_DIR="data/qc"
OUT_DIR="results/admixture"
CHR=22
K_MIN=2
K_MAX=7
SEEDS="1 2 3"     # multiple seeds recommended by the fastmixture authors
THREADS=4

mkdir -p "$OUT_DIR"

for K in $(seq $K_MIN $K_MAX); do
    for SEED in $SEEDS; do
        echo "Running K=$K, seed=$SEED..."
        fastmixture \
            --bfile "$QC_DIR/chr${CHR}_pruned" \
            --K $K \
            --threads $THREADS \
            --seed $SEED \
            --cv 5 \
            --out "$OUT_DIR/chr${CHR}" \
            2>&1 | tee "$OUT_DIR/fastmixture_K${K}_s${SEED}.log"
    done
done

echo "Done. Q files: $OUT_DIR/chr${CHR}.K<k>.s<seed>.Q"
