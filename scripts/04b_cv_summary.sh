#!/bin/bash
# Collect cross-validation errors from fastmixture logs into one CSV
set -euo pipefail
OUT="results/admixture/cv_errors.csv"
echo "K,seed,cv_error" > "$OUT"
for f in results/admixture/fastmixture_K*_s*.log; do
    K=$(echo "$f" | sed -E 's/.*_K([0-9]+)_s([0-9]+)\.log/\1/')
    S=$(echo "$f" | sed -E 's/.*_K([0-9]+)_s([0-9]+)\.log/\2/')
    CV=$(grep "Cross-validation error" "$f" | awk '{print $4}')
    echo "$K,$S,$CV" >> "$OUT"
done
cat "$OUT"