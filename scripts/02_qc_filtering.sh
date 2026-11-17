#!/bin/bash
# =============================================================================
# Script 02: Quality Control and Variant Filtering
# Author: Niharika Joshi
# Description: Applies standard QC filters to the raw 1000 Genomes data
#              and performs LD pruning to produce a clean SNP set suitable
#              for PCA and ADMIXTURE analysis.
# =============================================================================

set -euo pipefail

# --- Configuration ---
PLINK_DIR="data/plink_files"
QC_DIR="data/qc"
CHR=22

# QC thresholds (standard population genetics filters)
MAF=0.05          # Minor allele frequency threshold
GENO=0.05         # Maximum per-SNP missingness
MIND=0.10         # Maximum per-individual missingness
HWE=0.000001      # Hardy-Weinberg equilibrium p-value threshold
LD_WIN=50         # LD pruning window size (kb)
LD_STEP=10        # LD pruning step size
LD_R2=0.2         # LD pruning r² threshold

mkdir -p "$QC_DIR"

echo "============================================"
echo "  QC and Filtering Pipeline"
echo "============================================"

# --- Step 1: Basic QC filters ---
echo "[1/3] Applying QC filters..."
echo "      MAF > $MAF | GENO < $GENO | MIND < $MIND | HWE p > $HWE"

plink \
    --bfile "$PLINK_DIR/chr${CHR}_raw" \
    --maf $MAF \
    --geno $GENO \
    --mind $MIND \
    --hwe $HWE \
    --make-bed \
    --out "$QC_DIR/chr${CHR}_qc" \
    --allow-extra-chr

echo "  ✓ QC filters applied"

# --- Step 2: LD pruning ---
# LD pruning removes correlated SNPs so that PCA/ADMIXTURE
# reflects population structure rather than local LD patterns.
echo "[2/3] LD pruning (window=${LD_WIN}kb, step=${LD_STEP}, r²<${LD_R2})..."

plink \
    --bfile "$QC_DIR/chr${CHR}_qc" \
    --indep-pairwise $LD_WIN $LD_STEP $LD_R2 \
    --out "$QC_DIR/chr${CHR}_ldprune" \
    --allow-extra-chr

# Extract the pruned-in SNP set
plink \
    --bfile "$QC_DIR/chr${CHR}_qc" \
    --extract "$QC_DIR/chr${CHR}_ldprune.prune.in" \
    --make-bed \
    --out "$QC_DIR/chr${CHR}_pruned" \
    --allow-extra-chr

echo "  ✓ LD pruning complete"

# --- Step 3: Report SNP counts at each stage ---
echo "[3/3] Summary of filtering steps:"
echo ""

RAW_SNPS=$(wc -l < "$PLINK_DIR/chr${CHR}_raw.bim")
QC_SNPS=$(wc -l < "$QC_DIR/chr${CHR}_qc.bim")
PRUNED_SNPS=$(wc -l < "$QC_DIR/chr${CHR}_pruned.bim")
SAMPLES=$(wc -l < "$QC_DIR/chr${CHR}_pruned.fam")

echo "  Raw SNPs:       $RAW_SNPS"
echo "  After QC:       $QC_SNPS"
echo "  After pruning:  $PRUNED_SNPS"
echo "  Samples:        $SAMPLES"
echo ""
echo "============================================"
echo "  Done! Clean dataset: $QC_DIR/chr${CHR}_pruned"
echo "  Next steps:"
echo "    PCA     → Rscript scripts/03_pca_analysis.R"
echo "    ADMIXTURE → bash scripts/04_admixture.sh"
echo "    ROH     → Rscript scripts/05_roh_analysis.R"
echo "============================================"
