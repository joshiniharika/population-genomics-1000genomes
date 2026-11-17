#!/bin/bash
# =============================================================================
# Script 01: Download 1000 Genomes Phase 3 Data
# Author: Niharika Joshi
# Description: Downloads chromosome 22 from 1000 Genomes Phase 3 as a
#              starting dataset for population structure analysis.
#              Using chr22 keeps file sizes manageable for demonstration.
# =============================================================================

set -euo pipefail  # Exit on error, undefined variable, or pipe failure

# --- Configuration ---
DATA_DIR="data"
PLINK_DIR="$DATA_DIR/plink_files"
BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
CHR=22

mkdir -p "$DATA_DIR" "$PLINK_DIR"

echo "============================================"
echo "  Downloading 1000 Genomes Phase 3 — Chr${CHR}"
echo "============================================"

# --- Step 1: Download VCF for chromosome 22 ---
echo "[1/3] Downloading VCF file..."
VCF_FILE="ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

if [ ! -f "$DATA_DIR/$VCF_FILE" ]; then
    wget -P "$DATA_DIR" "${BASE_URL}/${VCF_FILE}"
    wget -P "$DATA_DIR" "${BASE_URL}/${VCF_FILE}.tbi"  # tabix index
    echo "  ✓ VCF downloaded"
else
    echo "  ✓ VCF already exists — skipping download"
fi

# --- Step 2: Download population panel file ---
echo "[2/3] Downloading population panel..."
PANEL_FILE="integrated_call_samples_v3.20130502.ALL.panel"

if [ ! -f "$DATA_DIR/$PANEL_FILE" ]; then
    wget -P "$DATA_DIR" "${BASE_URL}/${PANEL_FILE}"
    echo "  ✓ Panel file downloaded"
else
    echo "  ✓ Panel file already exists — skipping"
fi

# --- Step 3: Convert VCF to PLINK binary format ---
echo "[3/3] Converting VCF to PLINK binary format (BED/BIM/FAM)..."

plink \
    --vcf "$DATA_DIR/$VCF_FILE" \
    --make-bed \
    --out "$PLINK_DIR/chr${CHR}_raw" \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:#\$1,\$2 \
    --vcf-half-call m

echo ""
echo "============================================"
echo "  Done! Raw PLINK files saved to: $PLINK_DIR"
echo "  Next step: run 02_qc_filtering.sh"
echo "============================================"
