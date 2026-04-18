#!/usr/bin/env bash
# ============================================================
# Cloud / Codespace setup for AS TWAS-COLOC pipeline
# Run once after creating a GitHub Codespace from this repo.
# ============================================================
set -euo pipefail

TISSUES=(
  Whole_Blood
  Adipose_Subcutaneous
  Muscle_Skeletal
  Colon_Sigmoid
  Colon_Transverse
  Small_Intestine_Terminal_Ileum
  Spleen
  Lung
)

# GTEx v8 all-pairs eQTL base URL (Google Cloud Storage, public)
EQTL_BASE="https://storage.googleapis.com/adult-gtex/bulk-gex/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL_all_associations"

# PredictDB mashr models (Zenodo)
MASHR_TAR_URL="https://zenodo.org/record/3518299/files/mashr_eqtl.tar"

echo "============================================================"
echo "AS TWAS-COLOC: Cloud Environment Setup"
echo "============================================================"

# ── 1. Python dependencies ──────────────────────────────────
echo ""
echo "[1/5] Installing Python dependencies..."
pip install -q -r requirements.txt

# ── 2. R + coloc package ────────────────────────────────────
echo ""
echo "[2/5] Installing R coloc package..."
if command -v R &>/dev/null; then
  Rscript -e 'if (!requireNamespace("coloc", quietly=TRUE)) install.packages("coloc", repos="https://cloud.r-project.org")'
  Rscript -e 'if (!requireNamespace("jsonlite", quietly=TRUE)) install.packages("jsonlite", repos="https://cloud.r-project.org")'
  echo "  R + coloc installed."
else
  echo "  WARNING: R not found. Install R first (apt install r-base) or use devcontainer."
fi

# ── 3. Clone MetaXcan (S-PrediXcan) ─────────────────────────
echo ""
echo "[3/5] Setting up MetaXcan..."
if [[ ! -d external/MetaXcan ]]; then
  mkdir -p external
  git clone --depth 1 https://github.com/hakyimlab/MetaXcan.git external/MetaXcan
  echo "  MetaXcan cloned."
else
  echo "  MetaXcan already present, skipping."
fi

# ── 4. Download PredictDB mashr models ───────────────────────
echo ""
echo "[4/5] Downloading GTEx v8 mashr models..."
MODEL_DIR="data/reference/gtex_v8_models"
mkdir -p "$MODEL_DIR"

# Check if models already exist (at least the 8 .db files we need)
MODELS_PRESENT=0
for tissue in "${TISSUES[@]}"; do
  if [[ -f "${MODEL_DIR}/gtex_v8_mashr_${tissue}.db" ]]; then
    MODELS_PRESENT=$((MODELS_PRESENT + 1))
  fi
done

if [[ $MODELS_PRESENT -lt ${#TISSUES[@]} ]]; then
  echo "  Downloading mashr models from Zenodo..."
  curl -L -o /tmp/mashr_eqtl.tar "$MASHR_TAR_URL"
  # Extract only the tissues we need (.db and .txt.gz covariance files)
  for tissue in "${TISSUES[@]}"; do
    tar -xf /tmp/mashr_eqtl.tar -C "$MODEL_DIR" --strip-components=1 \
      "eqtl/mashr/mashr_${tissue}.db" \
      "eqtl/mashr/mashr_${tissue}.txt.gz" 2>/dev/null || true
    # Rename to match our naming convention if needed
    if [[ -f "${MODEL_DIR}/mashr_${tissue}.db" ]] && [[ ! -f "${MODEL_DIR}/gtex_v8_mashr_${tissue}.db" ]]; then
      mv "${MODEL_DIR}/mashr_${tissue}.db" "${MODEL_DIR}/gtex_v8_mashr_${tissue}.db"
      mv "${MODEL_DIR}/mashr_${tissue}.txt.gz" "${MODEL_DIR}/gtex_v8_mashr_${tissue}.txt.gz"
    fi
  done
  rm -f /tmp/mashr_eqtl.tar
  echo "  Models downloaded and extracted."
else
  echo "  All ${#TISSUES[@]} tissue models already present, skipping."
fi

# ── 5. Download GTEx v8 all-pairs eQTL files ────────────────
echo ""
echo "[5/5] Downloading GTEx v8 all-pairs eQTL files..."
echo "  NOTE: Each file is ~4GB. Total ~32GB for 8 tissues."
echo "  Downloads go to: data/reference/gtex_v8_eqtl/"
EQTL_DIR="data/reference/gtex_v8_eqtl"
mkdir -p "$EQTL_DIR"

for tissue in "${TISSUES[@]}"; do
  DEST="${EQTL_DIR}/${tissue}.allpairs.txt.gz"
  if [[ -f "$DEST" ]]; then
    echo "  ${tissue}.allpairs.txt.gz already exists, skipping."
  else
    echo "  Downloading ${tissue}.allpairs.txt.gz ..."
    curl -L -o "$DEST" "${EQTL_BASE}/${tissue}.allpairs.txt.gz"
    echo "  Done: $(du -h "$DEST" | cut -f1)"
  fi
done

# ── 6. Verify liftover chain ────────────────────────────────
echo ""
CHAIN="data/reference/chains/hg19ToHg38.over.chain.gz"
if [[ -f "$CHAIN" ]]; then
  echo "Liftover chain file present."
else
  echo "WARNING: Liftover chain missing at $CHAIN"
  echo "  Download from: https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz"
  mkdir -p "data/reference/chains"
  curl -L -o "$CHAIN" "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz"
fi

echo ""
echo "============================================================"
echo "Setup complete! Run the pipeline with:"
echo "  python scripts/run_pipeline.py --no-mock 2>&1 | tee pipeline.log"
echo "============================================================"
