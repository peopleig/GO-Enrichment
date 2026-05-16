#!/usr/bin/env bash
set -euo pipefail

GENOME="data/hg38.fa.bgz"
PROMOTER_BED="data/promoters_500bp.bed"
PROMOTER_FA="data/promoters_500bp.fa"
DREG_OUT="data/nrf1_dreg_hits"
MOTIF="GCGC..GCGC"

echo "==> Extracting promoter sequences from genome..."
bedtools getfasta \
    -fi "$GENOME" \
    -bed "$PROMOTER_BED" \
    -name \
    -s \
    -fo "$PROMOTER_FA"

echo "==> Running dreg to scan for NRF1 motif: $MOTIF"
dreg \
    -sequence "$PROMOTER_FA" \
    -pattern "$MOTIF" \
    -outfile "${DREG_OUT}.dreg" \
    -rformat excel \
    -nowarning

echo "==> dreg finished. Output: ${DREG_OUT}.dreg"
echo "==> Next: run step3_parse_dreg_go.py"