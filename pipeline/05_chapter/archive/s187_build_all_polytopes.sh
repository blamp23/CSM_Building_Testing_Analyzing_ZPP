#!/bin/bash
# s187 -- build all 14 lipid-relay swap-corner polytopes locally.
# Serial, ~14 min each = ~3 h total. The already-cached corner (D_24 model
# + D_48 env) is detected and skipped.

set -euo pipefail
cd "/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models"
OUT=grace_v4/results/v7_final_addenda/lipid_relay_samples

CORNERS=(
    "D_24   D_48  8.48"
    "D_48   D_24  8.48"
    "LD_24  LD_48 8.48"
    "LD_48  LD_24 8.48"
    "D_48   D_72  8.36"
    "D_72   D_48  8.36"
    "LD_48  LD_72 8.43"
    "LD_72  LD_48 8.43"
    "BL_48  BL_72 8.23"
    "BL_72  BL_48 8.23"
    "D_72   D_96  10.11"
    "D_96   D_72  10.11"
    "LD_72  LD_96 11.41"
    "LD_96  LD_72 11.41"
)

I=0
for row in "${CORNERS[@]}"; do
    I=$((I+1))
    read -r MDL ENV FLOOR <<<"$row"
    echo "============ [$I/14] model=$MDL env=$ENV floor=$FLOOR ============"
    python grace_v4/scripts/s184_lipid_relay_swap_sample.py \
        --model_cell "$MDL" --env_cell "$ENV" --floor "$FLOOR" \
        --out_dir "$OUT" --preprocess_only
done
echo "============ all 14 polytopes built ============"
du -sh "$OUT/polytopes"
