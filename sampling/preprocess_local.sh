#!/usr/bin/env bash
# Local PolyRound preprocess, 10 concurrent cells × 1 thread each.
set -e
cd "/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models/grace_v4"

CELLS=(BL_24 BL_48 BL_72 BL_96 BL_120 D_24 D_48 D_72 D_96 D_120 LD_24 LD_48 LD_72 LD_96 LD_120)

export GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-$HOME/gurobi.lic}"
# 10 workers × 1 thread each; leaves 4 cores for OS/MATLAB.
# Memory ~3 GB per worker × 10 ≈ 30 GB (within 32 GB).
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export GUROBI_NUM_THREADS=1

run_one () {
    cell=$1
    log=logs/preprocess_v7_local_${cell}.log
    echo "[$(date +%H:%M:%S)] start $cell -> $log"
    python -u scripts/s41d_hopsy_sample_direct.py \
        --cell "$cell" \
        --models_dir v7_models \
        --outdir_name hopsy_v7_local \
        --preprocess_only > "$log" 2>&1
    echo "[$(date +%H:%M:%S)] done  $cell"
}
export -f run_one

mkdir -p logs results/hopsy_v7_local
printf '%s\n' "${CELLS[@]}" | xargs -n1 -P10 -I{} bash -c 'run_one "$@"' _ {}

echo "[$(date +%H:%M:%S)] ALL PREPROCESS DONE"
