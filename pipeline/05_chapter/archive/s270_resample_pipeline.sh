#!/bin/bash
# Full LD_72 resample under corrected bounds — one-shot driver.
# Assumes preprocess cache already exists at
# results/hopsy_v7_corrected_LD_72/LD_72_polytope_cache.pkl
# (produced by preprocess_only run of s41d_hopsy_sample_direct.py).
#
# Produces:
#   results/hopsy_v7_corrected_LD_72/flux_samples_LD_72_ch{1..4}.csv  (raw)
#   results/hopsy_v7_corrected_LD_72_loopless/flux_samples_LD_72_ch{1..4}.csv  (CFF)
#   results/pooled_v7_corrected/LD_72_pooled.csv  (final rep-level pool)

set -e
cd "/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models"
MODELS=reviewer_packet_v4_slim/v7_models
OUTDIR=hopsy_v7_corrected_LD_72
LOOPLESS=hopsy_v7_corrected_LD_72_loopless
POOLDIR=results/pooled_v7_corrected

echo "[1/4] Sampling 4 chains in parallel..."
for ch in 1 2 3 4; do
    python3 -u grace_v4/scripts/s41d_hopsy_sample_direct.py \
        --cell LD_72 --chain_only $ch \
        --n_samples 500 --thinning 300 --sampler billiard_walk \
        --models_dir $MODELS --outdir_name $OUTDIR \
        > /tmp/sample_LD72_ch${ch}.log 2>&1 &
done
wait
echo "[1/4] chains done"

echo "[2/4] CycleFreeFlux loop removal..."
# CFF uses stored lb/ub — but model needs same bounds the sampling used.
# s51b_cyclefreeflux_gurobi.py calls apply_context_bounds internally,
# which now picks up the demand_bounds table automatically.
python3 -u grace_v4/scripts/s51b_cyclefreeflux_gurobi.py \
    --cell LD_72 --outdir_in $OUTDIR --outdir_out $LOOPLESS \
    --models_dir $MODELS \
    > /tmp/cff_LD72.log 2>&1
echo "[2/4] CFF done"

echo "[3/4] Pool chains at representative level..."
mkdir -p $POOLDIR
python3 << 'PYEOF'
import pandas as pd, numpy as np
from pathlib import Path
REPO = Path(".")
LOOPLESS = REPO / "results/hopsy_v7_corrected_LD_72_loopless"
SUBSETS = REPO / "grace_v4/results/subsets_v7"
POOL = REPO / "results/pooled_v7_corrected"
POOL.mkdir(parents=True, exist_ok=True)

cell = "LD_72"
memb = pd.read_csv(SUBSETS / f"{cell}_subset_membership.csv")
rep_ids = memb.drop_duplicates("representative").representative.tolist()
parts = []
for ch in [1,2,3,4]:
    fp = LOOPLESS / f"flux_samples_{cell}_ch{ch}.csv"
    df = pd.read_csv(fp, index_col=0)
    parts.append(df.reindex(rep_ids).values)
pooled = np.concatenate(parts, axis=1)
cols = [f"ch{ch}_s{s+1}" for ch in [1,2,3,4] for s in range(500)]
dfo = pd.DataFrame(pooled, index=rep_ids, columns=cols)
dfo.index.name = "rxn"
dfo.to_csv(POOL / f"{cell}_pooled.csv", float_format="%.5g")
print(f"[write] {POOL}/{cell}_pooled.csv  shape={dfo.shape}")
PYEOF
echo "[3/4] pooled"

echo "[4/4] W1 contrast old vs new..."
python3 -u grace_v4/scripts/s270_LD72_old_vs_new_contrast.py \
    --new results/pooled_v7_corrected/LD_72_pooled.csv \
    > /tmp/contrast_LD72.log 2>&1
echo "[4/4] contrast done. Report:"
cat /tmp/contrast_LD72.log
