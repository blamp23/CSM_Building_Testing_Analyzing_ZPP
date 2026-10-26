#!/bin/bash
#SBATCH --job-name=hopsy_v7_twotier
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-15%15
#SBATCH --output=logs/v7_twotier_%A_%a.out
#SBATCH --error=logs/v7_twotier_%A_%a.err

set -eo pipefail
module purge
module load GCCcore/14.3.0 Python/3.13.5
source $SCRATCH/grace_v4_venv/bin/activate
source setup/env.sh
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# 12 tasks = 4 hpfs (48/72/96/120) × 3 contrasts (D-LD, BL-LD, BL-D)
TASK=$SLURM_ARRAY_TASK_ID
HPFS=(48 72 96 120)
if   [[ $TASK -le 4 ]]; then A="D";  B="LD"; h=${HPFS[$((TASK-1))]}
elif [[ $TASK -le 8 ]]; then A="BL"; B="LD"; h=${HPFS[$((TASK-5))]}
else                          A="BL"; B="D";  h=${HPFS[$((TASK-9))]}
fi

REPS=results/hopsy_loopless_v6_reps
mkdir -p "$REPS"

# Pool reps for each cell (skip if already done)
for cell in "${A}_${h}" "${B}_${h}"; do
    out="$REPS/reps_pooled_${cell}.csv"
    if [[ -f "$out" ]]; then continue; fi
    python -u <<PY
import pandas as pd, numpy as np
from pathlib import Path
cell = "$cell"
memb = pd.read_csv(f"results/subsets_v6/{cell}_subset_membership.csv")
rep_ids = memb.drop_duplicates("representative").representative.tolist()
chains = []
for ch in range(1, 5):
    fp = Path(f"results/hopsy_loopless_v6/flux_samples_{cell}_ch{ch}.csv")
    if not fp.exists(): raise SystemExit(f"missing {fp}")
    df = pd.read_csv(fp, index_col=0)
    chains.append(df.reindex(rep_ids).values)
pooled = np.concatenate(chains, axis=1)
dfo = pd.DataFrame(pooled, index=rep_ids,
                   columns=[f"s{i+1}" for i in range(pooled.shape[1])])
dfo.index.name = "rxn"
dfo.to_csv("$out", float_format="%.5g")
print(f"pooled v6 {cell}: {dfo.shape}")
PY
done

mkdir -p results/two_tier_v6
OUT=results/two_tier_v6/two_tier_${A}${h}_vs_${B}${h}.csv
echo "==== v6 TWO-TIER: ${A}_${h} vs ${B}_${h} ===="
python -u scripts/two_tier_subsystem_test_fast.py \
    $REPS/reps_pooled_${A}_${h}.csv \
    $REPS/reps_pooled_${B}_${h}.csv \
    results/subsets_v6/${B}_${h}_subset_membership.csv \
    "$OUT" --ndraw 2000
