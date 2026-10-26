#!/bin/bash
#SBATCH --job-name=lipid_relay
#SBATCH --array=1-56%12
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/lipid_relay_%A_%a.out
#SBATCH --error=logs/lipid_relay_%A_%a.err

# 56 tasks = 14 corners × 4 chains. Throttle to 12 concurrent per
# §11.20 convention. Pre-tar the polytope caches + scripts on local Mac;
# scp to Grace $SCRATCH/lipid_relay/; then `sbatch run_lipid_relay.slurm`.
#
# Each task reads its (model_cell, env_cell, floor, chain_id) tuple from
# the manifest below (tasks numbered 1..56).

source "$HOME/.bashrc"
module load python  # or your venv activation
export V7_REPO=$SCRATCH/lipid_relay

# ------- manifest (regenerate from s188_build_manifest.py) -----------------
#   Each line: model_cell env_cell floor chain
MANIFEST=(
  "D_24 D_48 8.48 1"  "D_24 D_48 8.48 2"  "D_24 D_48 8.48 3"  "D_24 D_48 8.48 4"
  "D_48 D_24 8.48 1"  "D_48 D_24 8.48 2"  "D_48 D_24 8.48 3"  "D_48 D_24 8.48 4"
  "LD_24 LD_48 8.48 1" "LD_24 LD_48 8.48 2" "LD_24 LD_48 8.48 3" "LD_24 LD_48 8.48 4"
  "LD_48 LD_24 8.48 1" "LD_48 LD_24 8.48 2" "LD_48 LD_24 8.48 3" "LD_48 LD_24 8.48 4"
  "D_48 D_72 8.36 1"  "D_48 D_72 8.36 2"  "D_48 D_72 8.36 3"  "D_48 D_72 8.36 4"
  "D_72 D_48 8.36 1"  "D_72 D_48 8.36 2"  "D_72 D_48 8.36 3"  "D_72 D_48 8.36 4"
  "LD_48 LD_72 8.43 1" "LD_48 LD_72 8.43 2" "LD_48 LD_72 8.43 3" "LD_48 LD_72 8.43 4"
  "LD_72 LD_48 8.43 1" "LD_72 LD_48 8.43 2" "LD_72 LD_48 8.43 3" "LD_72 LD_48 8.43 4"
  "BL_48 BL_72 8.23 1" "BL_48 BL_72 8.23 2" "BL_48 BL_72 8.23 3" "BL_48 BL_72 8.23 4"
  "BL_72 BL_48 8.23 1" "BL_72 BL_48 8.23 2" "BL_72 BL_48 8.23 3" "BL_72 BL_48 8.23 4"
  "D_72 D_96 10.11 1" "D_72 D_96 10.11 2" "D_72 D_96 10.11 3" "D_72 D_96 10.11 4"
  "D_96 D_72 10.11 1" "D_96 D_72 10.11 2" "D_96 D_72 10.11 3" "D_96 D_72 10.11 4"
  "LD_72 LD_96 11.41 1" "LD_72 LD_96 11.41 2" "LD_72 LD_96 11.41 3" "LD_72 LD_96 11.41 4"
  "LD_96 LD_72 11.41 1" "LD_96 LD_72 11.41 2" "LD_96 LD_72 11.41 3" "LD_96 LD_72 11.41 4"
)
# ------- end manifest ------------------------------------------------------

IDX=$((SLURM_ARRAY_TASK_ID - 1))
read -r MDL ENV FLOOR CHAIN <<<"${MANIFEST[$IDX]}"
echo "[task $SLURM_ARRAY_TASK_ID] model=$MDL env=$ENV floor=$FLOOR chain=$CHAIN"

python $V7_REPO/grace_v4/scripts/s184_lipid_relay_swap_sample.py \
    --model_cell "$MDL" --env_cell "$ENV" --floor "$FLOOR" \
    --out_dir "$V7_REPO/grace_v4/results/v7_final_addenda/lipid_relay_samples" \
    --chain_only "$CHAIN" --n_chains 4 --n_samples 500 --thinning 300
