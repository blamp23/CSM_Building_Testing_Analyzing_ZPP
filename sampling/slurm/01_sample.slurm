#!/bin/bash
#SBATCH --job-name=hopsy_v7_chain
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --array=1-60%60
#SBATCH --output=logs/v7_chain_%A_%a.out
#SBATCH --error=logs/v7_chain_%A_%a.err

set -eo pipefail
module purge
module load GCCcore/14.3.0 Python/3.13.5
source $SCRATCH/grace_v4_venv/bin/activate
source setup/env.sh
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# 60 tasks = 15 cells × 4 chains
TASK=$SLURM_ARRAY_TASK_ID
CELL_INDEX=$(( (TASK - 1) / 4 + 1 ))
CHAIN_ID=$(( (TASK - 1) % 4 + 1 ))
CELL=$(sed -n "${CELL_INDEX}p" cells_v7.txt)

echo "==== v7 CHAIN: $CELL chain=$CHAIN_ID (task $TASK) ===="
echo "started $(date)"
python -u scripts/s41d_hopsy_sample_direct.py \
    --cell "$CELL" \
    --models_dir v7_models \
    --outdir_name hopsy_v7 \
    --chain_only $CHAIN_ID \
    --n_samples 500 --thinning 300 --seed 1
echo "ended $(date)"
