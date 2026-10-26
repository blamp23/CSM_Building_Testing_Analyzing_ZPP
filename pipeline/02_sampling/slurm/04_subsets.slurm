#!/bin/bash
#SBATCH --job-name=hopsy_v7_subsets
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-15%15
#SBATCH --output=logs/v7_subsets_%A_%a.out
#SBATCH --error=logs/v7_subsets_%A_%a.err

set -eo pipefail
module purge
module load GCCcore/14.3.0 Python/3.13.5
source $SCRATCH/grace_v4_venv/bin/activate
source setup/env.sh
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

CELL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" cells_v7.txt)
echo "==== v6 SUBSETS: $CELL ===="
python -u scripts/s52_post_sample_subsets.py \
    --cell "$CELL" \
    --models_dir v7_models \
    --outdir_name subsets_v6
