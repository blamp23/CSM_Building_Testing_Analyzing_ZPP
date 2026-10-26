#!/bin/bash
#SBATCH --job-name=hopsy_v7_agg
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=1-15%15
#SBATCH --output=logs/v7_agg_%A_%a.out
#SBATCH --error=logs/v7_agg_%A_%a.err

set -eo pipefail
module purge
module load GCCcore/14.3.0 Python/3.13.5
source $SCRATCH/grace_v4_venv/bin/activate
source setup/env.sh

CELL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" cells_v7.txt)
echo "==== v6 AGGREGATE: $CELL ===="
python -u scripts/s41d_aggregate.py --cell "$CELL" --outdir_name hopsy_v7
