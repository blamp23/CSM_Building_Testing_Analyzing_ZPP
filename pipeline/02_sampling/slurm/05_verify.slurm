#!/bin/bash
#SBATCH --job-name=hopsy_v7_verify
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-15%15
#SBATCH --output=logs/v7_verify_%A_%a.out
#SBATCH --error=logs/v7_verify_%A_%a.err

set -eo pipefail
module purge
module load GCCcore/14.3.0 Python/3.13.5
source $SCRATCH/grace_v4_venv/bin/activate
source setup/env.sh

CELL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" cells_v7.txt)
echo "==== v6 VERIFY: $CELL ===="
mkdir -p results/cff_verify_v6
python -u scripts/s70_cff_verify.py \
    --cell "$CELL" \
    --models_dir v7_models \
    --pre_outdir hopsy_v7 \
    --post_outdir hopsy_loopless_v6
# also move the per-cell verify csv into the v6 dir (script writes to
# results/cff_verify/ hardcoded)
mv -f results/cff_verify/cff_verify_${CELL}.csv results/cff_verify_v6/ 2>/dev/null || true
