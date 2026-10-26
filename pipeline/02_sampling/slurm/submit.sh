#!/usr/bin/env bash
# Submit v7 pipeline (frozen 2026-10-01):
#   preprocess -> sample -> aggregate -> CFF -> subsets -> verify -> twotier
set -e
mkdir -p logs results

JID_PRE=$(sbatch --parsable slurm/s41d_v7_preprocess.slurm)
echo "preprocess:  $JID_PRE  (15 cells, Gurobi)"

JID_SAMP=$(sbatch --parsable --dependency=afterok:$JID_PRE slurm/s41d_v7_sample.slurm)
echo "sample:      $JID_SAMP  (15 cells x 4 chains = 60 tasks, 4 x 500 x thin 300)"

JID_AGG=$(sbatch --parsable --dependency=afterok:$JID_SAMP slurm/s41d_v7_aggregate.slurm)
echo "aggregate:   $JID_AGG"

JID_CFF=$(sbatch --parsable --dependency=afterok:$JID_AGG slurm/s51_v7_cff.slurm)
echo "CFF:         $JID_CFF"

JID_SUB=$(sbatch --parsable --dependency=afterok:$JID_CFF slurm/s52_v7_subsets.slurm)
echo "subsets:     $JID_SUB"

JID_VER=$(sbatch --parsable --dependency=afterok:$JID_SUB slurm/s70_v7_verify.slurm)
echo "verify:      $JID_VER"

JID_TT=$(sbatch --parsable --dependency=afterok:$JID_VER slurm/s72_v7_twotier.slurm)
echo "twotier:     $JID_TT"

echo ""
echo "Watch: squeue -u \$USER"
