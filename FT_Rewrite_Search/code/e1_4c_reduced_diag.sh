#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=e1_4c_reduced_diag
#SBATCH --array=0-120
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x-%A_%a.out
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python e1_4c_reduced_diag.py $SLURM_ARRAY_TASK_ID ../results/e1_4c_reduced_diag
