#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=e1_4e_idle_scan
#SBATCH --array=0-20
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%A_%a.out
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python e1_4e_idle_scan.py $SLURM_ARRAY_TASK_ID ../results/e1_4e_idle_scan
