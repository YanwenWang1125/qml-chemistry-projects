#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=e1_3_fix
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%A_%a.out
# 0: E1.3 Steane search (L=2, X-part checker); 1: E1.2 Shor regression with the new rules.
SCRIPTS=(e1_3_steane_search e1_2_shor_search)
S=${SCRIPTS[$SLURM_ARRAY_TASK_ID]}
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u $S.py > ../results/${S}_output.txt 2>&1
