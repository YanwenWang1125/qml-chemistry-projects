#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=e1_t2_search
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%A_%a.out
# E1.2 (0) and E1.3 (1) searches with the t = 2 checker (E1.2 regression failure, 2026-10-29).
SCRIPTS=(e1_2_shor_search e1_3_steane_search)
S=${SCRIPTS[$SLURM_ARRAY_TASK_ID]}
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u $S.py 2 > ../results/${S}_t2_output.txt 2>&1
