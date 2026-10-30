#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=e1_jeon_search
#SBATCH --array=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%A_%a.out
# E1.2 (0) and E1.3 (1) with the t = 2 checker and Jeon's literal restricted commutation (A2).
SCRIPTS=(e1_2_shor_search e1_3_steane_search)
S=${SCRIPTS[$SLURM_ARRAY_TASK_ID]}
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u $S.py 2 jeon > ../results/${S}_t2_jeon_output.txt 2>&1
