#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=e2_1_per_gadget
#SBATCH --array=0-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%A_%a.out
# 0: w4 jeon, 1: w8 jeon, 2: w8 fusion (BFS capped at 30000 nodes; 64G ran out at 11 min)
WS=(4 8 8); MS=(jeon jeon fusion); NS=(100000 100000 30000)
W=${WS[$SLURM_ARRAY_TASK_ID]}; M=${MS[$SLURM_ARRAY_TASK_ID]}; N=${NS[$SLURM_ARRAY_TASK_ID]}
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u e2_1_per_gadget.py $W $M $N > ../results/e2_1_per_gadget_w${W}_${M}_output.txt 2>&1
