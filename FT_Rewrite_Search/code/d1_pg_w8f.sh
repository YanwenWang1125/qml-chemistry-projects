#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=d1_per_gadget_w8f
#SBATCH --array=3
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%A_%a.out
# D1: per-gadget reachable sets from the distance-5 start, post-filtered by s + flags >= weight.
WS=(4 4 8 8); MS=(jeon fusion jeon fusion); NS=(100000 100000 100000 3000)
W=${WS[$SLURM_ARRAY_TASK_ID]}; M=${MS[$SLURM_ARRAY_TASK_ID]}; N=${NS[$SLURM_ARRAY_TASK_ID]}
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u e2_1_per_gadget.py $W $M $N d5 > ../results/e2_1_per_gadget_w${W}_${M}_d5_output.txt 2>&1
