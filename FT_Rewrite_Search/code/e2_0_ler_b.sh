#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=e2_0_ler_b
#SBATCH --array=0-5
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%A_%a.out
# Liou [[17,1,5]] memory, setting B (default noise, BP+OSD): p in {1,2,5}e-3 x memory {Z,X}
PS=(0.001 0.002 0.005 0.001 0.002 0.005); MS=(Z Z Z X X X)
P=${PS[$SLURM_ARRAY_TASK_ID]}; M=${MS[$SLURM_ARRAY_TASK_ID]}
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u e2_0_ler_b.py $P $M 400 20000000 > ../results/e2_0_ler_b_p${P}_${M}_output.txt 2>&1
