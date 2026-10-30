#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=e2_0_ler_b_lowp
#SBATCH --array=0-3
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%A_%a.out
# Liou [[17,1,5]] memory, setting B (default noise, BP+OSD): p in low p for the decoder slope check (plan 0.7)
PS=(0.00025 0.0005 0.00025 0.0005); MS=(Z Z X X)
P=${PS[$SLURM_ARRAY_TASK_ID]}; M=${MS[$SLURM_ARRAY_TASK_ID]}
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u e2_0_ler_b.py $P $M 200 200000000 > ../results/e2_0_ler_b_p${P}_${M}_output.txt 2>&1
