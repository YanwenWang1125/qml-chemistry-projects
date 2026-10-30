#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=e2_0_ler_a
#SBATCH --array=0-7
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%A_%a.out
# Setting A (simplified): Liou round under Liou's noise (gamma = beta = 1), p sweep around their 6e-5 pseudothreshold.
PS=(0.00002 0.00005 0.0001 0.0002 0.00002 0.00005 0.0001 0.0002); MS=(Z Z Z Z X X X X)
i=$SLURM_ARRAY_TASK_ID
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u e2_0_ler_a.py ${PS[$i]} ${MS[$i]} 300 > ../results/e2_0_ler_a_p${PS[$i]}_${MS[$i]}_output.txt 2>&1
