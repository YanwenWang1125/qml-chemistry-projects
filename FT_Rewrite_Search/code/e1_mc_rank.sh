#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=e1_mc_rank
#SBATCH --array=0-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=%x-%A_%a.out
PS=(1e-3 2e-3 5e-3)
P=${PS[$SLURM_ARRAY_TASK_ID]}
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python e1_mc_rank.py $P 2e7 ../results/e1_mc_rank_p$P.json
