#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=e1_3_recount
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u e1_3_recount.py > ../results/e1_3_recount_output.txt 2>&1
