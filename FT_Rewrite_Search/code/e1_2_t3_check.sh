#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=e1_2_t3_check
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%j.out
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u e1_2_t3_check.py > ../results/e1_2_t3_check_output.txt 2>&1
