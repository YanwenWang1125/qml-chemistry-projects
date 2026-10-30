#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=e1_5b
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u e1_5b_t2_negative_control.py > ../results/e1_5b_t2_negative_control_output.txt 2>&1
