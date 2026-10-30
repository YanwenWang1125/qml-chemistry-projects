#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=a2_flag_start
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%j.out
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u a2_flag_start.py 400 > ../results/a2_flag_start_output.txt 2>&1
