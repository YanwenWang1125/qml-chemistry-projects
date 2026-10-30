#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=d1_shor_start
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%j.out
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u d1_shor_start.py > ../results/d1_shor_start_output.txt 2>&1
