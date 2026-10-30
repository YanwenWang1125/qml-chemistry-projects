#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=a5_liou_g8
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%j.out
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u a5_liou_g8.py 400 > ../results/a5_liou_g8_output.txt 2>&1
