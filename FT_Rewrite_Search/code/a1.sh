#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=a1_flag_merge
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%j.out
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u a1_flag_merge.py fusion_d5-jeon_d5 > ../results/a1_flag_merge_fusion_d5-jeon_d5_output.txt 2>&1
