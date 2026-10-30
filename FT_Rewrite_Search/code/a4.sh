#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=a4_merge_predict
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%j.out
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u a4_merge_predict.py ../data/circuits/a2_flag_start_fix_sorted.json sorted > ../results/a4_merge_predict_sorted_output.txt 2>&1
