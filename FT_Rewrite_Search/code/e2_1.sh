#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=e2_1_single_step
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%j.out
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u e2_1_single_step.py fusion > ../results/e2_1_single_step_fusion_output.txt 2>&1
