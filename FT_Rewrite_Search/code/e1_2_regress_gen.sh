#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=e1_2_regress_gen
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%j.out
# E1.2 regression under the new rules, dedup only (the no-dedup run OOMed in job 1158940_1).
cd ~/FTRS_gen/code
~/envs/qec/bin/python -u e1_2_shor_search.py > ../results/e1_2_shor_search_output.txt 2>&1
