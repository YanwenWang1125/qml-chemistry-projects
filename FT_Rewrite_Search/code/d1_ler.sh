#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=d1_round_ler
#SBATCH --array=0-7
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%A_%a.out
# D1: setting-B LER of rounds built from the distance-5 per-gadget bests, vs Liou (e2_0_ler_b).
VS=(fusion_d5-jeon_d5 fusion_d5-jeon_d5 fusion_d5-jeon_d5 fusion_d5-jeon_d5 jeon_d5-jeon_d5 jeon_d5-jeon_d5 jeon_d5-jeon_d5 jeon_d5-jeon_d5)
PS=(0.001 0.002 0.001 0.002 0.001 0.002 0.001 0.002); MS=(Z Z X X Z Z X X)
i=$SLURM_ARRAY_TASK_ID
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u e2_7_round_ler.py ${VS[$i]} ${PS[$i]} ${MS[$i]} 400 > ../results/e2_7_round_ler_${VS[$i]}_p${PS[$i]}_${MS[$i]}_output.txt 2>&1
