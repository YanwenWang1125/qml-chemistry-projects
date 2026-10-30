#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=e2_1_shor_start
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%A_%a.out
# E2.1 part 2: Shor-type [[17,1,5]] start, single-step reachable set; 0 = jeon, 1 = fusion.
MODES=(jeon fusion)
M=${MODES[$SLURM_ARRAY_TASK_ID]}
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u e2_1_shor_start.py $M > ../results/e2_1_shor_start_${M}_output.txt 2>&1
