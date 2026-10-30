#!/bin/bash
#SBATCH --account=def-ffang26_cpu
#SBATCH --job-name=a2fix_interleave
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%A_%a.out
# Item 2: corrected flag merge (shared flag read out after its last CNOT), weight-8 data order
# sorted (0) or Liou's (1); then interleave the final round (a3).
ORD=(sorted liou); O=${ORD[$SLURM_ARRAY_TASK_ID]}; T=_fix_$O
cd ~/FT_Rewrite_Search/code
~/envs/qec/bin/python -u a2_flag_start.py 400 $O $T > ../results/a2_flag_start${T}_output.txt 2>&1
~/envs/qec/bin/python -u a3_interleave.py ../data/circuits/a2_flag_final${T}.json final${T} 400 > ../results/a3_interleave_final${T}_output.txt 2>&1
~/envs/qec/bin/python -u a3_interleave.py ../data/circuits/a2_flag_start${T}.json start${T} 400 > ../results/a3_interleave_start${T}_output.txt 2>&1
