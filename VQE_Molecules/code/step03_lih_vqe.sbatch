#!/bin/bash
#SBATCH --job-name=lih_vqe
#SBATCH --account=def-CHANGEME
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=slurm-%x-%j.out

# Submit from VQE_Molecules/code/:  sbatch step03_lih_vqe.sbatch
# If the job hits the time limit, submit it again: runs already in the CSV are skipped.
source ~/envs/qml/bin/activate
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
python -u step03_lih_vqe.py
