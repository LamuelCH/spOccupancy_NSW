#!/bin/bash --login
#SBATCH --job-name=NSW_spOccupancy_modelFitting
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --qos=normal
#SBATCH --partition=general
#SBATCH --account=a_reside
#SBATCH --constraint=epyc3
#SBATCH --batch=epyc3
#SBATCH --array=3-5
#SBATCH --output=log_NSW_spOccupancy_%A_%a.out
#SBATCH --error=log_NSW_spOccupancy_%A_%a.err

----------------------------------------------------
# Load the right module (here R 4.4.0)
module load r



---------------------------------------------------
# Execution
srun Rscript scripts/S2_Model_Fitting.R $SLURM_ARRAY_TASK_ID

