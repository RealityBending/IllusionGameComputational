#!/bin/bash
#SBATCH --job-name=igcomputational_test_combine
#SBATCH --array=1-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/mnt/lustre/scratch/psych/oc236/IGComputational/combine_%A_%a.out
#SBATCH --error=/mnt/lustre/scratch/psych/oc236/IGComputational/combine_%A_%a.err
#SBATCH --partition=long

# Set working directory
# cd /mnt/lustre/users/psych/oc236/IGComputational/models/ || exit 1

echo "Running on $(hostname)"
echo "Current directory: $(/bin/pwd)"
echo "Job started at: $(date '+%Y-%m-%d %H:%M:%S')"

module load R/4.4.1-gfbf-2023b
Rscript combine_models.R

echo "COMPLETED at: $(date '+%Y-%m-%d %H:%M:%S')"