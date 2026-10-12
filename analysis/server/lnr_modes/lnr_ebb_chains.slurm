#!/bin/bash
#SBATCH --job-name=igc_lnr_chains
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --partition=short
#SBATCH --output=lnr_ebb_chains_%j.out
#SBATCH --error=lnr_ebb_chains_%j.err
# Diagnostic for gam_lnr_ebbinghaus: see lnr_ebb_chains.R. Submitted from the
# project dir on the cluster with IGC_R_LIBS exported.
[ -n "${LMOD_CMD:-}" ] || [ ! -f /etc/profile.d/lmod.sh ] || . /etc/profile.d/lmod.sh
module use /mnt/shared/easybuild/modules/all
module load CmdStanR/0.7.1-foss-2023a-R-4.3.2
export R_LIBS="$IGC_R_LIBS"
M=/mnt/lustre/users/psych/dmm56/IGComputational/models
Rscript lnr_ebb_chains.R "$M" ""
echo; echo "=================== OLD FIT (cogmod 0.3.3, _20261002) ==================="; echo
Rscript lnr_ebb_chains.R "$M" "_20261002"
