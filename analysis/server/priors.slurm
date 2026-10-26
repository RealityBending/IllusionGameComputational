#!/bin/bash
#SBATCH --job-name=igc_priors
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --output=priors_%j.out
#SBATCH --error=priors_%j.err
#SBATCH --mail-type=FAIL

# One job for every combined fit (or those of IGC_MODEL): priors against
# posteriors, see prior_check.R and priors.R. It only reads the population-level
# draws, one fit at a time, so one CPU and the memory of one fit (~4 GB for the
# 13k-parameter fits) are enough; ~1 min a fit.

PROJECT_DIR="${IGC_USERS_DIR:-$SLURM_SUBMIT_DIR}"
export IGC_MODELS_DIR="${IGC_MODELS_DIR:-$PROJECT_DIR/models}"

cd "$PROJECT_DIR" || exit 1

echo "Running on $(hostname)"
echo "Models: ${IGC_MODEL:-all}"
echo "Job started at: $(date '+%Y-%m-%d %H:%M:%S')"

# SLURM runs this non-interactively, where 'module' is not yet a function.
[ -n "${LMOD_CMD:-}" ] || [ ! -f /etc/profile.d/lmod.sh ] || . /etc/profile.d/lmod.sh
module use "${IGC_EB_MODULES:-/mnt/shared/easybuild/modules/all}"
module load "${IGC_R_MODULE:-CmdStanR/0.7.1-foss-2023a-R-4.3.2}"
if [ -n "${IGC_R_LIBS:-}" ]; then
  export R_LIBS="${IGC_R_LIBS}${R_LIBS:+:$R_LIBS}"
fi
# brms 2.21.0, from cr291's shared library: see fit.slurm and AGENT.md 3.8.
export R_LIBS_USER="${IGC_R_LIBS_USER:-/research/cisc2/shared/cluster_R_libs/x86_64-pc-linux-gnu-library/4.2/}"
Rscript -e 'cat("brms", format(packageVersion("brms")), "from", find.package("brms"), "\n")'

Rscript prior_check.R
status=$?

echo "FINISHED at: $(date '+%Y-%m-%d %H:%M:%S') with status $status"
exit $status
