#!/bin/bash
#SBATCH --job-name=igc_combine
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --partition=general
#SBATCH --time=06:00:00
#SBATCH --output=combine_%j.out
#SBATCH --error=combine_%j.err
#SBATCH --mail-type=END,FAIL

# One job = one model (IGC_MODEL): merge its shards, add waic, delete them.
# Not an array.
#
# Memory, not CPU, is the binding constraint: every shard is read into RAM at
# once, and add_criterion("waic", ndraws = 1500) then builds a pointwise
# log-likelihood matrix of 1500 draws x 323,981 responses, which is ~3.9 GB
# before R makes a copy of it. Hence 128G and only 4 CPUs -- the custom
# family's log_lik is R code called once per response and is not parallelised,
# so extra cores buy nothing and a slimmer request dispatches sooner. That
# serial loop over 324k responses is also why this gets 6 h rather than a
# slot on `short`. See the note in fit.slurm about the command-line overrides.

PROJECT_DIR="${IGC_USERS_DIR:-$SLURM_SUBMIT_DIR}"
export IGC_MODELS_DIR="${IGC_MODELS_DIR:-$PROJECT_DIR/models}"

cd "$PROJECT_DIR" || exit 1
mkdir -p "$IGC_MODELS_DIR/combined"

echo "Running on $(hostname)"
echo "Model: ${IGC_MODEL:-<unset>}"
echo "Models directory: $IGC_MODELS_DIR"
echo "Job started at: $(date '+%Y-%m-%d %H:%M:%S')"

# SLURM runs this non-interactively, where 'module' is not yet a function.
[ -n "${LMOD_CMD:-}" ] || [ ! -f /etc/profile.d/lmod.sh ] || . /etc/profile.d/lmod.sh

# Compute nodes only get /opt/ohpc/pub/modulefiles; the EasyBuild tree that
# actually holds R/CmdStanR is on the login node's MODULEPATH only.
module use "${IGC_EB_MODULES:-/mnt/shared/easybuild/modules/all}"

# CmdStanR/0.7.1-foss-2023a-R-4.3.2 brings R 4.3.2 + mgcv + brms + cmdstanr.
module load "${IGC_R_MODULE:-CmdStanR/0.7.1-foss-2023a-R-4.3.2}"

if [ -n "${IGC_R_LIBS:-}" ]; then
  export R_LIBS="${IGC_R_LIBS}${R_LIBS:+:$R_LIBS}"
fi

echo "Rscript: $(command -v Rscript)"

Rscript combine_model.R
status=$?

echo "FINISHED at: $(date '+%Y-%m-%d %H:%M:%S') with status $status"
exit $status
