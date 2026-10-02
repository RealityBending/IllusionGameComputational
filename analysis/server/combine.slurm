#!/bin/bash
#SBATCH --job-name=igc_combine
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=192G
#SBATCH --partition=general
#SBATCH --output=combine_%j.out
#SBATCH --error=combine_%j.err
#SBATCH --mail-type=END,FAIL

# One job = one model (IGC_MODEL): merge its shards, add loo over every draw,
# keep the shards. Not an array.
#
# Memory, not CPU, is the binding constraint, and it scales with draws x rows:
# log_lik() and loo() hold several draws x 323,981 matrices at once (9.7 GiB
# each at 4,000 draws). Measured peaks for 4-shard (4,000-draw) combines:
# 100-106 GiB for the 5-6-parameter models, 122 GiB for gam_lnr6, which was
# OOM-killed at 128G on 2026-09-29 (job 11417254) and finished at 192G
# (11417275). Hence 192G (AGENT.md 4.8). Only 4 CPUs: more cores do not speed
# it up (4.8), and a slimmer request dispatches sooner. No --time is set, so
# the job takes general's 8 h maximum; loo over 4,000 draws takes 11-18 min.
# See the note in fit.slurm about the command-line overrides.

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
