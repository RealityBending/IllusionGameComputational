#!/bin/bash
#SBATCH --job-name=igc_fit
#SBATCH --array=1-4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=long
#SBATCH --time=2-00:00:00
#SBATCH --output=fit_%A_%a.out
#SBATCH --error=fit_%A_%a.err
#SBATCH --mail-type=END,FAIL

# One array task = one shard of one model (IGC_MODEL) = 2 chains x 8 threads.
# Sizing (measured, see AGENT.md 4.3/4.4): ~6.5 GB and 12-20 h per chain at
# full data. 4 tasks x 16 CPUs per model lets two models run at once under
# long's 140-CPU per-user cap.
#
# --job-name/--array/--output/--error/--mail-user/--chdir and the IGC_*
# variables are all overridden on the sbatch command line by `./hpc fit`, so
# nothing here is hard-coded to one account. Submitting this file by hand still
# works: it then falls back to the directory it was submitted from.

PROJECT_DIR="${IGC_USERS_DIR:-$SLURM_SUBMIT_DIR}"
export IGC_MODELS_DIR="${IGC_MODELS_DIR:-$PROJECT_DIR/models}"

cd "$PROJECT_DIR" || exit 1
mkdir -p "$IGC_MODELS_DIR"

echo "Running on $(hostname)"
echo "Model: ${IGC_MODEL:-<unset>}  shard: ${SLURM_ARRAY_TASK_ID:-1}"
echo "Models directory: $IGC_MODELS_DIR"
echo "Job started at: $(date '+%Y-%m-%d %H:%M:%S')"

# SLURM runs this non-interactively, where 'module' is not yet a function.
[ -n "${LMOD_CMD:-}" ] || [ ! -f /etc/profile.d/lmod.sh ] || . /etc/profile.d/lmod.sh

# Compute nodes only get /opt/ohpc/pub/modulefiles; the EasyBuild tree that
# actually holds R/CmdStanR is on the login node's MODULEPATH only.
module use "${IGC_EB_MODULES:-/mnt/shared/easybuild/modules/all}"

# CmdStanR/0.7.1-foss-2023a-R-4.3.2 brings R 4.3.2 + mgcv + brms + cmdstanr.
module load "${IGC_R_MODULE:-CmdStanR/0.7.1-foss-2023a-R-4.3.2}"

# Project library holds what the module doesn't ship (datawizard, cogmod, and
# a newer cmdstanr); it precedes the module on R_LIBS and so shadows it.
if [ -n "${IGC_R_LIBS:-}" ]; then
  export R_LIBS="${IGC_R_LIBS}${R_LIBS:+:$R_LIBS}"
fi

echo "Rscript: $(command -v Rscript)"

# Propagate R's exit status. Without this the trailing echo below makes every
# task report COMPLETED/0:0 even when Rscript died, which hides failures in
# sacct. Do not add a command after Rscript without preserving $?.
Rscript fit_model.R
status=$?

echo "FINISHED at: $(date '+%Y-%m-%d %H:%M:%S') with status $status"
exit $status
