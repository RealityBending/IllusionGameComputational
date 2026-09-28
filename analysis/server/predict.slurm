#!/bin/bash
#SBATCH --job-name=igc_predict
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=general
#SBATCH --output=predict_%j.out
#SBATCH --error=predict_%j.err
#SBATCH --mail-type=END,FAIL

# One job = one model (IGC_MODEL): posterior predictions of its combined fit,
# for the model comparison. Not an array.
#
# Unlike combine.slurm this builds no pointwise log-likelihood, so it needs a
# fraction of the memory: the fit itself, plus predictions over 10k trials and
# small grids. The 8 CPUs are used -- predictions.R forks over its components
# and over the parameters of the curves and heatmaps -- and each fork shares
# the fit copy-on-write, hence 64G rather than 8x the fit. No --time is set, so
# the job takes general's 8 h maximum. Queue it behind a combine with
#   ./hpc predict <model> --dependency=afterok:<combine job id>

PROJECT_DIR="${IGC_USERS_DIR:-$SLURM_SUBMIT_DIR}"
export IGC_MODELS_DIR="${IGC_MODELS_DIR:-$PROJECT_DIR/models}"

cd "$PROJECT_DIR" || exit 1
mkdir -p "$IGC_MODELS_DIR/predictions"

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

# Single-threaded BLAS, or the job never finishes. The module's OpenBLAS runs
# its matrix products on an OpenMP thread pool, which a forked child does not
# inherit: once a process has used it, the first %*% in any child it forks
# spins on sched_yield() forever, at 100% CPU. predictions.R forks (mclapply)
# after estimate_prediction() has done exactly that, so the curves never
# came back (2026-09-28: >65 min against ~1 min locally, where Windows cannot
# fork). The parallelism here is the forking, so BLAS threads buy nothing.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FLEXIBLAS_NUM_THREADS=1

echo "Rscript: $(command -v Rscript)"

Rscript predict_model.R
status=$?

echo "FINISHED at: $(date '+%Y-%m-%d %H:%M:%S') with status $status"
exit $status
