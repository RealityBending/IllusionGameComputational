#!/bin/bash
#SBATCH --job-name=igc_stanli_install
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
# Build stanli into this account's cogmod_stanli/R_libs_<version>, on the
# release's runtime in cogmod_stanli/runtime_v<version>, as the fits load it
# (AGENT.md 6.10, 6.14). Run by `./hpc install-stanli`, which first fetches
# the release files on the login node and checks their sha256: this script
# only builds and checks. stanli_install() is never called (AGENT.md 6.6:
# it overwrites the runtime a running job has loaded). nutpieR, needed only
# for IGC_SAMPLER=nutpie, is not built here (it wants Rust >= 1.88, 6.14).
set -euo pipefail

D="${IGC_STANLI_HOME:?IGC_STANLI_HOME is not set (run this through ./hpc install-stanli)}"
V="${IGC_STANLI_VERSION:?IGC_STANLI_VERSION is not set}"
L="$D/R_libs_$V"

[ -n "${LMOD_CMD:-}" ] || [ ! -f /etc/profile.d/lmod.sh ] || . /etc/profile.d/lmod.sh
module use "${IGC_EB_MODULES:-/mnt/shared/easybuild/modules/all}"
# The fits' module, so that stanli's C bridge is built with their toolchain
# (GCCcore 12.3.0)
module load "${IGC_R_MODULE:-CmdStanR/0.7.1-foss-2023a-R-4.3.2}"
[ -z "${IGC_R_LIBS:-}" ] || export R_LIBS="${IGC_R_LIBS}"
export R_LIBS_USER="${IGC_R_LIBS_USER:-/research/cisc2/shared/cluster_R_libs/x86_64-pc-linux-gnu-library/4.2/}"
export STANLI_RUNTIME="$D/runtime_v$V/libstanli.so"
# -march=x86-64-v3 rather than the default -march=native: a build on one node
# type died with SIGILL on another (cogmod's benchmarks/stanli/hpc/run.sh)
export R_MAKEVARS_USER="$D/Makevars_$V.portable"
export MAKEFLAGS="-j${SLURM_CPUS_PER_TASK:-4}"

echo "node: $(hostname) | $(gcc --version | head -1)"
mkdir -p "$L"
R CMD INSTALL -l "$L" "$D/src_$V/stanli_$V.tar.gz"

# The check, as a fit job sees stanli: its version and runtime, a toy model
# sampled with NUTS, and the one-element `vector[1]` that fit_model.R boxes
# before unconstrain() (AGENT.md 6.19)
STANLI_NEWLIB="$L" Rscript -e '
L <- Sys.getenv("STANLI_NEWLIB")
invisible(loadNamespace("stanli", lib.loc = L))
cat("stanli", format(packageVersion("stanli", lib.loc = L)), "from", L, "\n")
cat("runtime:", stanli::stanli_runtime_path(), "| thread safe:", stanli::stanli_thread_safe(), "\n")
m <- stanli::stanli_model(code = "parameters { real mu; vector<lower=0>[1] tau; } model { mu ~ normal(1, 1); tau ~ lognormal(0, 0.5); }")
u <- stanli::unconstrain(m, list(mu = 0.1, tau = array(0.5, dim = 1)), partial = TRUE)
cat("unconstrain of a boxed vector[1]:", signif(u, 4), "\n")
f <- stanli::sample_model(m, chains = 2, warmup = 500, samples = 500, seed = 1, refresh = 0)
s <- rstan::summary(methods::as(stanli::as_stanfit(f), "stanfit", strict = TRUE))$summary
cat("NUTS: mu", round(s["mu", "mean"], 2), "(expected 1), max Rhat", round(max(s[, "Rhat"]), 3), "\n")
if (abs(s["mu", "mean"] - 1) > 0.3 || max(s[, "Rhat"]) > 1.05) stop("the toy model is off")
cat("STANLI INSTALL OK\n")
'
