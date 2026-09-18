# Build the CmdStan precompiled header ONCE, with the exact cpp_options that
# fit_model.R uses. Run by `./hpc precompile`.
#
# Nothing to do with warm starts, and nothing to do with Stan's warmup
# iterations -- this file used to be called warmup.R, which confused all three.
#
# Why it exists: the PCH lives at a shared path inside ~/.cmdstan keyed by
# compiler version and flags (e.g. model_header_threads_nochecks_12_3.hpp.gch,
# ~800 MB). If an array of N tasks starts cold, every task races to build that
# same file and all but one die with
#   "while reading precompiled header: No such file or directory".
# The losers leave a corrupt variant behind, after which nothing compiles at
# all until cmdstanr::rebuild_cmdstan(). Compiling one trivial model first
# makes the array start warm. Run it after any change to the toolchain, the
# cmdstan version, or fit_model.R's stan_model_args.

suppressMessages(library(cmdstanr))
cat("cmdstan:", cmdstan_path(), "|", as.character(cmdstan_version()), "\n")
cat("node   :", Sys.info()[["nodename"]], "\n")

gch <- list.files(file.path(cmdstan_path(), "stan/src/stan/model"),
  pattern = "[.]gch$", full.names = TRUE
)
cat("PCH before:\n")
if (length(gch)) cat(paste0("  ", basename(gch), collapse = "\n"), "\n") else cat("  (none)\n")

stan <- write_stan_file("
data { int<lower=0> N; vector[N] y; }
parameters { real mu; }
model { y ~ normal(mu, 1); }
")

# Must match fit_model.R's stan_model_args exactly, or a different variant is
# keyed and the race comes back.
m <- cmdstan_model(stan,
  cpp_options = list(
    stan_threads = TRUE, # brms sets this whenever threads = threading(n)
    STAN_CPP_OPTIMS = TRUE,
    STAN_NO_RANGE_CHECKS = TRUE
  ),
  stanc_options = list("O1"),
  force_recompile = TRUE
)
cat("compiled OK:", basename(m$exe_file()), "\n")

gch <- list.files(file.path(cmdstan_path(), "stan/src/stan/model"),
  pattern = "[.]gch$", full.names = TRUE
)
cat("PCH after:\n")
for (g in gch) cat(sprintf("  %-45s %.0f MB\n", basename(g), file.size(g) / 1e6))
