# Build / refresh the project R library on Artemis.
#
#   ./hpc install          install anything missing
#   ./hpc install cogmod   force-reinstall cogmod at IGC_COGMOD_REF's latest commit
#   ./hpc install all      force-reinstall everything
#
# The CmdStanR module already ships mgcv / brms / cmdstanr / dplyr; this only
# adds what it doesn't (datawizard, cogmod) plus their dependencies, and a
# newer cmdstanr.

lib <- Sys.getenv("IGC_R_LIBS", unset = file.path(
  "/mnt/lustre/users/psych", Sys.getenv("USER"),
  "cluster_R_libs/x86_64-pc-linux-gnu-library",
  paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
))
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
cat("target library:", lib, "\n")
cat("running on:", Sys.info()[["nodename"]], "\n")

`%||%` <- function(a, b) if (is.null(a) || !nzchar(a)) b else a # base R >= 4.4 only
repos <- "https://cloud.r-project.org"
need <- function(p) !requireNamespace(p, quietly = TRUE)

# IGC_FORCE="cogmod" (comma/space separated, or "all") reinstalls a package
# even when it is already present -- this is how cogmod is pulled up to the
# latest commit.
force <- strsplit(Sys.getenv("IGC_FORCE", unset = ""), "[, ]+")[[1]]
force <- force[nzchar(force)]
forced <- function(p) "all" %in% force || p %in% force

# Which cogmod to install, and the floor the fits require. 0.3.3 fixed the LNR
# tail gradient (chains rejecting their initial value) and the init jitter;
# fit_model.R refuses to start below it. Until 0.3.3 is merged, that is the
# dev branch -- set IGC_COGMOD_REF=main once it lands.
cogmod_ref <- Sys.getenv("IGC_COGMOD_REF", unset = "dev")
cogmod_min <- Sys.getenv("IGC_COGMOD_MIN", unset = "0.3.3")

# cmdstanr is NOT on CRAN and the CmdStanR module pins 0.7.1, which cannot read
# the CSV metadata written by CmdStan >= 2.36 (fails with
# "invalid 'y' type in 'x && y'" on save_warmup). Install a current one into the
# project library, which precedes the module on R_LIBS and so shadows it.
if (need("cmdstanr") || forced("cmdstanr") ||
  utils::packageVersion("cmdstanr") < "0.8.0") {
  cat("installing cmdstanr from stan-dev r-universe\n")
  install.packages("cmdstanr",
    lib = lib,
    repos = c("https://stan-dev.r-universe.dev", repos)
  )
}

# CmdStan itself, not just the R interface. An established account already has
# one (dmm56's is ~/.cmdstan/cmdstan-2.39.0, rebuilt 2026-09-17), but a fresh
# account has nothing and every fit dies at compile time. Building it takes
# ~15-25 min and is a once-per-account cost.
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  have_cmdstan <- tryCatch(
    {
      cmdstanr::cmdstan_version()
      TRUE
    },
    error = function(e) FALSE
  )
  if (!have_cmdstan) {
    cat("no CmdStan found -- building one (~15-25 min, once per account)
")
    cmdstanr::install_cmdstan(
      cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4")),
      overwrite = FALSE
    )
  }
  cat("cmdstan:", tryCatch(
    paste(cmdstanr::cmdstan_path(), as.character(cmdstanr::cmdstan_version())),
    error = function(e) "MISSING"
  ), "
")
}

for (p in c("remotes", "insight", "datawizard", "bayestestR")) {
  if (need(p) || forced(p)) {
    cat("installing", p, "\n")
    install.packages(p, lib = lib, repos = repos)
  }
}

# Also reinstall when the installed cogmod is below the floor, so that a run
# is never blocked by a version nobody remembered to refresh.
cogmod_version <- tryCatch(as.character(utils::packageVersion("cogmod")),
  error = function(e) NA_character_
)
if (need("cogmod") || forced("cogmod") ||
  (!is.na(cogmod_version) && package_version(cogmod_version) < cogmod_min)) {
  cat("installing cogmod@", cogmod_ref, " from GitHub (was: ",
      cogmod_version %||% "none", ")\n", sep = "")
  remotes::install_github(paste0("DominiqueMakowski/cogmod@", cogmod_ref),
    lib = lib, upgrade = "never", dependencies = TRUE,
    force = TRUE
  )
}

cat("=== FINAL CHECK ===\n")
for (p in c("mgcv", "brms", "cmdstanr", "dplyr", "datawizard", "cogmod")) {
  cat(sprintf(
    "%-11s %s", p,
    tryCatch(paste("OK", packageVersion(p)), error = function(e) "MISSING")
  ), "\n")
}
d <- tryCatch(utils::packageDescription("cogmod"), error = function(e) NULL)
if (!is.null(d)) {
  cat("cogmod source:", d$RemoteUrl %||% d$URL %||% "?", "\n")
  cat("cogmod ref   :", d$RemoteRef %||% "?", "\n")
  cat("cogmod commit:", substr(d$RemoteSha %||% "?", 1, 10), "\n")
  cat("cogmod built :", d$Built %||% "?", "\n")
}

# Fail the install rather than let a job discover this 20 hours in.
installed <- tryCatch(as.character(utils::packageVersion("cogmod")),
  error = function(e) NA_character_
)
if (is.na(installed) || package_version(installed) < cogmod_min) {
  stop("cogmod ", if (is.na(installed)) "none" else installed,
       " is below the required ", cogmod_min,
       " -- check IGC_COGMOD_REF (currently '", cogmod_ref, "')", call. = FALSE)
}
cat("cogmod >=", cogmod_min, "OK\n")
