# =========================================================================
# Combine the shards of ONE model into a single fit
# =========================================================================
# fit_model.R runs as an array; shard N writes <model>_<illusion>_N.rds into
# IGC_MODELS_DIR. This merges all of one model's shards into
# IGC_MODELS_DIR/combined/<model>_<illusion>.rds and adds a criterion.
# Which model is IGC_MODEL. Submitted by `./hpc combine <model>`.
#
#   IGC_CRITERION         "loo" (default), "waic", or "none"
#   IGC_CRITERION_NDRAWS  a draw count, or unset / "all" (default) for all
#   IGC_DELETE_SHARDS     set to 1 to remove the shards once combined
#
# All three defaults changed on 2026-09-20, and all three are measured
# (AGENT.md 4.8): loo over every draw took 9.5 min on the 3-shard gam_lnr,
# against an 8 h wall, so there is no reason to spend accuracy on speed. Fall
# back to a draw count, or to waic, only for a model that actually overruns.
#
# Shards are KEPT now. A combined fit is minutes to rebuild from them -- a
# different criterion, a re-run after a cogmod change -- and each shard is
# ~43 h of compute, against disk that is free.
#
# But the hazard that deleting used to guard against is still live:
# fit_model.R uses file_refit = "never", so a shard left on disk is silently
# reused even when the formula or the data changed (3.6). With shards kept
# that guard is now manual -- give a changed parametrisation its own
# IGC_MODELS_DIR, or pass IGC_FILE_REFIT=always.

library(brms)
library(cmdstanr)
# cogmod is needed here too, not just to fit: add_criterion() calls log_lik(),
# which for a custom family resolves log_lik_cogmod_<family>().
library(cogmod)
library(loo)

options(mc.cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "2")))

source("models.R")

spec <- igc_model(Sys.getenv("IGC_MODEL", unset = ""))
name <- paste0(spec$name, "_", spec$illusion)

models_dir <- Sys.getenv("IGC_MODELS_DIR", unset = "models")
combined_dir <- file.path(models_dir, "combined")
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)

cat("**", name, ":", format(Sys.time()), "\n")

# Shards are named <model>_<illusion>_<shard>.rds, e.g. gam_lnr_MullerLyer_3.rds
pattern <- paste0("^", name, "_[0-9]+[.]rds$")
files <- list.files(models_dir, pattern = pattern, full.names = TRUE)
cat("** found", length(files), "shards in", models_dir, "\n")
if (length(files) == 0) {
  stop("no shards matching ", pattern, " in ", models_dir)
}

# A shard was written by brm(file = ...), so it carries that path in $file --
# and add_criterion() writes the fit back there when it does. combine_models()
# keeps the first shard's $file, so without this the combined fit is silently
# saved *over shard 1*, and the mini below then re-reads that path and merges
# an already-combined fit with shard 2. That is how a "mini" came out with more
# draws than the full fit. Drop the slot on the way in.
read_shard <- function(f) {
  fit <- readRDS(f)
  fit$file <- NULL
  fit
}

criterion <- tolower(Sys.getenv("IGC_CRITERION", unset = "loo"))
crit_ndraws <- Sys.getenv("IGC_CRITERION_NDRAWS", unset = "")

# Full
out <- file.path(combined_dir, paste0(name, ".rds"))
m <- brms::combine_models(mlist = lapply(files, read_shard))
m$file <- NULL
if (identical(criterion, "none")) {
  cat("** IGC_CRITERION=none, no criterion added\n")
} else {
  t0 <- Sys.time()
  if (nzchar(crit_ndraws) && !identical(tolower(crit_ndraws), "all")) {
    # Subsampling: cap at what the fit has. add_criterion() errors outright
    # above it ("should be between 1 and the maximum number of draws"), which
    # a short test run or a part-finished array would hit.
    nd <- max(1L, min(as.integer(crit_ndraws), brms::ndraws(m)))
    m <- brms::add_criterion(m, criterion, ndraws = nd)
  } else {
    # No ndraws argument at all, rather than ndraws = ndraws(m): that keeps
    # the draws in their chains, which is what loo's r_eff needs to discount
    # autocorrelation. Subsampling flattens them and r_eff degrades.
    nd <- brms::ndraws(m)
    m <- brms::add_criterion(m, criterion)
  }
  cat(sprintf("** %s over %d draws in %.1f min\n", criterion, nd,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
saveRDS(m, out)
cat("** wrote", out, "with", brms::ndraws(m), "draws\n")


# Shards are kept unless deletion is asked for, and then only once the
# combined fit is on disk and reads back.
if (!nzchar(Sys.getenv("IGC_DELETE_SHARDS", unset = ""))) {
  cat("** keeping", length(files), "shards (set IGC_DELETE_SHARDS=1 to drop)\n")
} else {
  readable <- tryCatch(
    {
      chk <- readRDS(out)
      inherits(chk, "brmsfit") && brms::ndraws(chk) > 0
    },
    error = function(e) FALSE
  )
  if (isTRUE(readable)) {
    removed <- file.remove(files)
    cat("** removed", sum(removed), "of", length(files), "shards\n")
  } else {
    warning("combined fit at ", out, " did not read back; keeping shards")
  }
}

cat("** finished:", name, "at", format(Sys.time()), "\n")
