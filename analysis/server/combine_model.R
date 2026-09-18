# =========================================================================
# Combine the shards of ONE model into a single fit
# =========================================================================
# fit_model.R runs as an array; shard N writes <model>_<illusion>_N.rds into
# IGC_MODELS_DIR. This merges all of one model's shards into
# IGC_MODELS_DIR/combined/<model>_<illusion>.rds, adds waic, then deletes the
# shards. Which model is IGC_MODEL. Submitted by `./hpc combine <model>`.
#
# Deleting matters: fit_model.R uses file_refit = "never", so a shard left on
# disk is silently reused even if the formula or data changed. Removing shards
# once they are safely combined keeps "the file exists" equivalent to "this fit
# is finished and banked". Set IGC_KEEP_SHARDS=1 to skip the cleanup.

library(brms)
library(cmdstanr)
# cogmod is needed here too, not just to fit: add_criterion("waic") calls
# log_lik(), which for a custom family resolves log_lik_cogmod_<family>().
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

# How many draws to subsample for waic. add_criterion() errors outright if
# ndraws exceeds what the fit has ("should be between 1 and the maximum number
# of draws"), which a short test run or a part-finished array will, so ask for
# the target or everything available, whichever is smaller. Production has
# 4 shards x 2 chains x 500 = 4000 draws, so the full fit gets the 1500.
waic_draws <- function(m, target) max(1L, min(as.integer(target), brms::ndraws(m)))

# Full
out <- file.path(combined_dir, paste0(name, ".rds"))
m <- brms::combine_models(mlist = lapply(files, read_shard))
m$file <- NULL
m <- brms::add_criterion(m, "waic", ndraws = waic_draws(m, 1500)) # waic is faster than loo
saveRDS(m, out)
cat("** wrote", out, "with", brms::ndraws(m), "draws\n")

# Mini (first two shards only) -- handy for quick local inspection
if (length(files) >= 2) {
  out_mini <- file.path(combined_dir, paste0(name, "_mini.rds"))
  mini <- brms::combine_models(mlist = lapply(files[1:2], read_shard))
  mini$file <- NULL
  mini <- brms::add_criterion(mini, "waic", ndraws = waic_draws(mini, 500))
  saveRDS(mini, out_mini)
  cat("** wrote", out_mini, "with", brms::ndraws(mini), "draws\n")
} else {
  cat("** skipping mini: needs >= 2 shards\n")
}

# Only drop the shards once the combined fit is on disk and reads back.
if (nzchar(Sys.getenv("IGC_KEEP_SHARDS", unset = ""))) {
  cat("** IGC_KEEP_SHARDS set, keeping", length(files), "shards\n")
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
