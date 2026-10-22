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
#   IGC_LOO_BLOCK         most draws per log_lik() block for loo (default: as
#                         many as one dense matrix holds, see loo_in_blocks())
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

# The first stanli shards (job 11603743, AGENT.md 6.6) hold a stanli_stanfit,
# stanli's subclass of rstan's stanfit, which readRDS() restores only with
# stanli on the library path. ./hpc passes IGC_STANLI_LIB under
# IGC_BACKEND=stanli. It goes last on the path, so that the brms 2.23.1 in that
# library cannot shadow the production one. read_shard() turns the fit into a
# plain stanfit, as fit_model.R now saves it, so the combined fit reads anywhere.
stanli_lib <- Sys.getenv("IGC_STANLI_LIB")
if (nzchar(stanli_lib)) .libPaths(c(.libPaths(), stanli_lib))

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

# A shard written before 2026-10-06 came from brm(file = ...), so it carries
# that path in $file (fit_model.R now saves shards itself, without it) --
# and add_criterion() writes the fit back there when it does. combine_models()
# keeps the first shard's $file, so without this the combined fit is silently
# saved *over shard 1*, and the mini below then re-reads that path and merges
# an already-combined fit with shard 2. That is how a "mini" came out with more
# draws than the full fit. Drop the slot on the way in.
read_shard <- function(f) {
  fit <- readRDS(f)
  fit$file <- NULL
  if (methods::is(fit$fit, "stanli_stanfit")) {
    fit$fit <- methods::as(fit$fit, "stanfit", strict = TRUE)
  }
  fit
}

criterion <- tolower(Sys.getenv("IGC_CRITERION", unset = "loo"))
crit_ndraws <- Sys.getenv("IGC_CRITERION_NDRAWS", unset = "")
loo_block <- Sys.getenv("IGC_LOO_BLOCK", unset = "")

# loo over every draw, with the log-likelihood computed a block of draws at a
# time.
#
# brms's log_lik() builds each participant-level predictor as one dense
# draws x rows matrix through Matrix, which cannot allocate more than 2^31 - 1
# cells. 8,000 draws x 321,113 rows is 2.6e9, and the first 8,000-draw combine
# failed there: "Cholmod error 'problem too large'" (AGENT.md 4.8). Each block
# stays under the limit, and the blocks are stacked in draw order into one
# base R matrix, which has no such limit.
#
# loo::loo() on that matrix hits the same limit one step later: loo.matrix()
# builds a draws x rows matrix of log weights through vapply(), which
# segfaults in R 4.3.2 above 2^31 cells (jobs 11608068-9, 2026-10-19). PSIS-LOO
# treats each row on its own, so the rows go through loo's function interface
# instead, which runs the same importance_sampling.matrix() and
# pointwise_loo_calcs() one row at a time and never builds that matrix.
#
# The rest is what add_criterion(m, "loo") does in brms 2.21.0 (compute_loo(),
# .loo(), prepare_loo_args()): r_eff from the chains, loo, and the response
# hash and model name that brms's loo_compare() reads back. r_eff_log_lik()
# and hash_response() are brms internals, so check them again if the
# cluster's brms moves from 2.21.0. Tested against add_criterion() on two
# smoke fits, in blocks that cut across chains (AGENT.md 4.8).
loo_in_blocks <- function(m, block) {
  nd <- brms::ndraws(m)
  ll <- NULL
  for (s in seq(1L, nd, by = block)) {
    ids <- s:min(s + block - 1L, nd)
    b <- brms::log_lik(m, draw_ids = ids)
    if (is.null(ll)) ll <- matrix(NA_real_, nd, ncol(b))
    ll[ids, ] <- b
    rm(b)
    invisible(gc())
  }
  r_eff <- brms:::r_eff_log_lik(ll, fit = m)
  out <- suppressWarnings(loo::loo(
    function(data_i, draws) draws[, data_i$i],
    data = data.frame(i = seq_len(ncol(ll))), draws = ll,
    r_eff = r_eff, save_psis = FALSE
  ))
  attr(out, "yhash") <- brms:::hash_response(m)
  # add_criterion() names it after the object it was given, which is m here.
  attr(out, "model_name") <- "m"
  out
}

# Full
out <- file.path(combined_dir, paste0(name, ".rds"))
shards <- lapply(files, read_shard)

# Every shard must come from the same cogmod: since 0.3.4 the priors depend on
# it (AGENT.md 2), and combine_models() does not compare priors. fit_model.R
# stamps a shard with its cogmod; one fitted before the stamp existed (all
# 0.3.3 or older) reads "unrecorded".
cogmod_of <- vapply(shards, function(s) {
  if (is.null(s$cogmod$version)) "unrecorded" else s$cogmod$version
}, "")
cat("** cogmod per shard:", paste0(basename(files), " ", cogmod_of, collapse = ", "), "\n")
if (length(unique(cogmod_of)) > 1) {
  stop("the shards of ", name, " were fitted with different cogmod versions (",
       paste(unique(cogmod_of), collapse = ", "), "), so on different priors. ",
       "Refit the odd ones, or combine each set in its own IGC_MODELS_DIR.", call. = FALSE)
}

# And from the same sampler: cmdstanr, or one stanli version (fit_model.R
# stamps m$stanli). The stanli tests compare samplers, so a mixed set would
# measure neither.
sampler_of <- vapply(shards, function(s) {
  if (is.null(s$stanli$version)) "cmdstanr" else paste("stanli", s$stanli$version)
}, "")
cat("** sampler per shard:", paste0(basename(files), " ", sampler_of, collapse = ", "), "\n")
if (length(unique(sampler_of)) > 1) {
  stop("the shards of ", name, " come from different samplers (",
       paste(unique(sampler_of), collapse = ", "), "). ",
       "Refit the odd ones, or combine each set in its own IGC_MODELS_DIR.", call. = FALSE)
}

# One line per chain, so that a frozen chain or a chain in another mode is
# seen here rather than as a bad Rhat after the combine (AGENT.md 4.10). A
# frozen chain has a step size orders of magnitude below the others: lp__
# barely moves, far below the rest. A chain in another mode samples normally,
# but its lp__ is off by hundreds or thousands. The first needs a restart of
# its shard. The second will not be fixed by one.
chains <- do.call(rbind, Map(function(s, f) {
  sp <- rstan::get_sampler_params(s$fit, inc_warmup = FALSE)
  lp <- rstan::extract(s$fit, pars = "lp__", permuted = FALSE)[, , 1, drop = FALSE]
  data.frame(shard = basename(f), chain = seq_along(sp),
             stepsize = vapply(sp, \(x) mean(x[, "stepsize__"]), 0),
             treedepth_hit = vapply(sp, \(x) mean(x[, "treedepth__"] >= 10), 0),
             divergent = vapply(sp, \(x) mean(x[, "divergent__"]), 0),
             lp_mean = apply(lp, 2, mean), lp_sd = apply(lp, 2, sd))
}, shards, files))
chains$lp_gap <- chains$lp_mean - stats::median(chains$lp_mean)
chains$frozen <- chains$stepsize < 1e-4 | chains$stepsize < stats::median(chains$stepsize) / 20
cat("** chains:\n")
print(format(chains, digits = 3), row.names = FALSE)
if (any(chains$frozen)) {
  warning("frozen chain(s) in ", paste(unique(chains$shard[chains$frozen]), collapse = ", "),
          ": refit that shard (AGENT.md 4.10)", call. = FALSE)
}
off <- abs(chains$lp_gap) > 5 * stats::median(chains$lp_sd) & !chains$frozen
if (any(off)) {
  warning(sum(off), " chain(s) with lp__ more than 5 SDs from the others: ",
          "likely separate modes (AGENT.md 4.10)", call. = FALSE)
}

# Every shard must also have the same warmup and number of draws. The default
# went from 500 to 1,000 draws on 2026-10-09, and fit_model.R reads it when a
# task starts, so a push while an array still has pending tasks gives that
# array shards of both lengths.
run_of <- vapply(shards, function(s) {
  sprintf("%d warmup + %d draws", s$fit@sim$warmup, s$fit@sim$iter - s$fit@sim$warmup)
}, "")
if (length(unique(run_of)) > 1) {
  stop("the shards of ", name, " differ in run length (",
       paste0(basename(files), ": ", run_of, collapse = ", "), "). ",
       "Refit the odd ones with IGC_SAMPLES/IGC_WARMUP set to match.", call. = FALSE)
}

stamp <- shards[[1]]$cogmod
m <- brms::combine_models(mlist = shards)
rm(shards)
m$file <- NULL
m$cogmod <- stamp
if (identical(criterion, "none")) {
  cat("** IGC_CRITERION=none, no criterion added\n")
} else {
  t0 <- Sys.time()
  n_obs <- nrow(m$data)
  max_cells <- .Machine$integer.max
  if (nzchar(crit_ndraws) && !identical(tolower(crit_ndraws), "all")) {
    # Subsampling: cap at what the fit has. add_criterion() errors outright
    # above it ("should be between 1 and the maximum number of draws"), which
    # a short test run or a part-finished array would hit.
    nd <- max(1L, min(as.integer(crit_ndraws), brms::ndraws(m)))
    if (as.numeric(nd) * n_obs > max_cells) {
      stop(nd, " draws x ", n_obs, " rows is over the 2^31 cells log_lik() can ",
           "hold at once. Use fewer draws, or loo over every draw.", call. = FALSE)
    }
    m <- brms::add_criterion(m, criterion, ndraws = nd)
  } else {
    # No ndraws argument at all, rather than ndraws = ndraws(m): that keeps
    # the draws in their chains, which is what loo's r_eff needs to discount
    # autocorrelation. Subsampling flattens them and r_eff degrades.
    nd <- brms::ndraws(m)
    # Equal blocks of at most max_block draws. One block is add_criterion()
    # itself, as before.
    max_block <- if (nzchar(loo_block)) as.integer(loo_block) else floor(max_cells / n_obs)
    n_blocks <- ceiling(nd / max_block)
    if (n_blocks == 1) {
      m <- brms::add_criterion(m, criterion)
    } else if (identical(criterion, "loo")) {
      block <- ceiling(nd / n_blocks)
      cat(sprintf("** log_lik() in %d blocks of <= %d draws x %d rows\n", n_blocks, block, n_obs))
      m$criteria$loo <- loo_in_blocks(m, block)
    } else {
      stop(nd, " draws x ", n_obs, " rows is over the 2^31 cells log_lik() can ",
           "hold at once, and only loo is computed in blocks. Use loo.", call. = FALSE)
    }
  }
  cat(sprintf("** %s over %d draws in %.1f min\n", criterion, nd,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  print(m$criteria[[criterion]])
}
# Written to <out>.tmp and renamed, as fit_model.R does for shards: on
# 2026-10-06 a combine died in saveRDS ("error writing to connection", Lustre
# briefly unwritable) after an hour of loo, which a direct write would leave
# as a half-written file under the combined fit's name.
tmp <- paste0(out, ".tmp")
ok <- tryCatch({
  saveRDS(m, tmp)
  file.rename(tmp, out)
}, error = function(e) {
  cat("saveRDS failed:", conditionMessage(e), "\n")
  FALSE
})
if (!isTRUE(ok)) {
  unlink(tmp)
  stop("could not save ", out, ". The shards are untouched; re-run the combine.", call. = FALSE)
}
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
