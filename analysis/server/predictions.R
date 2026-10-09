# =========================================================================
# Posterior predictions for the model comparison -- the compute half
# =========================================================================
# Everything in 1_modelcomparison.qmd that needs a fitted model lives here, so
# that it can run on the cluster, where the fits are, rather than on a laptop
# that has to hold every 500 MB fit in memory at once. Sourced by
#
#   - predict_model.R (on Artemis, via `./hpc predict <model>`), which writes
#     IGC_MODELS_DIR/predictions/<model>_<illusion>.rds, and
#   - 1_modelcomparison.qmd, which reads those files and only relabels and
#     plots -- or, for a fit that has no predictions file yet, calls
#     run_predictions() itself as a fallback.
#
# What goes in here is what depends on the fit; what stays in the qmd is what
# depends on taste (labels, orders, colours, ranges drawn). The one exception
# is the prediction grid (plot_ranges and the grid sizes): it has to be decided
# before predicting, so it is here, and changing it means re-running this.
#
# Written with namespaced calls and no library() so that sourcing it on the
# cluster needs nothing attached. Needs: brms, cogmod, modelbased, insight,
# posterior, rstan, dplyr.

# Bump when the structure of what run_predictions() returns changes. The qmd
# refuses a file with another version rather than failing somewhere downstream
# on a missing column.
#
# 2 (2026-09-28): meta gains fit_brms, and run_predictions() refuses a brms
# other than the one that fitted the model (see there). Every version-1 file
# was made on a laptop with brms 2.23.1 from fits made with 2.21.0, and is
# wrong, so the bump is also what retires them.
igc_predictions_version <- 2L

# The settings every prediction file is made with. Stored in the file, so what
# a plot shows can always be traced to how it was computed. Override any of
# them per run, e.g. run_predictions(m, settings = list(ppc_ndraws = 1000)), or
# on the cluster through the IGC_PRED_* variables (see predict_model.R).
igc_prediction_settings <- function(...) {
  s <- list(
    seed = 123,
    # Posterior predictive check: draws per trial, and how many trials to
    # predict (a random subset estimates the same marginal distribution; see
    # compute_ppc()), in vectorised chunks of this many trials
    ppc_ndraws = 500,
    ppc_ntrials = 10000,
    ppc_chunk = 50,
    # Width of the bins the predicted RTs are counted in, in seconds. Fine
    # enough to be far below any density bandwidth, so the qmd can redraw the
    # density over any range without recomputing.
    ppc_bin = 0.001,
    # Parameter curves: posterior predictive draws per grid point (for the
    # error rate and mean RT), and the grid
    curves_iterations = 500,
    curves_n_strength = 31,
    curves_n_difference = 5,
    # Heatmaps: an n x n grid
    heatmap_n = 40,
    # Components run in parallel (forked, so not on Windows)
    cores = 1
  )
  utils::modifyList(s, list(...))
}


# Units ------------------------------------------------------------------

restore_units <- function(df, illusion = NULL) {
  if(!is.null(illusion)) df$Illusion_Type <- illusion
  df <- df |>
    dplyr::mutate(
      # Inverse of fit_model.R's normalisation. Difference is on [-1, 1] with 0
      # at the per-illusion mid-point difficulty (2026-09-20); a fit made before
      # that is on [0, 1] and cannot be read with this function. The bounds are
      # the observed per-illusion min/max of abs(Illusion_Difference) -- the
      # Ebbinghaus min (0.05, not 0.07) and the VerticalHorizontal max (0.30,
      # not 0.24) were wrong until 2026-09-20.
      Illusion_Difference = dplyr::case_when(
        Illusion_Type == "MullerLyer" ~ (Illusion_DifferenceZ + 1) / 2 * (0.46 - 0.04) + 0.04,
        Illusion_Type == "Ebbinghaus" ~ (Illusion_DifferenceZ + 1) / 2 * (0.70 - 0.05) + 0.05,
        Illusion_Type == "VerticalHorizontal" ~ (Illusion_DifferenceZ + 1) / 2 * (0.30 - 0.03) + 0.03,
        .default = NA
      ),
      Illusion_Strength = dplyr::case_when(
        Illusion_Type == "MullerLyer" ~ Illusion_StrengthZ * 49,
        Illusion_Type == "Ebbinghaus" ~ Illusion_StrengthZ * 2.03,
        Illusion_Type == "VerticalHorizontal" ~ Illusion_StrengthZ * 66.50,
        .default = NA
      )
  )
}

# Inverse of restore_units(): raw units -> the Z scale the models were fitted
# on, e.g. to build a prediction grid over a custom range. The bounds must match
# those of restore_units(). A lookup by name rather than case_when(), which no
# longer accepts a single `type` against a vector `x` (deprecated in dplyr
# 1.2.0); this way `type` can be one illusion or one per value of `x`, and an
# unknown illusion gives NA.
normalize_units <- function(x, type = "MullerLyer", what = "Strength") {
  what <- match.arg(what, c("Strength", "Difference"))
  if (what == "Strength") {
    scale <- c(MullerLyer = 49, Ebbinghaus = 2.03, VerticalHorizontal = 66.50)[type]
    unname(x / scale)
  } else {
    # The models only see the size of the difference, not its side
    # (fit_model.R takes abs() before normalising), so neither does this
    x <- abs(x)
    lo <- c(MullerLyer = 0.04, Ebbinghaus = 0.05, VerticalHorizontal = 0.03)[type]
    hi <- c(MullerLyer = 0.46, Ebbinghaus = 0.70, VerticalHorizontal = 0.30)[type]
    unname((x - lo) / (hi - lo) * 2 - 1)
  }
}


# Grid -------------------------------------------------------------------

# The illusion a fit was made on. A fit does not know it -- its data only has
# the Z-scored predictors -- so it is attached as an attribute by whoever loads
# the fit: the qmd's load_fit() from the file name, predict_model.R from
# models.R.
model_illusion <- function(gam) {
  illusion <- attr(gam, "illusion")
  if (is.null(illusion)) stop("The model has no illusion recorded: set attr(m, \"illusion\").")
  illusion
}

# The range of the space the curves and heatmaps are drawn over, per illusion
# and in raw units: strength signed, difference absolute (see normalize_units).
# An illusion not listed here is drawn over its full fitted range.
# TODO: custom ranges for the other illusions
plot_ranges <- list(
  MullerLyer = list(strength = c(-30, 30), difference = c(0.05, 0.45))
)

# Prediction grid of n_strength x n_difference points over that range,
# converted to the Z scale the models were fitted on
make_datagrid <- function(gam, n_strength = 31, n_difference = 5) {
  illusion <- model_illusion(gam)
  range <- plot_ranges[[illusion]]
  if (is.null(range)) {
    warning("No plot range set for ", illusion, ", using its full fitted range.")
    strength <- seq(-1, 1, length.out = n_strength)
    difference <- seq(-1, 1, length.out = n_difference)
  } else {
    strength <- normalize_units(seq(range$strength[1], range$strength[2], length.out = n_strength),
                                type = illusion, what = "Strength")
    difference <- normalize_units(seq(range$difference[1], range$difference[2], length.out = n_difference),
                                  type = illusion, what = "Difference")
  }
  insight::get_datagrid(gam, by = list(Illusion_DifferenceZ = difference,
                                       Illusion_StrengthZ = strength))
}

# A plain data frame. modelbased's outputs carry attributes (the call, the
# model's data, ...) that would otherwise be saved into the file with them.
as_plain_df <- function(d) {
  d <- as.data.frame(d)
  attributes(d) <- attributes(d)[c("names", "row.names")]
  class(d) <- "data.frame"
  d
}


# Computational cost ------------------------------------------------------

# Sampling efficiency: bulk effective sample size per hour of sampling wall
# clock, computed chain by chain so that the spread across chains is visible.
# Only the population-level parameters are counted (the smooth coefficients,
# their SDs and the group-level SDs): the thousands of participant offsets would
# otherwise dominate the average and say more about the data than the sampler.
compute_efficiency <- function(m, pars = "^b_|^bs_|^sd_|^sds_") {
  draws <- posterior::as_draws_array(m, variable = pars, regex = TRUE)
  # Time per chain. Efficiency is per hour of *sampling*, the usual
  # convention, but what one actually waits for is warmup + sampling, so both
  # are kept. This is wall clock, so it is only comparable between fits of the
  # same data with the same number of threads per chain - hence those are
  # carried along too and printed under the plot.
  time <- rstan::get_elapsed_time(m$fit)
  threads <- as.numeric(m$fit@stan_args[[1]]$num_threads)

  do.call(rbind, lapply(seq_len(posterior::nchains(draws)), function(chain) {
    # Parameters that were fixed rather than estimated (sigmadrift = 0 and the
    # like) have no effective sample size, hence the na.rm
    ess <- posterior::summarise_draws(
      posterior::subset_draws(draws, chain = chain),
      ess_bulk = posterior::ess_bulk
    )$ess_bulk
    data.frame(
      Chain = chain,
      ESS = mean(ess, na.rm = TRUE),
      Time_warmup = time[[chain, "warmup"]],
      Time_sample = time[[chain, "sample"]],
      Time_total = time[[chain, "warmup"]] + time[[chain, "sample"]],
      # get_elapsed_time() is in seconds
      ESS_per_hour = mean(ess, na.rm = TRUE) / (time[[chain, "sample"]] / 3600),
      Threads = if (length(threads) == 1) threads else NA_real_,
      N = stats::nobs(m)
    )
  }))
}


# Convergence ---------------------------------------------------------------

# Whether the chains of the combined fit agree, for the convergence table of the
# qmd. Rhat is the rank-normalised split Rhat of the posterior package, over
# every parameter, and over the population-level ones on their own: the
# participant offsets reach Rhat 1.3-1.9 in a few participants in models that
# are otherwise fine (AGENT.md 6.1), so the two answer different questions. The
# effective sample sizes are for the population-level parameters only, which is
# also what the efficiency plot counts. lp__ is kept apart because a chain
# stuck somewhere else (gam_lnr6's shard 6, gam_lba) shows in it first. The
# sampler's own diagnostics are the divergent transitions and the share of
# iterations that hit the treedepth limit.
#
# Added after version 2 of the prediction files without a bump: the qmd shows
# a model whose file has no `convergence` as not computed, and nothing else in
# it reads this.
compute_convergence <- function(m, pars = "^b_|^bs_|^sd_|^sds_") {
  draws <- posterior::as_draws_array(m)
  rhat <- posterior::summarise_draws(draws, rhat = posterior::rhat)
  rhat <- rhat[!is.na(rhat$rhat), ]
  is_pop <- grepl(pars, rhat$variable)
  is_lp <- rhat$variable == "lp__"
  all_rhat <- rhat[!is_lp, ]
  pop_rhat <- rhat[is_pop, ]

  pop <- posterior::subset_draws(draws, variable = pars, regex = TRUE)
  ess <- posterior::summarise_draws(pop, ess_bulk = posterior::ess_bulk,
                                    ess_tail = posterior::ess_tail)

  lp <- posterior::extract_variable_matrix(draws, "lp__")

  # Sampler diagnostics. Not every fit keeps them in a form brms can hand
  # back, so a failure here is a missing value rather than a failed job.
  np <- tryCatch(brms::nuts_params(m), error = function(e) NULL)
  max_treedepth <- tryCatch(m$fit@stan_args[[1]]$control$max_treedepth, error = function(e) NULL)
  if (is.null(max_treedepth)) max_treedepth <- 10
  nuts <- function(name) if (is.null(np)) NA_real_ else np$Value[np$Parameter == name]

  list(
    n_chains = posterior::nchains(draws),
    n_draws = posterior::ndraws(draws),
    n_params = nrow(all_rhat),
    n_pop_params = nrow(pop_rhat),
    rhat_pop_max = max(pop_rhat$rhat),
    rhat_pop_worst = pop_rhat$variable[which.max(pop_rhat$rhat)],
    rhat_all_max = max(all_rhat$rhat),
    rhat_all_worst = all_rhat$variable[which.max(all_rhat$rhat)],
    rhat_lp = rhat$rhat[is_lp],
    n_rhat_101_pop = sum(pop_rhat$rhat > 1.01),
    n_rhat_101_all = sum(all_rhat$rhat > 1.01),
    n_rhat_105_all = sum(all_rhat$rhat > 1.05),
    ess_bulk_pop_min = min(ess$ess_bulk, na.rm = TRUE),
    ess_tail_pop_min = min(ess$ess_tail, na.rm = TRUE),
    divergent = sum(nuts("divergent__")),
    divergent_pct = 100 * mean(nuts("divergent__")),
    treedepth_hit_pct = 100 * mean(nuts("treedepth__") >= max_treedepth),
    lp_chain_mean = colMeans(lp),
    stepsize_chain = if (is.null(np)) NA_real_ else
      tapply(nuts("stepsize__"), np$Chain[np$Parameter == "stepsize__"], mean)
  )
}


# Posterior predictive check ---------------------------------------------

# brms::posterior_predict() calls the family's sampler once per trial
# (thousands of calls of a few dozen draws each), and most of that is fixed
# cost. cogmod's posterior_predict_<family>() methods take a vector of
# observations, so the same draws come out of ~50 vectorised calls instead,
# about 3x faster for the DDM - see ?posterior_predict_cogmod_ddm.
predict_chunked <- function(m, newdata, ndraws, chunk = 50) {
  prep <- brms::prepare_predictions(m, newdata = newdata, ndraws = ndraws,
                                    allow_new_levels = TRUE)
  for (dp in names(prep$dpars)) prep$dpars[[dp]] <- brms::get_dpar(prep, dp)
  method <- get(paste0("posterior_predict_", m$family$name),
                envir = asNamespace("cogmod"))
  chunks <- split(seq_len(prep$nobs), ceiling(seq_len(prep$nobs) / chunk))
  do.call(rbind, lapply(chunks, method, prep = prep))   # (draws * N) x 2
}

# Counts per (response, bin of `bin` seconds) -- a few thousand rows instead
# of millions of RTs -- and, per response, the count and the bandwidth that
# density() would have picked on the raw values. From these the qmd redraws the
# same defective densities as from the raw RTs, over any range.
bin_rts <- function(rt, response, bin) {
  keep <- is.finite(rt) & is.finite(response)
  rt <- rt[keep]
  response <- response[keep]
  counts <- dplyr::count(data.frame(Response = response, Bin = floor(rt / bin)),
                         Response, Bin, name = "n")
  bw <- vapply(c(0, 1), function(r) {
    x <- rt[response == r]
    if (length(x) < 2) NA_real_ else stats::bw.nrd0(x)
  }, numeric(1))
  list(counts = counts,
       bw = stats::setNames(bw, c("0", "1")),
       n = c("0" = sum(response == 0), "1" = sum(response == 1)),
       bin = bin)
}

# The predicted curve is the predictive RT distribution marginal over trials,
# estimated by simulation. A random subset of trials estimates the same
# distribution as all of them - each trial contributes ndraws simulated RTs
# either way - so there is no reason to predict 320,000 trials x 50 draws
# (16 million RTs per model, minutes per model) to draw a density curve that
# 250,000 RTs already pin down. The observed histogram still uses every trial:
# it is stored as the count of each distinct (RT, Error), which is exact.
compute_ppc <- function(m, s) {
  observed <- insight::get_data(m)
  pred_obs <- if (nrow(observed) > s$ppc_ntrials) {
    observed[sample.int(nrow(observed), s$ppc_ntrials), , drop = FALSE]
  } else {
    observed
  }
  pp <- predict_chunked(m, pred_obs, s$ppc_ndraws, chunk = s$ppc_chunk)
  stopifnot(all(pp[, 2] %in% c(0, 1)))
  list(
    observed = dplyr::count(data.frame(RT = observed$RT, Error = observed$Error),
                            RT, Error, name = "n"),
    predicted = bin_rts(pp[, 1], pp[, 2], s$ppc_bin)
  )
}


# Parameter curves and heatmaps -------------------------------------------

# Predictions for each of the model's own parameters (mu and the auxiliaries,
# minus the outlier weight) over whatever grid is passed. The curves and the
# heatmaps differ only in that grid, so they share this. Each parameter is
# deterministic (summaries over all draws), so they can run in parallel.
#
# The list comes from brms, not insight::find_auxiliary(): before insight
# 1.5.4.1 that adds a "sigma" to every cogmod family, none of which has one,
# and estimate_relation() then fails on it. Parameters fixed in bf() (e.g.
# sigmabias = 0) are not in $dpars, so they are not predicted.
predict_parameters <- function(gam, datagrid, cores = 1) {
  params <- names(brms::brmsterms(stats::formula(gam))$dpars)
  params <- params[params != "poutlier"]
  out <- par_lapply(params, function(p) {
    d <- modelbased::estimate_relation(gam, data = datagrid, predict = p)
    d <- as_plain_df(d)
    d$Parameter <- p
    d
  }, cores = cores)
  do.call(rbind, out)
}

# Parameters against illusion strength, one curve per illusion difficulty. On
# the Z scale and with the family's own parameter names: the qmd relabels them
# and restores the units.
compute_curves <- function(gam, s) {
  datagrid <- make_datagrid(gam, n_strength = s$curves_n_strength,
                            n_difference = s$curves_n_difference)
  pred <- modelbased::estimate_prediction(gam, data = datagrid, centrality = "mean",
                                          iterations = s$curves_iterations,
                                          keep_iterations = TRUE)
  pred <- as_plain_df(pred)
  pred$Parameter <- ifelse(pred$Component == "response", "Error Rate", "Mean RT")
  pred$Type <- "GAM"
  pred$Component <- NULL
  pred$Row <- NULL

  # Only keep iterations of correct responses to compute mean
  iter_cols <- startsWith(names(pred), "iter_")
  is_rt <- pred$Parameter == "Mean RT"
  is_er <- pred$Parameter == "Error Rate"
  iters_rt <- pred[is_rt, iter_cols]
  iters_rt[pred[is_er, iter_cols] == 1] <- NA
  pred[is_rt, "Predicted"] <- rowMeans(iters_rt, na.rm = TRUE)
  # Remove mean RTs when only a few correct trials as it's noisy
  pred[is_rt, "Predicted"][pred[is_er, "Predicted"] > 0.9] <- NA
  pred <- pred[, !iter_cols]

  # Add the predictions for the parameters of the model itself
  params <- predict_parameters(gam, datagrid, cores = s$cores)
  params$Type <- "GAM"
  rbind(pred, params[, names(pred)])
}

# The same parameters as the curves, but over the full 2-D grid of illusion
# strength x illusion difference rather than five difficulty levels
compute_heatmaps <- function(gam, s) {
  datagrid <- make_datagrid(gam, n_strength = s$heatmap_n, n_difference = s$heatmap_n)
  predict_parameters(gam, datagrid, cores = s$cores)
}


# Everything, for one model -----------------------------------------------

# lapply, forked over `cores` where forking exists
par_lapply <- function(x, f, cores = 1) {
  if (cores > 1 && .Platform$OS.type != "windows") {
    out <- parallel::mclapply(x, f, mc.cores = min(cores, length(x)),
                              mc.preschedule = FALSE)
    failed <- vapply(out, inherits, logical(1), what = "try-error")
    if (any(failed)) stop(out[[which(failed)[1]]], call. = FALSE)
    out
  } else {
    lapply(x, f)
  }
}

# All the model-dependent results of the qmd for one fit. `fit_file` is only
# recorded (its size and time), so the qmd can tell a prediction file made
# from another fit than the one it has. Each component re-seeds, so its result
# does not depend on which others run or in what order.
#
# Only with the brms that fitted the model. brms 2.23.1 reads a fit made with
# 2.21.0 without complaint and predicts something else from the same draws --
# mu on the training rows correlates 0.67 with what 2.21.0 gives, and the
# error-rate curves come out roughly inverted along illusion difference.
# Measured 2026-09-28 on gam_rdm; the same fit under 2.21.0 on the same laptop
# matched the cluster exactly. The cause (found 2026-10-02, AGENT.md 3.9):
# 2.23's restructure() replaces an older fit's stored smooth bases with ones
# rebuilt by the local mgcv, and on the laptop's LAPACK those have
# sign-flipped columns. Nothing fails, so this check is the only thing that
# catches it.
run_predictions <- function(m, name, fit_file = NULL, settings = list()) {
  fit_brms <- as.character(m$version$brms)
  if (length(fit_brms) && utils::packageVersion("brms") != fit_brms) {
    stop(name, " was fitted with brms ", fit_brms, " but this is brms ",
         utils::packageVersion("brms"), ", which predicts such fits wrongly. ",
         "Run `./hpc predict` on the cluster, or use brms ", fit_brms, ".", call. = FALSE)
  }
  s <- do.call(igc_prediction_settings, settings)
  components <- list(
    efficiency = function(s) compute_efficiency(m),
    convergence = function(s) compute_convergence(m),
    loo = function(s) m$criteria$loo,
    ppc = function(s) compute_ppc(m, s),
    curves = function(s) compute_curves(m, s),
    heatmaps = function(s) compute_heatmaps(m, s)
  )
  # The components run in parallel, and the curves and heatmaps then fork
  # again over their parameters: give those half the cores each, which keeps
  # the total near `cores` (the other three components are one process each
  # and finish early)
  inner <- s
  inner$cores <- max(1, s$cores %/% 2)
  out <- par_lapply(names(components), function(k) {
    set.seed(s$seed)
    t0 <- Sys.time()
    res <- components[[k]](inner)
    message(sprintf("** %s: %s in %.1f min", name, k,
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    res
  }, cores = s$cores)
  names(out) <- names(components)
  if (is.null(out$loo)) warning(name, ": the fit has no loo criterion (add it at combine time).")

  c(list(meta = list(
    version = igc_predictions_version,
    name = name,
    illusion = model_illusion(m),
    family = m$family$name,
    nobs = stats::nobs(m),
    ndraws = brms::ndraws(m),
    fit_file = if (!is.null(fit_file)) basename(fit_file),
    fit_size = if (!is.null(fit_file)) file.size(fit_file),
    fit_mtime = if (!is.null(fit_file)) file.mtime(fit_file),
    created = Sys.time(),
    settings = s,
    # What the grid was built over, so the qmd can tell when plot_ranges has
    # changed since
    plot_range = plot_ranges[[model_illusion(m)]],
    packages = vapply(c("brms", "cogmod", "modelbased", "insight"),
                      \(p) as.character(utils::packageVersion(p)), character(1)),
    fit_brms = fit_brms,
    R = R.version.string
  )), out)
}
