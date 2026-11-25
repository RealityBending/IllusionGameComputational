# =========================================================================
# Posterior predictions for the model comparison -- the compute half
# =========================================================================
# Everything in 2_modelcomparison.qmd that needs a fitted model lives here, so
# that it can run on the cluster, where the fits are, rather than on a laptop
# that has to hold every 500 MB fit in memory at once. Sourced by
#
#   - predict_model.R (on Artemis, via `./hpc predict <model>`), which writes
#     IGC_MODELS_DIR/predictions/<model>_<illusion>.rds, and
#   - 2_modelcomparison.qmd, which reads those files and only relabels and
#     plots -- or, for a fit that has no predictions file yet, calls
#     run_predictions() itself as a fallback, and
#   - 3_linearmodels.qmd, which reads those of Study 1's DDM-5 and of Study
#     2's linear models (lin_*), whose grid is make_linear_datagrid() below,
#     and whose files also have the effect of the ISI (compute_isi_curves()).
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
    # error rate and the mean and median RT), and the grid. 9 difficulty levels
    # rather than 5 since 2026-11-16: the error rate and RT panels of the qmd's
    # Best Model section draw all of them, and the parameter curves every other
    # one, which are the 5 levels of before. Simulated responses per grid point:
    # 32,000 rather than 500, since the median RT rests on the correct ones
    # alone, and at 500 and still at 2000 it zigzags by more than the gap
    # between neighbouring levels wherever errors are frequent. That is more
    # than a fit's draws, so compute_curves() simulates several per draw.
    curves_iterations = 32000,
    curves_n_strength = 31,
    curves_n_difference = 9,
    # Heatmaps: an n x n grid
    heatmap_n = 40,
    # The ISI curves of a linear model (compute_isi_curves()): points along
    # each study's window, from its start (-1) to its end (1). The responses
    # per point are curves_iterations, as for the curves.
    isi_n = 21,
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

# The reliable subrange ("the reliable predictor space" in the paper): the
# range the curves and heatmaps are drawn over, per illusion and in raw units,
# strength signed and difference absolute (see normalize_units). It is the
# second of the three ranges of the predictors (AGENT.md 6.17). It lies within
# the full range (the bounds of restore_units() above: every trial, what the
# Study 1 fits cover). It contains Study 2's linear subrange (igc_linear_ranges
# in models.R: what the linear models are fitted on, still to be refined, for
# now this whole subrange). The rule (2026-11-25, the same for the three
# illusions): at each end of each predictor, leave out the most extreme values
# that together hold less than 5% of the illusion's trials, i.e. keep the
# difference from its 5th to its 95th percentile and the absolute strength up
# to its 95th percentile (quantile type 1, over every trial of the illusion).
# What is left out is, in every illusion, the levels that only the Validation
# study showed (~250 participants, ~4,000 trials per level, against 10,000 to
# 27,000 elsewhere), and what is kept is exactly the range of the main design
# (International, Suggestibility, PCS, Self-Association). Any threshold from 4%
# to 8% gives the same space. 3_linearmodels.qmd ("The Reliable Predictor
# Space") recomputes it from the data and stops if it differs from this list.
#
# The fits cover each illusion's whole range; only what is drawn is limited.
# Until 2026-11-25 only Müller-Lyer had a range here (0.05-0.45 x +-30), and
# the other two were drawn over their whole range. A prediction file made over
# a wider range is cropped to this one when read (crop_to_plot_range()).
plot_ranges <- list(
  MullerLyer = list(strength = c(-30, 30), difference = c(0.045, 0.45)),
  VerticalHorizontal = list(strength = c(-45, 45), difference = c(0.03, 0.30)),
  Ebbinghaus = list(strength = c(-1.5, 1.5), difference = c(0.05, 0.50))
)

# Whether raw (difference, strength) points are within an illusion's
# plot_ranges, with a tolerance for grid values computed in floating point
in_plot_range <- function(illusion, difference, strength) {
  r <- plot_ranges[[illusion]]
  if (is.null(r)) return(rep(TRUE, length(difference)))
  tol <- 1e-6
  abs(difference) >= r$difference[1] - tol & abs(difference) <= r$difference[2] + tol &
    strength >= r$strength[1] - tol & strength <= r$strength[2] + tol
}

# The rows of a prediction file's curves or heatmaps within plot_ranges. A
# GAM's are on the Z scale, a linear model's in raw units as well.
crop_to_plot_range <- function(df, illusion) {
  if (is.null(df) || nrow(df) == 0) return(df)
  raw <- if ("Illusion_Difference" %in% names(df)) df else restore_units(df, illusion)
  df[in_plot_range(illusion, raw$Illusion_Difference, raw$Illusion_Strength), , drop = FALSE]
}

# Whether a file's grid, made over `made` (its meta$plot_range; NULL for the
# whole range), covers the current plot_ranges of its illusion
covers_plot_range <- function(made, illusion) {
  now <- plot_ranges[[illusion]]
  if (is.null(made)) return(TRUE)
  if (is.null(now)) return(FALSE)
  tol <- 1e-6
  made$difference[1] <= now$difference[1] + tol && made$difference[2] >= now$difference[2] - tol &&
    made$strength[1] <= now$strength[1] + tol && made$strength[2] >= now$strength[2] - tol
}

# Prediction grid of n_strength x n_difference points over that range,
# converted to the Z scale the models were fitted on. A linear model of Study 2
# has predictors of its own, and a grid of its own (below).
make_datagrid <- function(gam, n_strength = 31, n_difference = 5) {
  if (!is.null(attr(gam, "linear"))) return(make_linear_datagrid(gam, n_strength, n_difference))
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

# The linear models of Study 2 (lin_* in models.R) are fitted on their own
# predictors (igc_linear_data() there): the difference rescaled from -1 to 1
# over the range fitted, the absolute strength from 0 to 1, the side of the
# illusion as Effect = +-0.5 (which the formula has only multiplied by the
# strength), and the ISI. A fit does not know its range or
# transformations: like the illusion, they are attached by whoever loads it,
# as attr(m, "linear") (the registry entry's `linear`, by predict_model.R).
#
# Its grid has the same points in raw units as the GAM's of the same illusion,
# so that the two can be drawn together, restricted to the range the linear
# model was fitted on. Each strength is predicted on its own side, and zero on
# both, so its curves are two half-lines, one per side, that meet at zero
# strength. (In the fits of the coding before 2026-11-25, with Effect and
# Difference:Effect, the model jumped there by the Effect coefficient, and the
# two half-lines need not meet.) The ISI is at the centre of its study's window
# (0), and the participant and the study at the population level (NA), as the
# GAMs' participant. Needs igc_linear_data(), so models.R sourced first.
make_linear_datagrid <- function(gam, n_strength = 31, n_difference = 5) {
  illusion <- model_illusion(gam)
  linear <- attr(gam, "linear")
  if (!exists("igc_linear_data", mode = "function")) {
    stop("make_linear_datagrid() needs igc_linear_data(): source models.R first.", call. = FALSE)
  }
  # The GAM's grid (make_datagrid() above), in raw units
  range <- plot_ranges[[illusion]]
  if (is.null(range)) {
    difference <- restore_units(data.frame(Illusion_DifferenceZ = seq(-1, 1, length.out = n_difference),
                                           Illusion_StrengthZ = 0), illusion)$Illusion_Difference
    strength <- restore_units(data.frame(Illusion_DifferenceZ = 0,
                                         Illusion_StrengthZ = seq(-1, 1, length.out = n_strength)),
                              illusion)$Illusion_Strength
  } else {
    difference <- seq(range$difference[1], range$difference[2], length.out = n_difference)
    strength <- seq(range$strength[1], range$strength[2], length.out = n_strength)
  }
  # Coded as the fit's data were. igc_linear_data() drops what is outside the
  # range fitted. Rounded before unique(): mirrored points of a seq() can
  # differ in their last bits, and would then be predicted twice.
  raw <- expand.grid(Illusion_Difference = difference,
                     Illusion_Strength = unique(round(abs(strength), 10)))
  raw$ISI <- 0
  raw$ISI_centre <- 0
  raw$Dataset <- NA
  coded <- igc_linear_data(raw, linear)
  if (nrow(coded) == 0) stop("No point of the grid is in the range of the linear model.", call. = FALSE)

  grid <- insight::get_datagrid(gam, by = list(Difference = sort(unique(coded$Difference)),
                                               Strength = sort(unique(coded$Strength)),
                                               Effect = c(-0.5, 0.5),
                                               ISI = 0))
  # The raw values of the coded ones, for restore_linear_units()
  attr(grid, "raw_units") <- list(
    difference = unique(coded[c("Difference", "Illusion_Difference")]),
    strength = unique(data.frame(Strength = coded$Strength, Illusion_Strength = abs(coded$Illusion_Strength)))
  )
  grid
}

# The raw units of a linear model's predictions, from the grid they were made
# on: Illusion_Difference, and Illusion_Strength signed by its side (Effect).
# The predictions keep the coded predictors too, and Effect tells the two
# sides of zero strength apart. Values are matched to the nearest grid value,
# which they are copies of.
restore_linear_units <- function(pred, datagrid) {
  units <- attr(datagrid, "raw_units")
  if (is.null(units)) stop("The grid has no raw units: not a make_linear_datagrid() grid.", call. = FALSE)
  nearest <- function(x, table) {
    i <- vapply(x, function(v) which.min(abs(table - v)), integer(1))
    if (any(abs(table[i] - x) > 1e-8)) stop("A prediction is not on the grid.", call. = FALSE)
    i
  }
  pred$Illusion_Difference <- units$difference$Illusion_Difference[nearest(pred$Difference, units$difference$Difference)]
  pred$Illusion_Strength <- sign(pred$Effect) *
    units$strength$Illusion_Strength[nearest(pred$Strength, units$strength$Strength)]
  pred
}

# "GAM" for the Study 1 models, "Linear" for those of Study 2, in the Type
# column of the curves
model_type <- function(gam) if (is.null(attr(gam, "linear"))) "GAM" else "Linear"

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
  # A stanli fit's stan_args do not have it (fit_model.R, "Backend"), but the
  # brm(empty = TRUE) it was built from keeps the threading it was given
  if (length(threads) == 0) threads <- as.numeric(m$threads$threads)

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
  # One simulated response per posterior draw and grid point, so one call
  # cannot simulate more responses than the fit has draws (8000 for every fit
  # so far). curves_iterations beyond that repeats the call, each time with a
  # fresh response for every draw, and pools the responses: 32,000 is 4 per
  # draw. The point summaries are then recomputed over the pool.
  per_call <- min(s$curves_iterations, brms::ndraws(gam))
  n_calls <- ceiling(s$curves_iterations / per_call)
  runs <- lapply(seq_len(n_calls), function(i) {
    n <- min(per_call, s$curves_iterations - (i - 1) * per_call)
    as_plain_df(modelbased::estimate_prediction(gam, data = datagrid, centrality = "mean",
                                                iterations = n, keep_iterations = TRUE))
  })
  iters <- do.call(cbind, lapply(runs, \(r) as.matrix(r[, startsWith(names(r), "iter_")])))
  pred <- runs[[1]][, !startsWith(names(runs[[1]]), "iter_")]
  # Every call on the same rows in the same order: the component and the
  # grid's predictors (the Z-scored two of a GAM, the coded ones of a linear
  # model)
  keys <- intersect(c("Component", names(datagrid)), names(pred))
  for (r in runs[-1]) for (k in keys) stopifnot(identical(r[[k]], pred[[k]]))
  rm(runs)
  if (n_calls > 1) {
    pred$Predicted <- rowMeans(iters)
    if ("SE" %in% names(pred)) pred$SE <- apply(iters, 1, stats::sd)
    if ("CI_low" %in% names(pred)) pred$CI_low <- apply(iters, 1, stats::quantile, probs = 0.025, names = FALSE)
    if ("CI_high" %in% names(pred)) pred$CI_high <- apply(iters, 1, stats::quantile, probs = 0.975, names = FALSE)
  }
  pred$Parameter <- ifelse(pred$Component == "response", "Error Rate", "Mean RT")
  pred$Type <- model_type(gam)
  pred$Component <- NULL
  pred$Row <- NULL

  # Only keep iterations of correct responses to compute the mean and the
  # median RT. The median was added on 2026-11-16 as rows of its own, without a
  # version bump: a file without them still has the mean, which the qmd falls
  # back on. Its SE and CIs are left empty, since modelbased's are over every
  # draw, errors included.
  is_rt <- pred$Parameter == "Mean RT"
  is_er <- pred$Parameter == "Error Rate"
  iters_rt <- iters[is_rt, , drop = FALSE]
  iters_rt[iters[is_er, , drop = FALSE] == 1] <- NA
  rm(iters)
  pred[is_rt, "Predicted"] <- rowMeans(iters_rt, na.rm = TRUE)
  median_rt <- pred[is_rt, ]
  median_rt$Parameter <- "Median RT"
  median_rt$Predicted <- apply(iters_rt, 1, stats::median, na.rm = TRUE)
  median_rt[intersect(c("SE", "CI_low", "CI_high"), names(median_rt))] <- NA_real_
  # Remove RTs when only a few correct trials as it's noisy
  too_few <- pred[is_er, "Predicted"] > 0.9
  pred[is_rt, "Predicted"][too_few] <- NA
  median_rt$Predicted[too_few] <- NA
  pred <- rbind(pred, median_rt)

  # Add the predictions for the parameters of the model itself
  params <- predict_parameters(gam, datagrid, cores = s$cores)
  params$Type <- model_type(gam)
  out <- rbind(pred, params[, names(pred)])
  # A linear model's in raw units too (a GAM's are restored by the qmd)
  if (is.null(attr(gam, "linear"))) out else restore_linear_units(out, datagrid)
}

# The same parameters as the curves, but over the full 2-D grid of illusion
# strength x illusion difference rather than five difficulty levels
compute_heatmaps <- function(gam, s) {
  datagrid <- make_datagrid(gam, n_strength = s$heatmap_n, n_difference = s$heatmap_n)
  out <- predict_parameters(gam, datagrid, cores = s$cores)
  if (is.null(attr(gam, "linear"))) out else restore_linear_units(out, datagrid)
}


# The interval before the stimulus (linear models) -----------------------------

# The effect of the ISI in a linear model of Study 2, for 3_linearmodels.qmd:
# the parameters and the observable outcomes along each study's window, at
# three stimuli. Only linear models have the ISI (a GAM gets NULL).
#
# The grid: the ISI from the start of the window (-1) to its end (1), in each
# study the fit has, for no illusion (strength 0) at the hardest, the middle
# and the easiest difficulty of the range fitted. The ISI is the same
# quadratic on every parameter at every stimulus, on the link scale, so the
# parameters' curves only move up or down from one stimulus to the next; the
# observable outcomes depend on the stimulus too. The ISI stays coded: its
# seconds are the study's window centre + 0.25 x ISI (fit_model.R), which the
# qmd adds.
make_isi_datagrid <- function(gam, n_isi = 21) {
  linear <- attr(gam, "linear")
  if (!exists("igc_linear_data", mode = "function")) {
    stop("make_isi_datagrid() needs igc_linear_data(): source models.R first.", call. = FALSE)
  }
  raw <- data.frame(Illusion_Difference = c(linear$difference[1], mean(linear$difference),
                                            linear$difference[2]),
                    Illusion_Strength = 0, ISI = 0, ISI_centre = 0, Dataset = NA)
  coded <- igc_linear_data(raw, linear)
  studies <- sort(unique(as.character(insight::get_data(gam)$Study)))
  grid <- expand.grid(Difference = coded$Difference, Strength = 0, Effect = 0,
                      ISI = seq(-1, 1, length.out = n_isi), Study = studies,
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  attr(grid, "raw_units") <- unique(coded[c("Difference", "Illusion_Difference")])
  grid
}

# Each study is predicted with its own random intercept and the participant at
# the population level: brms's re_formula = ~ (1 | Study), on every parameter.
# brms is called directly, not through modelbased as for the curves: insight
# 1.5.4 cannot pass a formula as include_random (`is.na()` of a formula fails in
# .get_predicted_args(), seen on the cluster 2026-11-25). The parameters are the
# posterior median and 95% interval at each point. The outcomes are simulated,
# curves_iterations responses per point as for the curves (all the draws, as
# many times as needed, each time with fresh responses), one point at a time
# so that nothing depends on how cogmod orders draws and observations.
compute_isi_curves <- function(gam, s) {
  if (is.null(attr(gam, "linear"))) return(NULL)
  grid <- make_isi_datagrid(gam, n_isi = s$isi_n)
  re <- ~ (1 | Study)
  # The response is never used for predictions, but brms reads it
  newdata <- grid
  newdata$RT <- 1
  newdata$Error <- 0
  keys <- grid[c("Difference", "ISI", "Study")]

  params <- names(brms::brmsterms(stats::formula(gam))$dpars)
  params <- params[params != "poutlier"]
  par_rows <- do.call(rbind, par_lapply(params, function(p) {
    draws <- brms::posterior_epred(gam, newdata = newdata, dpar = p, re_formula = re)
    data.frame(keys, Parameter = p,
               Predicted = apply(draws, 2, stats::median),
               CI_low = apply(draws, 2, stats::quantile, probs = 0.025, names = FALSE),
               CI_high = apply(draws, 2, stats::quantile, probs = 0.975, names = FALSE))
  }, cores = s$cores))

  simulate <- get(paste0("posterior_predict_", gam$family$name), envir = asNamespace("cogmod"))
  per_call <- min(s$curves_iterations, brms::ndraws(gam))
  n_calls <- ceiling(s$curves_iterations / per_call)
  sims <- lapply(seq_len(n_calls), function(k) {
    n <- min(per_call, s$curves_iterations - (k - 1) * per_call)
    prep <- brms::prepare_predictions(gam, newdata = newdata, ndraws = n, re_formula = re)
    for (dp in names(prep$dpars)) prep$dpars[[dp]] <- brms::get_dpar(prep, dp)
    # (rt, response) per simulated trial, response 1 the error boundary
    lapply(seq_len(prep$nobs), function(j) simulate(j, prep = prep))
  })
  out_rows <- do.call(rbind, lapply(seq_len(nrow(grid)), function(j) {
    x <- do.call(rbind, lapply(sims, `[[`, j))
    stopifnot(all(x[, 2] %in% c(0, 1)))
    error <- mean(x[, 2] == 1)
    rt <- x[x[, 2] == 0, 1]
    # No RT where few responses are correct, as in the curves
    few <- error > 0.9
    data.frame(keys[j, ], Parameter = c("Error Rate", "Mean RT", "Median RT"),
               Predicted = c(error, if (few) NA else mean(rt), if (few) NA else stats::median(rt)),
               CI_low = NA_real_, CI_high = NA_real_, row.names = NULL)
  }))
  rm(sims)

  curves <- rbind(out_rows, par_rows)
  curves$Type <- model_type(gam)
  units <- attr(grid, "raw_units")
  curves$Illusion_Difference <- units$Illusion_Difference[
    match(round(curves$Difference, 8), round(units$Difference, 8))]

  # The study intercepts must be in: if they were dropped, every study would
  # have the same curve. A check rather than a silent population-level file.
  mu <- curves[curves$Parameter == "mu", ]
  spread <- tapply(mu$Predicted, paste(mu$Difference, mu$ISI), function(x) diff(range(x)))
  if (length(unique(curves$Study)) > 1 && all(spread < 1e-8)) {
    stop("compute_isi_curves(): the studies' curves are identical, so their intercepts were not included.",
         call. = FALSE)
  }

  # The ISI coefficients themselves (linear and quadratic, per parameter) and
  # the SDs of the study intercepts, as posterior summaries
  draws <- posterior::as_draws_df(gam, variable = c("polyISI", "^sd_Study__"), regex = TRUE)
  coefficients <- as.data.frame(posterior::summarise_draws(
    draws, "median", ~ stats::quantile(.x, probs = c(0.025, 0.975)), "rhat", "ess_bulk"))

  list(curves = curves, coefficients = coefficients)
}


# The participants' effects (linear models) ----------------------------------

# The participant-level effects of a linear model of Study 2, for
# 4_reliability.qmd, which computes their reliability
# (performance::performance_dvour()) without loading the fit. Two tables in the
# format of modelbased::estimate_grouplevel(), one row per participant, model
# parameter and effect:
#   Component  the model parameter ("conditional" for the drift, as
#              estimate_grouplevel() has it)
#   Group      "Participant"
#   Level      the participant
#   Parameter  the effect
#   Median, MAD, CI, CI_low, CI_high   the posterior median, its MAD and its
#              95% equal-tailed interval
# with class estimate_grouplevel and coef_name "Median", which is all
# performance_dvour() reads. estimate_grouplevel() itself is not used: its
# result carries the whole model in its attributes (55 MB for the
# 50-participant test), and it summarises every random effect at once, which
# for a full fit means all ~70,000 of them in memory together. Here they are
# done one model parameter at a time, from the draws. The numbers are the same:
# checked against estimate_grouplevel() on the 50-participant test (2026-11-25).
#
#   terms    the participants' random effects as fitted (the deviations from
#            the population, `r_Participant`), every term of every parameter
#   effects  the interpretable ones, made from them draw by draw. The
#            strength is |strength| scaled from 0 to 1 over the range fitted,
#            and the difference -1 to 1, so a slope is the change over the
#            range (half of it for the difference):
#     Baseline                     Intercept: the parameter with no illusion,
#                                  at the middle difficulty of the range and
#                                  the centre of the ISI window, relative to
#                                  the participant's study (added 2026-11-25;
#                                  the only one the outlier proportion has)
#     Difficulty                   Difference: the effect of the difference
#                                  with no illusion
#     Strength (incongruent)       Strength + 0.5 Strength:Effect: the slope of
#                                  |strength| when the illusion is incongruent
#                                  (Effect = +0.5)
#     Strength (congruent)         Strength - 0.5 Strength:Effect
#     Difficulty x Strength (...)  the same two for Difference:Strength and
#                                  Difference:Strength:Effect: how the
#                                  difference changes that slope
igc_grouplevel_effects <- list(
  "Baseline" = c("Intercept" = 1),
  "Difficulty" = c("Difference" = 1),
  "Strength (incongruent)" = c("Strength" = 1, "Strength:Effect" = 0.5),
  "Strength (congruent)" = c("Strength" = 1, "Strength:Effect" = -0.5),
  "Difficulty x Strength (incongruent)" = c("Difference:Strength" = 1, "Difference:Strength:Effect" = 0.5),
  "Difficulty x Strength (congruent)" = c("Difference:Strength" = 1, "Difference:Strength:Effect" = -0.5)
)

# And one that is not a sum of random effects: the model-based point of
# subjective equality (PSE) of the drift (added 2026-11-25, at the user's
# request). It is the physical difference between the targets at which the
# participant's drift is zero (the evidence favours neither target) under an
# incongruent illusion of `igc_pse_strength` times the strongest strength of
# the range, the strength at which 4_reliability.qmd takes the behavioural
# PSE. It is computed draw by draw from the participant's *total*
# coefficients (population + study intercept + participant), at the centre of
# the ISI window:
#   a = Intercept + (Strength + 0.5 Strength:Effect) s     the drift at D = 0
#   g = Difference + (Difference:Strength + 0.5 Difference:Strength:Effect) s
#   PSE = the raw difference of the coded D = -a / g
# A gain common to all of a participant's drift coefficients cancels out of
# it, so unlike the raw strength slope it does not carry how precisely the
# participant sees differences. In raw units of the difference, like the
# behavioural PSE. Under an illusion the drift reaches zero beyond the
# hardest difference of the range for most participants, so it is
# extrapolated from the linear model, as the behavioural one is from the
# logistic.
igc_pse_strength <- 0.5
igc_grouplevel_effect_names <- c(names(igc_grouplevel_effects), "PSE")

# The inverse of each transformation in igc_transforms (models.R)
igc_transforms_inverse <- list(identity = identity, sqrt = \(x) x^2, log = exp, log1p = expm1)

compute_grouplevel <- function(m, group = "Participant") {
  if (is.null(attr(m, "linear"))) return(NULL)
  vars <- grep(paste0("^r_", group, "(__[A-Za-z]+)?[[]"), brms::variables(m), value = TRUE)
  info <- data.frame(
    variable = vars,
    Component = sub(paste0("^r_", group, "(__([A-Za-z]+))?[[].*$"), "\\2", vars),
    # The level is everything up to the last comma, the term after it
    Level = sub("^[^[]*[[](.*),[^,]*[]]$", "\\1", vars),
    Term = sub("^.*,([^,]*)[]]$", "\\1", vars)
  )
  info$Component[info$Component == ""] <- "conditional"

  # The summaries of each column of a draws x levels matrix
  summarise_columns <- function(x, levels, component, parameter) {
    data.frame(Component = component, Group = group, Level = levels, Parameter = parameter,
               Median = apply(x, 2, stats::median),
               MAD = apply(x, 2, stats::mad),
               CI = 0.95,
               CI_low = apply(x, 2, stats::quantile, probs = 0.025, names = FALSE),
               CI_high = apply(x, 2, stats::quantile, probs = 0.975, names = FALSE))
  }

  terms <- list()
  effects <- list()
  for (component in unique(info$Component)) {
    here <- info[info$Component == component, ]
    levels <- unique(here$Level)
    draws <- posterior::as_draws_matrix(m, variable = here$variable)
    # One draws x levels matrix per term, its columns in the order of `levels`
    by_term <- lapply(stats::setNames(nm = unique(here$Term)), function(term) {
      v <- here$variable[here$Term == term][match(levels, here$Level[here$Term == term])]
      unclass(draws[, v, drop = FALSE])
    })
    rm(draws)
    for (term in names(by_term)) {
      terms[[length(terms) + 1]] <- summarise_columns(by_term[[term]], levels, component, term)
    }
    for (effect in names(igc_grouplevel_effects)) {
      weights <- igc_grouplevel_effects[[effect]]
      if (!all(names(weights) %in% names(by_term))) next   # e.g. poutlier: an intercept only
      x <- Reduce(`+`, lapply(names(weights), \(t) weights[[t]] * by_term[[t]]))
      effects[[length(effects) + 1]] <- summarise_columns(x, levels, component, effect)
    }
    if (component == "conditional" && group == "Participant") {
      effects[[length(effects) + 1]] <- summarise_columns(drift_pse(m, by_term, levels), levels, component, "PSE")
    }
  }
  as_grouplevel <- function(d) {
    structure(do.call(rbind, d), class = c("estimate_grouplevel", "data.frame"), coef_name = "Median")
  }
  list(terms = as_grouplevel(terms), effects = as_grouplevel(effects))
}

# The model-based PSE (above igc_pse_strength), a draws x participants matrix,
# from the participants' drift deviations (`by_term`, draws x `levels`, as in
# compute_grouplevel()) and the model's population and study draws.
drift_pse <- function(m, by_term, levels) {
  linear <- attr(m, "linear")
  terms <- c("Intercept", "Difference", "Strength", "Strength:Effect", "Difference:Strength",
             "Difference:Strength:Effect")
  if (!all(terms %in% names(by_term))) stop("drift_pse(): the drift lacks a participant term.", call. = FALSE)
  fixed <- unclass(posterior::as_draws_matrix(m, variable = paste0("b_", terms)))
  # Each participant's total coefficient: population + participant (a vector
  # of draws added to every column of the draws x participants matrix)
  total <- function(term) by_term[[term]] + as.vector(fixed[, paste0("b_", term)])
  # ... and the intercept of the participant's study
  data <- insight::get_data(m, verbose = FALSE)
  study_of <- as.character(data$Study[match(levels, as.character(data$Participant))])
  study <- unclass(posterior::as_draws_matrix(m, variable = paste0("r_Study[", unique(study_of), ",Intercept]")))
  study <- study[, paste0("r_Study[", study_of, ",Intercept]"), drop = FALSE]

  ts <- igc_transforms[[linear$transform_strength]]
  s <- ts(igc_pse_strength * linear$strength) / ts(linear$strength)
  a <- total("Intercept") + study + (total("Strength") + 0.5 * total("Strength:Effect")) * s
  g <- total("Difference") + (total("Difference:Strength") + 0.5 * total("Difference:Strength:Effect")) * s
  coded <- -a / g
  # Back to the raw difference: the coded D runs from -1 to 1 over the
  # transformed range (igc_linear_data())
  td <- igc_transforms[[linear$transform_difference]]
  lo <- td(linear$difference[1])
  hi <- td(linear$difference[2])
  pse <- igc_transforms_inverse[[linear$transform_difference]](lo + (coded + 1) / 2 * (hi - lo))
  dimnames(pse) <- NULL
  pse
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
  # A linear model's ISI (compute_isi_curves()). Added on 2026-11-25 without a
  # version bump: only linear models have it, and the qmd checks for it.
  if (!is.null(attr(m, "linear"))) components$isi <- function(s) compute_isi_curves(m, s)
  # And its participants' effects, for 4_reliability.qmd (compute_grouplevel()),
  # also without a bump
  if (!is.null(attr(m, "linear"))) components$grouplevel <- function(s) compute_grouplevel(m)
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
    # A linear model's range and transformations (NULL for a GAM). Added
    # without a version bump: nothing reads it from a GAM's file.
    linear = attr(m, "linear"),
    packages = vapply(c("brms", "cogmod", "modelbased", "insight"),
                      \(p) as.character(utils::packageVersion(p)), character(1)),
    fit_brms = fit_brms,
    R = R.version.string
  )), out)
}
