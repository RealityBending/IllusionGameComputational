# =========================================================================
# Model registry -- the one place a model is defined
# =========================================================================
# Sourced by fit_model.R (which fits exactly one of these per job) and by
# combine_model.R (which merges that model's shards). `./hpc` reads the model
# *names* straight out of this file, so keep the declaration lines in the form
#
#     <name> = list(
#
# at two-space indent, one per model, and use that form for nothing else.
# `./hpc models` prints what it found.
#
# Two lists:
#   igc_formulas  one brms bf() per model family, as a function so that
#                 sourcing this file costs nothing and only the model actually
#                 being fitted is built.
#   igc_models    the registry: one entry per family x illusion, named
#                 <family>_<illusion in lower case>, e.g. gam_lnr_mullerlyer.
#                 Each entry is
#                   illusion  the Illusion_Type the model is fitted to
#                   formula   the family's function from igc_formulas
#                   linear    (Study 2's linear models only) the range of
#                             stimuli fitted (its illusion's linear subrange,
#                             igc_linear_ranges below) and the
#                             transformations of the predictors, which
#                             igc_linear_data() below turns into the model's
#                             data
#
# Until 2026-10-13 the MullerLyer models were named after the bare family
# (gam_lnr, gam_lba, ...) and held the formulas, which the other illusions
# borrowed. Older notes, job logs and the fits in models/old_20261013/ use
# those names: gam_lnr_MullerLyer.rds there is gam_lnr_mullerlyer_MullerLyer.rds
# here.

# Every distributional parameter gets the same 2-D tensor smooth of illusion
# difference x illusion strength, plus a participant intercept. Built from a
# string so the five identical blocks are written once; the environment is the
# global one, as it would be for a formula typed at the top of a script.
igc_t2 <- function(lhs) {
  stats::as.formula(
    paste0(
      lhs, " ~ t2(Illusion_DifferenceZ, Illusion_StrengthZ, ",
      "k = c(5, 5), bs = c('cr', 'cr')) + (1 | Participant)"
    ),
    env = globalenv()
  )
}

# Response side of every formula: RT, with the accuracy coding brms needs for
# the two-accumulator families.
igc_rt <- function() igc_t2("RT | dec(Error)")

# Study 2's linear models (analysis/3_linearmodels.qmd): every distributional
# parameter gets the same linear effects of difficulty x absolute strength,
# with a slope of its own on each side of the illusion, a quadratic in the ISI,
# all of them varying by participant, and an intercept varying by study. The
# side (Effect) enters only multiplied by Strength, so the two sides meet at
# zero strength: Difference * Strength * Effect without Effect and
# Difference:Effect, which would let them jump there (dropped 2026-11-25,
# AGENT.md 6.15). The predictors are made by igc_linear_data() below.
igc_lin <- function(lhs) {
  terms <- "Difference * Strength + Strength:Effect + Difference:Strength:Effect"
  stats::as.formula(
    paste0(
      lhs, " ~ ", terms, " + poly(ISI, 2, raw = TRUE) + ",
      "(1 + ", terms, " | Participant) + (1 | Study)"
    ),
    env = globalenv()
  )
}


igc_formulas <- list(

  # LNR ---------------------------------------------------------------------
  # sigmabias = 0 is the plain LNR. Omitting it does not error: it silently
  # becomes a freely estimated start-point range. Watch for this class of bug
  # whenever cogmod adds a distributional parameter.
  gam_lnr = function() {
    brms::bf(
      igc_rt(),
      igc_t2("nuone"),
      igc_t2("sigmazero"),
      igc_t2("sigmaone"),
      sigmabias = 0,
      igc_t2("ndt"),
      poutlier ~ 1 + (1 | Participant),
      family = cogmod::cogmod_lnr()
    )
  },

  # LNR-6 -------------------------------------------------------------------
  # The LNR above plus a freely estimated sigmabias -- the between-trial
  # start-point range -- smoothed like every other parameter. Six smoothed
  # dpars against the plain LNR's five, hence the name.
  #
  # Expect it to cost more per gradient than gam_lnr, not less: with
  # sigmabias > 0 cogmod takes the erfc-based ldiff_Phi path for the tails
  # rather than the single-tail shortcut, which cogmod 0.3.3 measured at about
  # 15% dearer per gradient (and 20% cheaper for the sigmabias = 0 case).
  gam_lnr6 = function() {
    brms::bf(
      igc_rt(),
      igc_t2("nuone"),
      igc_t2("sigmazero"),
      igc_t2("sigmaone"),
      igc_t2("sigmabias"),
      igc_t2("ndt"),
      poutlier ~ 1 + (1 | Participant),
      family = cogmod::cogmod_lnr()
    )
  },

  # DDM-4 -------------------------------------------------------------------
  # The four-parameter DDM: drift, boundary, bias, ndt. The three between-trial
  # variabilities are fixed at 0.
  gam_ddm4 = function() {
    brms::bf(
      igc_rt(),
      igc_t2("boundary"),
      igc_t2("bias"),
      igc_t2("ndt"),
      sigmadrift = 0,
      sigmabias = 0,
      sigmandt = 0,
      poutlier ~ 1 + (1 | Participant),
      family = cogmod::cogmod_ddm()
    )
  },

  # DDM-5 -------------------------------------------------------------------
  # DDM-4 plus between-trial drift variability, itself smoothed. One more 2-D
  # smooth than the others, so expect it to be the straggler; give it its own
  # --time if it is ever submitted.
  gam_ddm5 = function() {
    brms::bf(
      igc_rt(),
      igc_t2("boundary"),
      igc_t2("bias"),
      igc_t2("ndt"),
      igc_t2("sigmadrift"),
      sigmabias = 0,
      sigmandt = 0,
      poutlier ~ 1 + (1 | Participant),
      family = cogmod::cogmod_ddm()
    )
  },

  # DDM-7 -------------------------------------------------------------------
  # The full DDM: nothing fixed. Every one of cogmod_ddm()'s seven parameters
  # gets its own 2-D smooth, including all three between-trial variabilities.
  #
  # DO NOT SUBMIT THIS AT FULL DATA (measured 2026-09-18, AGENT.md 4.7/4.7.1).
  # At 30 participants it did not reach iteration 100 in 50 minutes, where
  # every other model in the registry finished all 400 in 6-11.
  #
  # The reason is NOT geometry, which is what this comment said first and what
  # the smooth-heavy formula invites you to assume. cogmod_ddm_decision_lpdf()
  # branches on `sw == 0 && sigmandt == 0`: estimating sigmabias or sigmandt
  # (as opposed to fixing them at 0, as gam_ddm4 and gam_ddm5 do) leaves the
  # analytic Wiener density for Stan's adaptive numerical quadrature, at 18x
  # the per-gradient cost with one freed and 55x with both. sigmadrift is
  # analytic and costs 2.8x, which is the whole reason gam_ddm5 is affordable.
  #
  # So: tighter priors do not rescue it (26x is the floor, because the test is
  # for EXACT zero); dropping the smooths on the variabilities does not either
  # (the cost is per observation, not per coefficient); and neither does a warm
  # start (it buys back warmup, not per-gradient cost). What does work is fewer
  # participants -- about 2.7-4.7 days per chain at 200, with 1 chain x 16
  # threads and its own IGC_MODELS_DIR. The real fix is a fixed-node quadrature
  # rule in cogmod; see AGENT.md 4.7.1.
  gam_ddm7 = function() {
    brms::bf(
      igc_rt(),
      igc_t2("boundary"),
      igc_t2("bias"),
      igc_t2("ndt"),
      igc_t2("sigmadrift"),
      igc_t2("sigmabias"),
      igc_t2("sigmandt"),
      poutlier ~ 1 + (1 | Participant),
      family = cogmod::cogmod_ddm()
    )
  },

  # RDM ---------------------------------------------------------------------
  # Racing diffusion: two Wald accumulators, one per response. mu is the drift
  # of accumulator 0 and driftone that of accumulator 1, so it has the same
  # two-drift shape as the LNR.
  #
  # sigmabias (the maximum starting point, A) is fixed at 0: that is the plain
  # racing diffusion, both accumulators starting from the same place. gam_rdm5
  # below is the same model with it freed. cogmod 0.3.3 also made this family's
  # gradient exact, where every normal tail previously went through
  # std_normal_lcdf()'s approximate partials; in a race those partials *are*
  # the gradient of the drifts.
  gam_rdm = function() {
    brms::bf(
      igc_rt(),
      igc_t2("driftone"),
      igc_t2("boundary"),
      igc_t2("ndt"),
      sigmabias = 0,
      poutlier ~ 1 + (1 | Participant),
      family = cogmod::cogmod_rdm()
    )
  },

  # RDM-5 -------------------------------------------------------------------
  # The RDM with the start-point range freed and smoothed like everything else,
  # so all five of cogmod_rdm()'s parameters are estimated. Stands to gam_rdm
  # as gam_lnr6 does to gam_lnr, and the waic comparison against gam_rdm is
  # what says whether the start-point range is identified on this data.
  #
  # Unlike the LBA, the RDM needs no scaling constraint -- its diffusion
  # coefficient is fixed internally -- so freeing sigmabias here does not open
  # the ridge that freeing sigmazero would in gam_lba.
  gam_rdm5 = function() {
    brms::bf(
      igc_rt(),
      igc_t2("driftone"),
      igc_t2("sigmabias"),
      igc_t2("boundary"),
      igc_t2("ndt"),
      poutlier ~ 1 + (1 | Participant),
      family = cogmod::cogmod_rdm()
    )
  },

  # RDM-5r ------------------------------------------------------------------
  # gam_rdm5 with the start-point range freed but not smoothed: a population
  # intercept and a participant intercept, no effect of the illusion. Added
  # 2026-11-15 because gam_rdm5 split into two stationary modes again on the
  # 2026-10-13 refit, 6 chains against 2 on both MullerLyer and
  # VerticalHorizontal (AGENT.md 4.10 B, option 1). The modes disagree most on
  # sigmabias's smooth and its trade-off with boundary, which varies with the
  # same predictors. If this still splits, the fallback is sigmabias ~ 1.
  gam_rdm5r = function() {
    brms::bf(
      igc_rt(),
      igc_t2("driftone"),
      sigmabias ~ 1 + (1 | Participant),
      igc_t2("boundary"),
      igc_t2("ndt"),
      poutlier ~ 1 + (1 | Participant),
      family = cogmod::cogmod_rdm()
    )
  },

  # LBA ---------------------------------------------------------------------
  # Linear ballistic accumulator, two accumulators (cogmod_lba2). Note this is
  # lba2 and not lba1: lba1 has a single drift and no second accumulator, so it
  # cannot model the choice that `dec(Error)` carries.
  #
  # sigmazero = 1 is the scaling constraint. The LBA is only identified up to a
  # scale -- multiply every drift, the start-point range and the boundary by
  # the same constant and the likelihood is unchanged -- so one parameter has
  # to be pinned, and the convention (cogmod's own note on lba1's `sigma`) is
  # the first accumulator's drift SD. sigmaone stays free and smoothed, which
  # is what lets the two accumulators differ. Do not free sigmazero without
  # replacing the constraint with another one, or the chains will wander along
  # that ridge and Rhat will show it.
  #
  # NOT CONVERGED at full data (MullerLyer, 2026-10-01): the chains split into
  # two modes that differ in the error accumulator and in how individual
  # differences divide between drift and boundary. Max Rhat is 13.8. See
  # AGENT.md 4.9 before refitting it on any illusion.
  gam_lba = function() {
    brms::bf(
      igc_rt(),
      igc_t2("driftone"),
      igc_t2("sigmaone"),
      igc_t2("sigmabias"),
      igc_t2("boundary"),
      igc_t2("ndt"),
      sigmazero = 1,
      poutlier ~ 1 + (1 | Participant),
      family = cogmod::cogmod_lba2()
    )
  },

  # Linear DDM-5 (Study 2) --------------------------------------------------
  # DDM-5, the model retained in Study 1 (AGENT.md 6.12), with linear effects
  # instead of the t2() smooths (igc_lin above). Specified in
  # analysis/3_linearmodels.qmd: the coding of strength (absolute value, with
  # a slope per side, Effect = +-0.5, the sides meeting at zero strength), the
  # range of stimuli, the ISI centred on its study's window, and the
  # transformations of the predictors, which are the registry entry's
  # (`linear`), not the formula's. Every effect varies by participant: 6 terms
  # on each of five parameters, with their correlations.
  lin_ddm5 = function() {
    brms::bf(
      igc_lin("RT | dec(Error)"),
      igc_lin("boundary"),
      igc_lin("bias"),
      igc_lin("ndt"),
      igc_lin("sigmadrift"),
      sigmabias = 0,
      sigmandt = 0,
      poutlier ~ poly(ISI, 2, raw = TRUE) + (1 | Participant) + (1 | Study),
      family = cogmod::cogmod_ddm()
    )
  }
)


# Study 2's linear subrange ---------------------------------------------------
# The range of stimuli the linear models are fitted on, per illusion, in raw
# units: the difference from-to and the largest absolute strength. It is the
# third of the project's three ranges of the predictors (AGENT.md 6.17):
#   - the full range: every trial of every study, what the Study 1 GAMs are
#     fitted on (the bounds of restore_units() in predictions.R);
#   - the reliable subrange: the full range less the values that only the
#     Validation study showed (the extremes holding < 5% of the trials),
#     over which every prediction and figure is made (plot_ranges in
#     predictions.R; "the reliable predictor space" in the paper);
#   - this linear subrange, within the reliable one: where the effects of the
#     two predictors are taken to be close enough to linear for the linear
#     models to stand in for the GAMs.
#
# STILL TO BE REFINED. For now it is the whole reliable subrange, on all three
# illusions (2026-11-25): with the coding of 3_linearmodels.qmd (|strength|
# with a slope per side), the least linear parameter keeps R2 .87-.89 over it.
# A narrower one would trade trials for linearity, and the criterion for that
# is not decided. The boxes drawn in Study 1 (linear_regions in
# 2_modelcomparison.qmd until 2026-11-25; Müller-Lyer 0.1-0.4 x +-30 by eye)
# were for the signed coding of strength, and are candidates, not this.
#
# The registry's lin_* entries take their range from here, so changing it here
# changes what they fit. 3_linearmodels.qmd stops if it is not within
# plot_ranges, and 2_modelcomparison.qmd outlines it on its heatmaps wherever
# it is narrower than the reliable subrange. One entry per line: ./hpc reads
# model names from lines that end in "= list(".
igc_linear_ranges <- list(
  MullerLyer = list(difference = c(0.045, 0.45), strength = 30),
  VerticalHorizontal = list(difference = c(0.03, 0.30), strength = 45),
  Ebbinghaus = list(difference = c(0.05, 0.50), strength = 1.5)
)


# Registry -----------------------------------------------------------------
# One entry per family x illusion, named <family>_<illusion in lower case>.
# Shard and fit names add the illusion again, e.g.
# gam_lnr_mullerlyer_MullerLyer_<n>.rds, as they always have.
#
# Write each entry out in full rather than generating them in a loop: ./hpc
# finds models by reading their declaration lines out of this file.
#
# The three illusions are about the same size (2026-10-12 data: MullerLyer
# 324,580 rows and 2,217 participants, VerticalHorizontal 321,113 and 2,202,
# Ebbinghaus 316,487 and 2,181), so each model's production settings and
# timings carry over between them.

igc_models <- list(

  # MullerLyer ----------------------------------------------------------------

  gam_lnr_mullerlyer = list(
    illusion = "MullerLyer",
    formula = igc_formulas$gam_lnr
  ),

  gam_lnr6_mullerlyer = list(
    illusion = "MullerLyer",
    formula = igc_formulas$gam_lnr6
  ),

  gam_ddm4_mullerlyer = list(
    illusion = "MullerLyer",
    formula = igc_formulas$gam_ddm4
  ),

  gam_ddm5_mullerlyer = list(
    illusion = "MullerLyer",
    formula = igc_formulas$gam_ddm5
  ),

  # Not at full data: see igc_formulas$gam_ddm7.
  gam_ddm7_mullerlyer = list(
    illusion = "MullerLyer",
    formula = igc_formulas$gam_ddm7
  ),

  gam_rdm_mullerlyer = list(
    illusion = "MullerLyer",
    formula = igc_formulas$gam_rdm
  ),

  gam_rdm5_mullerlyer = list(
    illusion = "MullerLyer",
    formula = igc_formulas$gam_rdm5
  ),

  gam_rdm5r_mullerlyer = list(
    illusion = "MullerLyer",
    formula = igc_formulas$gam_rdm5r
  ),

  gam_lba_mullerlyer = list(
    illusion = "MullerLyer",
    formula = igc_formulas$gam_lba
  ),

  # VerticalHorizontal --------------------------------------------------------

  gam_lnr_verticalhorizontal = list(
    illusion = "VerticalHorizontal",
    formula = igc_formulas$gam_lnr
  ),

  gam_ddm4_verticalhorizontal = list(
    illusion = "VerticalHorizontal",
    formula = igc_formulas$gam_ddm4
  ),

  gam_lnr6_verticalhorizontal = list(
    illusion = "VerticalHorizontal",
    formula = igc_formulas$gam_lnr6
  ),

  gam_ddm5_verticalhorizontal = list(
    illusion = "VerticalHorizontal",
    formula = igc_formulas$gam_ddm5
  ),

  gam_rdm_verticalhorizontal = list(
    illusion = "VerticalHorizontal",
    formula = igc_formulas$gam_rdm
  ),

  gam_rdm5_verticalhorizontal = list(
    illusion = "VerticalHorizontal",
    formula = igc_formulas$gam_rdm5
  ),

  gam_rdm5r_verticalhorizontal = list(
    illusion = "VerticalHorizontal",
    formula = igc_formulas$gam_rdm5r
  ),

  gam_lba_verticalhorizontal = list(
    illusion = "VerticalHorizontal",
    formula = igc_formulas$gam_lba
  ),

  # Ebbinghaus ----------------------------------------------------------------

  gam_lnr_ebbinghaus = list(
    illusion = "Ebbinghaus",
    formula = igc_formulas$gam_lnr
  ),

  gam_lnr6_ebbinghaus = list(
    illusion = "Ebbinghaus",
    formula = igc_formulas$gam_lnr6
  ),

  gam_ddm4_ebbinghaus = list(
    illusion = "Ebbinghaus",
    formula = igc_formulas$gam_ddm4
  ),

  gam_ddm5_ebbinghaus = list(
    illusion = "Ebbinghaus",
    formula = igc_formulas$gam_ddm5
  ),

  gam_rdm_ebbinghaus = list(
    illusion = "Ebbinghaus",
    formula = igc_formulas$gam_rdm
  ),

  gam_rdm5_ebbinghaus = list(
    illusion = "Ebbinghaus",
    formula = igc_formulas$gam_rdm5
  ),

  gam_rdm5r_ebbinghaus = list(
    illusion = "Ebbinghaus",
    formula = igc_formulas$gam_rdm5r
  ),

  gam_lba_ebbinghaus = list(
    illusion = "Ebbinghaus",
    formula = igc_formulas$gam_lba
  ),

  # Study 2: linear models ----------------------------------------------------
  # `linear`: the range of stimuli fitted (difference from-to and the largest
  # absolute strength, in raw units) and the transformations of the two
  # predictors (names in igc_transforms). Untransformed first; the transformed
  # variants are compared with it once fitted (3_linearmodels.qmd).
  # The range is the illusion's linear subrange (igc_linear_ranges above), for
  # now its whole reliable subrange. Until 2026-11-25 it was written out here,
  # and before that day Müller-Lyer started at 0.05, where its Study 1
  # surfaces did. A variant fitted on another range would write its own, which
  # must lie within plot_ranges (3_linearmodels.qmd checks).

  lin_ddm5_mullerlyer = list(
    illusion = "MullerLyer",
    formula = igc_formulas$lin_ddm5,
    linear = c(igc_linear_ranges$MullerLyer,
               list(transform_difference = "identity", transform_strength = "identity"))
  ),

  # The same with the square root of the absolute strength: the transformation
  # that recovers most of what the coding without a gap costs on Müller-Lyer
  # (worst R² .917 -> .935 on the Study 1 surfaces, 3_linearmodels.qmd), and
  # the answer to the V at zero strength of the untransformed test (AGENT.md
  # 6.15). "sqrt" is on the strength only.
  lin_ddm5sqrt_mullerlyer = list(
    illusion = "MullerLyer",
    formula = igc_formulas$lin_ddm5,
    linear = c(igc_linear_ranges$MullerLyer,
               list(transform_difference = "identity", transform_strength = "sqrt"))
  ),

  # The second illusion, untransformed like Müller-Lyer's retained test
  # (2026-11-25)
  lin_ddm5_verticalhorizontal = list(
    illusion = "VerticalHorizontal",
    formula = igc_formulas$lin_ddm5,
    linear = c(igc_linear_ranges$VerticalHorizontal,
               list(transform_difference = "identity", transform_strength = "identity"))
  )
)


# Study 2's predictors --------------------------------------------------------
# The transformations a `linear` entry can name
igc_transforms <- list(identity = identity, sqrt = sqrt, log = log, log1p = log1p)

# The data of a linear model, from the rows fit_model.R has read (one illusion,
# with ISI_centre already set from every trial of each study, before any
# subsetting). Keeps the stimuli in the entry's range and adds:
#   Difference  the (transformed) difference, from -1 to 1 over the range
#   Strength    the (transformed) absolute strength, from 0 to 1 over the range
#   Effect      +0.5 conflicting, -0.5 facilitating, 0 with no illusion (the
#               formula only has it multiplied by Strength, which is then 0, so
#               its value there does not matter)
#   ISI        its position in its study's window, from -1 to 1 (in steps of
#               0.25 s from the window's centre); the raw value is kept as
#               ISI_seconds
#   Study       the study (Dataset)
# The rescaling is over the entry's nominal range, not the values in the data,
# so that a subset (a test on a few participants) has the same predictors.
igc_linear_data <- function(data, linear) {
  tol <- 1e-9
  keep <- data$Illusion_Difference >= linear$difference[1] - tol &
    data$Illusion_Difference <= linear$difference[2] + tol &
    abs(data$Illusion_Strength) <= linear$strength + tol
  data <- data[keep, ]
  td <- igc_transforms[[linear$transform_difference]]
  ts <- igc_transforms[[linear$transform_strength]]
  if (is.null(td) || is.null(ts)) stop("unknown transformation in `linear`", call. = FALSE)
  lo <- td(linear$difference[1])
  hi <- td(linear$difference[2])
  data$Difference <- 2 * (td(data$Illusion_Difference) - lo) / (hi - lo) - 1
  data$Strength <- ts(abs(data$Illusion_Strength)) / ts(linear$strength)
  data$Effect <- 0.5 * sign(data$Illusion_Strength)
  if (is.null(data$ISI_centre)) stop("ISI_centre is not set", call. = FALSE)
  data$ISI_seconds <- data$ISI
  data$ISI <- (data$ISI - data$ISI_centre) / 0.25
  data$Study <- data$Dataset
  data
}


# Look one up, with an error that says what the alternatives are rather than
# "subscript out of bounds" three hours into a job.
igc_model <- function(name) {
  if (!nzchar(name)) {
    stop("no model requested: set IGC_MODEL to one of ",
         paste(names(igc_models), collapse = ", "), call. = FALSE)
  }
  if (!name %in% names(igc_models)) {
    stop("unknown model '", name, "'. Known models: ",
         paste(names(igc_models), collapse = ", "), call. = FALSE)
  }
  c(list(name = name), igc_models[[name]])
}
