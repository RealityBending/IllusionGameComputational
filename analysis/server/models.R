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
  }
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

  gam_lba_ebbinghaus = list(
    illusion = "Ebbinghaus",
    formula = igc_formulas$gam_lba
  )
)


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
