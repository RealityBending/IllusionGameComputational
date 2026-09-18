# =========================================================================
# Model registry -- the one place a model is defined
# =========================================================================
# Sourced by fit_model.R (which fits exactly one of these per job) and by
# combine_model.R (which merges that model's shards). `./hpc` reads the model
# *names* straight out of this file, so keep the declaration lines in the form
#
#     <name> = list(
#
# at two-space indent, one per model. `./hpc models` prints what it found.
#
# Each entry:
#   illusion  the Illusion_Type the model is fitted to
#   formula   a function returning the brms bf(). A function rather than the
#             object so that sourcing this file costs nothing and only the
#             model actually being fitted is built.

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


igc_models <- list(

  # LNR ---------------------------------------------------------------------
  # sigmabias = 0 is the plain LNR. Omitting it does not error: it silently
  # becomes a freely estimated start-point range. Watch for this class of bug
  # whenever cogmod adds a distributional parameter.
  gam_lnr = list(
    illusion = "MullerLyer",
    formula = function() {
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
    }
  ),

  # LNR-6 -------------------------------------------------------------------
  # The LNR above plus a freely estimated sigmabias -- the between-trial
  # start-point range -- smoothed like every other parameter. Six smoothed
  # dpars against the plain LNR's five, hence the name.
  #
  # Expect it to cost more per gradient than gam_lnr, not less: with
  # sigmabias > 0 cogmod takes the erfc-based ldiff_Phi path for the tails
  # rather than the single-tail shortcut, which cogmod 0.3.3 measured at about
  # 15% dearer per gradient (and 20% cheaper for the sigmabias = 0 case).
  gam_lnr6 = list(
    illusion = "MullerLyer",
    formula = function() {
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
    }
  ),

  # DDM-4 -------------------------------------------------------------------
  # The four-parameter DDM: drift, boundary, bias, ndt. The three between-trial
  # variabilities are fixed at 0.
  gam_ddm4 = list(
    illusion = "MullerLyer",
    formula = function() {
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
    }
  ),

  # DDM-5 -------------------------------------------------------------------
  # DDM-4 plus between-trial drift variability, itself smoothed. One more 2-D
  # smooth than the others, so expect it to be the straggler; give it its own
  # --time if it is ever submitted.
  gam_ddm5 = list(
    illusion = "MullerLyer",
    formula = function() {
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
    }
  ),

  # DDM-7 -------------------------------------------------------------------
  # The full DDM: nothing fixed. Every one of cogmod_ddm()'s seven parameters
  # gets its own 2-D smooth, including all three between-trial variabilities.
  #
  # DO NOT SUBMIT THIS AS IT STANDS (measured 2026-09-18, AGENT.md 4.7). At 30
  # participants it did not reach iteration 100 in 50 minutes, where every
  # other model in the registry finished all 400 in 6-11 -- upwards of 25x the
  # per-iteration cost on 1/74th of the production data, which scales to about
  # a fortnight per chain at full data, past long's 8-day ceiling.
  #
  # That is geometry, not compute, so a longer --time will not rescue it.
  # sigmadrift, sigmabias and sigmandt are identified through the shape of the
  # RT distribution rather than its location, weakly so even with flat
  # predictors, and here each carries 25 tensor coefficients plus a participant
  # intercept. The entry stays because it is where the investigation into
  # reparameterising it starts.
  gam_ddm7 = list(
    illusion = "MullerLyer",
    formula = function() {
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
    }
  ),

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
  gam_rdm = list(
    illusion = "MullerLyer",
    formula = function() {
      brms::bf(
        igc_rt(),
        igc_t2("driftone"),
        igc_t2("boundary"),
        igc_t2("ndt"),
        sigmabias = 0,
        poutlier ~ 1 + (1 | Participant),
        family = cogmod::cogmod_rdm()
      )
    }
  ),

  # RDM-5 -------------------------------------------------------------------
  # The RDM with the start-point range freed and smoothed like everything else,
  # so all five of cogmod_rdm()'s parameters are estimated. Stands to gam_rdm
  # as gam_lnr6 does to gam_lnr, and the waic comparison against gam_rdm is
  # what says whether the start-point range is identified on this data.
  #
  # Unlike the LBA, the RDM needs no scaling constraint -- its diffusion
  # coefficient is fixed internally -- so freeing sigmabias here does not open
  # the ridge that freeing sigmazero would in gam_lba.
  gam_rdm5 = list(
    illusion = "MullerLyer",
    formula = function() {
      brms::bf(
        igc_rt(),
        igc_t2("driftone"),
        igc_t2("sigmabias"),
        igc_t2("boundary"),
        igc_t2("ndt"),
        poutlier ~ 1 + (1 | Participant),
        family = cogmod::cogmod_rdm()
      )
    }
  ),

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
  gam_lba = list(
    illusion = "MullerLyer",
    formula = function() {
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
