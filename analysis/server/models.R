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
