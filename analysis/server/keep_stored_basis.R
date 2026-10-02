# =========================================================================
# Make brms >= 2.23 predict a fit made with brms < 2.23 from the fit's own
# smooth bases (AGENT.md 3.9)
# =========================================================================
# brms 2.23's restructure() throws away the smooth bases stored in an older
# fit and rebuilds them with the local mgcv/LAPACK. On the laptop, the rebuilt
# t2() bases have sign-flipped columns, so every prediction and log_lik() of a
# cluster fit is wrong, silently. This puts the stored smooths back and keeps
# the rebuilt group levels (which 2.23 needs under their new name,
# group_levels). Verified identical to brms 2.21.0 on gam_lnr6, gam_rdm and
# gam_ddm4 on 2026-10-02. It covers t2()/s() smooths and grouping factors,
# which is all these fits contain -- not gp() terms.
#
#   source("server/keep_stored_basis.R")
#   m <- keep_stored_basis(readRDS("models/gam_lnr6_MullerLyer.rds"))
#
# The cluster (./hpc predict) and a local brms 2.21.0 remain the first
# choices; run_predictions() refuses a patched fit on purpose.
keep_stored_basis <- function(m) {
  fitted_with <- m$version$brms
  if (fitted_with >= "2.23.0" || utils::packageVersion("brms") < "2.23.0") return(m)
  if (!is.null(m$version$restructure) && m$version$restructure >= "2.23.0") {
    stop("This fit was already restructured by brms >= 2.23, so its stored bases are gone. ",
         "Read it again from its file.", call. = FALSE)
  }
  if (any(lengths(lapply(m$basis$dpars, `[[`, "gp")) > 0)) {
    stop("The fit has gp() terms, whose basis format changed in brms 2.23; not covered.", call. = FALSE)
  }
  stored <- m$basis
  m <- brms::restructure(m)
  # The rebuilt levels must be the stored ones under their new name, or the
  # participants of a newdata would be matched to the wrong r_ columns
  if (!identical(unname(m$basis$group_levels), unname(stored$levels))) {
    stop("The rebuilt group levels differ from the stored ones.", call. = FALSE)
  }
  m$basis$dpars <- stored$dpars
  m
}
