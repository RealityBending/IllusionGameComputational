# =========================================================================
# Priors against posteriors, for 1_models.qmd
# =========================================================================
# For every population-level parameter of a fit (intercepts, smooth slopes
# bs_*, smooth SDs sds_*, participant SDs sd_*): the prior it was fitted with,
# as brms resolved it, and where its posterior sits in that prior. Two things
# are flagged:
#   - prior-data conflict: the posterior median beyond the prior's 1st or 99th
#     percentile, i.e. the prior is pulling against the data;
#   - weak data: the posterior barely narrower than the prior (contraction
#     below 0.1), i.e. the prior, not the data, sets the answer;
#   - flat prior: no prior at all (brms's improper default), so nothing to
#     compare with.
#
# The priors are resolved from prior_summary(m, all = TRUE) with brms's own
# inheritance (coefficient -> group -> class within a dpar), then checked by
# summing their log densities at every draw and comparing with the `lprior`
# Stan saved. A mismatch means a prior was matched to the wrong parameter, so
# the table is not to be trusted (meta$lprior_maxdiff, after removing a
# constant offset, meta$lprior_offset).
#
# Shared by prior_check.R (the cluster job, `./hpc priors`), which writes one
# file per fit to models/priors/<name>.rds; the qmd only reads those.

igc_priors_version <- "1"

# "student_t(3, 0, 2.5)" -> log density, CDF and quantile function, truncated
# below at lb as brms does for the SDs (its lprior subtracts the lccdf at lb).
igc_prior_dist <- function(prior, lb = -Inf) {
  s <- gsub("\\s", "", prior)
  fam <- sub("\\(.*$", "", s)
  a <- as.numeric(strsplit(sub("^[^(]*\\((.*)\\)$", "\\1", s), ",")[[1]])
  f <- switch(fam,
    normal = list(d = \(x) dnorm(x, a[1], a[2], log = TRUE), p = \(x) pnorm(x, a[1], a[2]),
                  q = \(u) qnorm(u, a[1], a[2])),
    student_t = list(d = \(x) dt((x - a[2]) / a[3], a[1], log = TRUE) - log(a[3]),
                     p = \(x) pt((x - a[2]) / a[3], a[1]), q = \(u) a[2] + a[3] * qt(u, a[1])),
    cauchy = list(d = \(x) dcauchy(x, a[1], a[2], log = TRUE), p = \(x) pcauchy(x, a[1], a[2]),
                  q = \(u) qcauchy(u, a[1], a[2])),
    exponential = list(d = \(x) dexp(x, a[1], log = TRUE), p = \(x) pexp(x, a[1]),
                       q = \(u) qexp(u, a[1])),
    gamma = list(d = \(x) dgamma(x, a[1], a[2], log = TRUE), p = \(x) pgamma(x, a[1], a[2]),
                 q = \(u) qgamma(u, a[1], a[2])),
    logistic = list(d = \(x) dlogis(x, a[1], a[2], log = TRUE), p = \(x) plogis(x, a[1], a[2]),
                    q = \(u) qlogis(u, a[1], a[2])),
    lognormal = list(d = \(x) dlnorm(x, a[1], a[2], log = TRUE), p = \(x) plnorm(x, a[1], a[2]),
                     q = \(u) qlnorm(u, a[1], a[2])),
    stop("priors.R does not know the prior '", prior, "'")
  )
  if (is.finite(lb)) {
    F0 <- f$p(lb)
    g <- f
    f$d <- \(x) g$d(x) - log1p(-F0)
    f$p <- \(x) (g$p(x) - F0) / (1 - F0)
    f$q <- \(u) g$q(F0 + u * (1 - F0))
  }
  f
}

# The brms class, dpar, group and coefficient of a parameter name, or NULL for
# what is not population-level (r_*, s_*, lp__, ...) or duplicates another
# (b_Intercept is the same draws as Intercept, whose prior it is).
igc_parse_par <- function(par, dpars) {
  dpars <- dpars[order(-nchar(dpars))]
  split_dpar <- function(rest) {
    for (d in dpars) if (startsWith(rest, paste0(d, "_"))) return(c(d, substring(rest, nchar(d) + 2)))
    c("", rest)
  }
  if (par == "Intercept") return(list(class = "Intercept", dpar = "", group = "", coef = ""))
  if (startsWith(par, "Intercept_")) {
    return(list(class = "Intercept", dpar = sub("^Intercept_", "", par), group = "", coef = ""))
  }
  if (startsWith(par, "bs_")) {
    x <- split_dpar(sub("^bs_", "", par))
    return(list(class = "b", dpar = x[1], group = "", coef = x[2]))
  }
  if (startsWith(par, "sds_")) {
    x <- split_dpar(sub("^sds_", "", par))
    return(list(class = "sds", dpar = x[1], group = "", coef = x[2]))
  }
  if (grepl("^sd_.+__", par)) {
    group <- sub("^sd_(.+?)__.*$", "\\1", par)
    x <- split_dpar(sub("^sd_.+?__", "", par))
    return(list(class = "sd", dpar = x[1], group = group, coef = x[2]))
  }
  if (startsWith(par, "b_")) {
    x <- split_dpar(sub("^b_", "", par))
    if (x[2] == "Intercept") return(NULL)
    return(list(class = "b", dpar = x[1], group = "", coef = x[2]))
  }
  NULL
}

# The prior brms used for one parameter: the most specific non-empty row of
# prior_summary(all = TRUE) within its class and dpar. Rows that exist but are
# all empty mean brms's flat (improper) default, which adds nothing to lprior:
# returned as prior "" (e.g. the boundary intercept of the cogmod <= 0.3.3
# DDMs). No row at all is a parsing error.
igc_resolve_prior <- function(ps, p) {
  cand <- ps[ps$class == p$class & ps$dpar == p$dpar, ]
  if (nrow(cand) == 0) return(NULL)
  for (k in list(c(p$group, p$coef), c(p$group, ""), c("", ""))) {
    r <- cand[cand$group == k[1] & cand$coef == k[2] & nzchar(cand$prior), ]
    if (nrow(r)) {
      lb <- suppressWarnings(as.numeric(cand$lb[nzchar(cand$lb)][1]))
      return(list(prior = r$prior[1], lb = if (is.na(lb)) -Inf else lb, source = r$source[1]))
    }
  }
  list(prior = "", lb = -Inf, source = cand$source[1])
}

igc_prior_check <- function(m, name, fit_file = NULL) {
  ps <- as.data.frame(brms::prior_summary(m, all = TRUE))
  for (k in c("class", "coef", "group", "dpar", "lb", "prior", "source")) ps[[k]][is.na(ps[[k]])] <- ""
  dpars <- unique(ps$dpar[nzchar(ps$dpar)])

  # a plain iterations x chains x variables array: the cluster's posterior has
  # no extract_variable_array()
  draws <- unclass(posterior::as_draws_array(m))
  vars <- dimnames(draws)[[3]]
  parsed <- lapply(vars, igc_parse_par, dpars = dpars)
  keep <- !vapply(parsed, is.null, logical(1))
  vars <- vars[keep]; parsed <- parsed[keep]

  # a fixed grid of quantiles stands in for prior draws (prior SD, truncated
  # or not, without randomness)
  u <- (seq_len(20000) - 0.5) / 20000
  lp_sum <- 0
  rows <- vector("list", length(vars))
  for (i in seq_along(vars)) {
    p <- parsed[[i]]
    pr <- igc_resolve_prior(ps, p)
    if (is.null(pr)) stop(name, ": no prior found for ", vars[i])
    flat <- !nzchar(pr$prior)
    # a flat prior: no density, CDF or SD, so NA for everything prior-relative
    f <- if (flat) list(d = \(x) 0, p = \(x) NA_real_, q = \(u) rep(NA_real_, length(u)))
         else igc_prior_dist(pr$prior, pr$lb)
    x <- draws[, , vars[i]]
    xv <- as.vector(x)
    lp_sum <- lp_sum + f$d(x)
    prior_sd <- sd(f$q(u))
    med <- median(xv)
    rows[[i]] <- data.frame(
      Parameter = vars[i], Class = p$class, Dpar = if (nzchar(p$dpar)) p$dpar else "mu",
      Coef = p$coef, Prior = if (flat) "flat (improper)" else pr$prior, Prior_LB = pr$lb,
      Prior_Source = pr$source,
      Prior_Median = f$q(0.5), Prior_Q01 = f$q(0.01), Prior_Q99 = f$q(0.99), Prior_SD = prior_sd,
      Post_Median = med, Post_Q05 = unname(quantile(xv, 0.05)), Post_Q95 = unname(quantile(xv, 0.95)),
      Post_SD = sd(xv), Rhat = posterior::rhat(x), ESS = posterior::ess_bulk(x),
      PQ_Median = f$p(med), PQ_Q05 = f$p(unname(quantile(xv, 0.05))), PQ_Q95 = f$p(unname(quantile(xv, 0.95))),
      Contraction = 1 - sd(xv)^2 / prior_sd^2
    )
  }
  tab <- do.call(rbind, rows)
  tab$Flag <- ifelse(tab$Prior == "flat (improper)", "flat prior",
              ifelse(tab$PQ_Median < 0.01 | tab$PQ_Median > 0.99, "conflict",
              ifelse(tab$PQ_Median < 0.05 | tab$PQ_Median > 0.95, "tail",
              ifelse(tab$Contraction < 0.1, "weak data", ""))))

  # Did we resolve the priors Stan actually used? Our sum of log densities
  # should reproduce lprior at every draw, up to the 6 significant digits
  # CmdStan writes, and up to a constant: brms 2.21 subtracts a truncated
  # prior's normalising lccdf once per vector, not once per element, so a
  # half-Student-t on the 3 sds of a t2 smooth is off by 2 log 2 (the 0.3.3
  # fits). A constant cannot come from a wrong prior; a deviation that varies
  # over draws does.
  lprior <- draws[, , "lprior"]
  diff <- as.vector(lprior) - as.vector(lp_sum)
  offset <- median(diff)

  list(
    meta = list(
      version = igc_priors_version, name = name, date = format(Sys.time()),
      cogmod = if (!is.null(m$cogmod)) m$cogmod$version else "<= 0.3.3",
      brms = as.character(m$version$brms), family = m$family$name,
      formula = paste(capture.output(print(m$formula)), collapse = "\n"),
      ndraws = prod(dim(draws)[1:2]), nchains = dim(draws)[2],
      fit_size = if (!is.null(fit_file)) file.size(fit_file) else NA,
      lprior_maxdiff = max(abs(diff - offset)), lprior_offset = offset
    ),
    priors = ps,
    table = tab
  )
}
