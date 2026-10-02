# Is there a barrier between the two modes? Interpolate linearly, on the
# unconstrained scale the sampler moves on (non-centred z = r / sd and
# zs = s / sds, log sd, log sds, raw intercepts and bs), from one draw of a
# mode-A chain to one draw of a mode-B chain, and evaluate the log-likelihood
# plus the std_normal terms on z and zs along the way. lprior (which differs by
# ~50 between the modes) is left out.
suppressMessages({library(brms); library(cogmod)})
scratch <- "lba_modes"
m <- readRDS("../models/gam_lba_MullerLyer.rds")
sd <- standata(m)
a <- as.array(m$fit)
nm <- dimnames(a)[[3]]
N <- sd$N
softplus <- function(x) ifelse(x > 30, x, log1p(exp(x)))
tag <- "t2Illusion_DifferenceZIllusion_StrengthZ_"
dps <- c("", "driftone", "sigmaone", "sigmabias", "boundary", "ndt", "poutlier")

# to and from the unconstrained (sampler) parametrisation
to_u <- function(x) {
  u <- x
  for (dp in dps) {
    sdn <- if (dp == "") "sd_Participant__Intercept" else paste0("sd_Participant__", dp, "_Intercept")
    rp <- if (dp == "") "r_Participant[" else paste0("r_Participant__", dp, "[")
    rc <- nm[startsWith(nm, rp)]
    u[rc] <- x[rc] / x[sdn]; u[sdn] <- log(x[sdn])
    if (dp != "poutlier") {
      t <- if (dp == "") "" else paste0("_", dp)
      for (k in 1:3) {
        sdsn <- paste0("sds", t, "_", tag, k)
        sc <- nm[startsWith(nm, paste0("s", t, "_", tag, k, "["))]
        u[sc] <- x[sc] / x[sdsn]; u[sdsn] <- log(x[sdsn])
      }
    }
  }
  u
}
from_u <- function(u) {
  x <- u; zterm <- 0
  for (dp in dps) {
    sdn <- if (dp == "") "sd_Participant__Intercept" else paste0("sd_Participant__", dp, "_Intercept")
    rp <- if (dp == "") "r_Participant[" else paste0("r_Participant__", dp, "[")
    rc <- nm[startsWith(nm, rp)]
    x[sdn] <- exp(u[sdn]); x[rc] <- u[rc] * x[sdn]; zterm <- zterm + sum(-0.5 * u[rc]^2)
    if (dp != "poutlier") {
      t <- if (dp == "") "" else paste0("_", dp)
      for (k in 1:3) {
        sdsn <- paste0("sds", t, "_", tag, k)
        sc <- nm[startsWith(nm, paste0("s", t, "_", tag, k, "["))]
        x[sdsn] <- exp(u[sdsn]); x[sc] <- u[sc] * x[sdsn]; zterm <- zterm + sum(-0.5 * u[sc]^2)
      }
    }
  }
  attr(x, "zterm") <- zterm
  x
}
eta_of <- function(x, dp, gi) {
  sfx <- if (dp == "mu") "" else paste0("_", dp)
  icpt <- if (dp == "mu") "Intercept" else paste0("Intercept_", dp)
  eta <- rep(x[icpt], N)
  if (dp != "poutlier") {
    eta <- eta + as.vector(sd[[paste0("Xs", sfx)]] %*% x[paste0("bs", sfx, "_", tag, 1:3)])
    for (k in 1:3) eta <- eta + as.vector(sd[[paste0("Zs", sfx, "_1_", k)]] %*% x[nm[startsWith(nm, paste0("s", sfx, "_", tag, k, "["))]])
  }
  rp <- if (dp == "mu") "r_Participant[" else paste0("r_Participant__", dp, "[")
  eta + x[nm[startsWith(nm, rp)]][sd[[paste0("J_", gi)]]]
}
ll_of <- function(x) {
  g <- c(mu = 1, driftone = 2, sigmaone = 3, sigmabias = 4, boundary = 5, ndt = 6, poutlier = 7)
  e <- lapply(names(g), function(dp) eta_of(x, dp, g[[dp]])); names(e) <- names(g)
  sum(dcogmod_lba2(sd$Y, driftzero = e$mu, driftone = e$driftone, sigmazero = 1,
                   sigmaone = softplus(e$sigmaone), sigmabias = softplus(e$sigmabias),
                   boundary = softplus(e$boundary), ndt = exp(e$ndt), response = sd$dec,
                   poutlier = plogis(e$poutlier), log = TRUE))
}

pairs <- list(c(1, 5), c(6, 3))   # (mode-A chain, mode-B chain)
out <- list()
for (pr in pairs) {
  uA <- to_u(a[500, pr[1], ]); uB <- to_u(a[500, pr[2], ])
  for (w in c(-0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)) {
    x <- from_u((1 - w) * uA + w * uB)
    ll <- ll_of(x)
    out[[length(out) + 1]] <- data.frame(from = pr[1], to = pr[2], w = w, loglik = ll,
                                         zterm = attr(x, "zterm"), target = ll + attr(x, "zterm"),
                                         driftone = x["Intercept_driftone"], sigmaone = x["Intercept_sigmaone"])
    cat(sprintf("%d->%d  w=%.1f  loglik=%.1f  target=%.1f\n", pr[1], pr[2], w, ll, ll + attr(x, "zterm")))
  }
}
saveRDS(do.call(rbind, out), file.path(scratch, "lba_path.rds"))
