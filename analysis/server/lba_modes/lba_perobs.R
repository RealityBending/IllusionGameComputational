# Per-observation log-likelihood and posterior predictions of the two modes of
# the combined gam_lba fit, computed by hand from the draws and the standata
# design matrices (not through brms's predict path, which the local 2.23.1 gets
# wrong for 2.21.0 fits). Validated against the per-draw log-likelihood implied
# by lp__ (lba_decompose.R): the sums must agree.
suppressMessages({library(brms); library(cogmod)})
scratch <- "lba_modes"
m <- readRDS("../models/gam_lba_MullerLyer.rds")
sd <- standata(m)
a <- as.array(m$fit)
nm <- dimnames(a)[[3]]
dat <- m$data
N <- nrow(dat)
modes <- c("A", "A", "B", "A", "B", "A", "B", "A")
softplus <- function(x) ifelse(x > 30, x, log1p(exp(x)))
smooth_tag <- "t2Illusion_DifferenceZIllusion_StrengthZ_"

# the participant levels, in the order of the r_ columns
plev <- sub("^r_Participant\\[(.*),Intercept\\]$", "\\1", nm[startsWith(nm, "r_Participant[")])

eta_of <- function(x, dp, gi) {
  sfx <- if (dp == "mu") "" else paste0("_", dp)
  icpt <- if (dp == "mu") "Intercept" else paste0("Intercept_", dp)
  eta <- rep(x[icpt], N)
  if (dp != "poutlier") {
    bs <- x[paste0("bs", sfx, "_", smooth_tag, 1:3)]
    eta <- eta + as.vector(sd[[paste0("Xs", sfx)]] %*% bs)
    for (k in 1:3) {
      s <- x[nm[startsWith(nm, paste0("s", sfx, "_", smooth_tag, k, "["))]]
      eta <- eta + as.vector(sd[[paste0("Zs", sfx, "_1_", k)]] %*% s)
    }
  }
  rpre <- if (dp == "mu") "r_Participant[" else paste0("r_Participant__", dp, "[")
  r <- x[nm[startsWith(nm, rpre)]]
  eta + r[sd[[paste0("J_", gi)]]]
}
dps <- c(mu = 1, driftone = 2, sigmaone = 3, sigmabias = 4, boundary = 5, ndt = 6, poutlier = 7)
pars_of <- function(x) {
  e <- lapply(names(dps), function(dp) eta_of(x, dp, dps[[dp]]))
  names(e) <- names(dps)
  list(mu = e$mu, driftone = e$driftone, sigmaone = softplus(e$sigmaone),
       sigmabias = softplus(e$sigmabias), boundary = softplus(e$boundary),
       ndt = exp(e$ndt), poutlier = plogis(e$poutlier))
}
ll_of <- function(p) {
  dcogmod_lba2(sd$Y, driftzero = p$mu, driftone = p$driftone, sigmazero = 1,
               sigmaone = p$sigmaone, sigmabias = p$sigmabias, boundary = p$boundary,
               ndt = p$ndt, response = sd$dec, poutlier = p$poutlier, log = TRUE)
}

# Validation: one draw per chain against the decomposition
dec_ll <- function(x) {
  # same arithmetic as lba_decompose.R, for one draw
  z <- 0; jac <- 0; zs <- 0
  for (dp in c("", "driftone", "sigmaone", "sigmabias", "boundary", "ndt", "poutlier")) {
    sdn <- if (dp == "") "sd_Participant__Intercept" else paste0("sd_Participant__", dp, "_Intercept")
    rpre <- if (dp == "") "r_Participant[" else paste0("r_Participant__", dp, "[")
    z <- z + sum(-0.5 * (x[nm[startsWith(nm, rpre)]] / x[sdn])^2); jac <- jac + log(x[sdn])
    if (dp != "poutlier") {
      tag <- if (dp == "") "" else paste0("_", dp)
      for (k in 1:3) {
        sdsn <- paste0("sds", tag, "_", smooth_tag, k)
        zs <- zs + sum(-0.5 * (x[nm[startsWith(nm, paste0("s", tag, "_", smooth_tag, k, "["))]] / x[sdsn])^2)
        jac <- jac + log(x[sdsn])
      }
    }
  }
  unname(x["lp__"] - x["lprior"] - z - zs - jac)
}
cat("validation (hand-computed sum vs lp__-implied, up to a constant):\n")
val <- t(sapply(1:8, function(ch) {
  x <- a[500, ch, ]
  c(chain = ch, hand = sum(ll_of(pars_of(x))), implied = dec_ll(x))
}))
val <- cbind(val, diff = val[, "hand"] - val[, "implied"])
print(val, digits = 8)

# Per-observation log-likelihood, averaged over 20 draws per chain, by mode;
# posterior-mean parameters by mode; one simulated data set per draw.
set.seed(1)
idx <- round(seq(25, 500, length.out = 20))
acc <- list(A = list(ll = 0, n = 0, p = NULL), B = list(ll = 0, n = 0, p = NULL))
sims <- list()
for (ch in 1:8) {
  md <- modes[ch]
  for (i in idx) {
    p <- pars_of(a[i, ch, ])
    acc[[md]]$ll <- acc[[md]]$ll + ll_of(p)
    acc[[md]]$n <- acc[[md]]$n + 1
    acc[[md]]$p <- if (is.null(acc[[md]]$p)) p else Map(`+`, acc[[md]]$p, p)
    if (i %in% idx[c(5, 10, 15, 20)]) {
      s <- rcogmod_lba2(N, driftzero = p$mu, driftone = p$driftone, sigmazero = 1,
                        sigmaone = p$sigmaone, sigmabias = p$sigmabias, boundary = p$boundary,
                        ndt = p$ndt, poutlier = p$poutlier)
      sims[[length(sims) + 1]] <- data.frame(mode = md, chain = ch, draw = i, row = seq_len(N),
                                             rt = s$rt, response = s$response)
    }
  }
  cat("chain", ch, "done\n")
}
res <- data.frame(dat, llA = acc$A$ll / acc$A$n, llB = acc$B$ll / acc$B$n)
for (k in names(acc$A$p)) {
  res[[paste0(k, "_A")]] <- acc$A$p[[k]] / acc$A$n
  res[[paste0(k, "_B")]] <- acc$B$p[[k]] / acc$B$n
}
saveRDS(list(val = val, res = res, sims = do.call(rbind, sims)), file.path(scratch, "lba_perobs.rds"))
cat("total ll A:", sum(res$llA), " B:", sum(res$llB), " B-A:", sum(res$llB - res$llA), "\n")
