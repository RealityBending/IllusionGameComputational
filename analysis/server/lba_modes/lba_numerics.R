# Is cogmod's LBA2 Stan density accurate, and smooth, where gam_lba lives?
# Reference: direct numerical integration over the start point, plus rtdists.
suppressMessages(library(cogmod))
sp <- function(x) log1p(exp(x))

# Reference defective density of the race (decision time t, winner k), each drift
# a Normal truncated at zero. Winner: f(t) = (1/A) int_0^A phiTN((b-z)/t) (b-z)/t^2 dz.
# Loser: S(t) = (1/A) int_0^A P(v < (b-z)/t | v > 0) dz.
ref_ldens <- function(t, vw, sw, vl, sl, A, B) {
  b <- A + B
  lq_w <- pnorm(vw / sw, log.p = TRUE)
  fw <- integrate(function(z) {
    v <- (b - z) / t
    exp(dnorm(v, vw, sw, log = TRUE) - lq_w) * (b - z) / t^2
  }, 0, A, rel.tol = 1e-12, abs.tol = 0)$value / A
  # survival of the loser: P(v < x | v > 0) = (Phi((x-vl)/sl) - Phi(-vl/sl)) / Phi(vl/sl)
  Sl <- integrate(function(z) {
    x <- (b - z) / t
    # P(0 < v < x) taken from the upper tail to keep digits: Phibar(-vl/sl) - Phibar((x-vl)/sl)
    p <- pnorm(-vl / sl, lower.tail = FALSE) - pnorm((x - vl) / sl, lower.tail = FALSE)
    p / pnorm(vl / sl)
  }, 0, A, rel.tol = 1e-12, abs.tol = 0)$value / A
  log(fw) + log(Sl)
}

modes <- list(
  A = list(mu = 2.537, driftone = -4.87, sigmaone = sp(1.996), sigmabias = sp(-0.040), boundary = sp(0.381)),
  B = list(mu = 2.665, driftone = -7.00, sigmaone = sp(2.60), sigmabias = sp(-0.46), boundary = sp(0.58))
)
ts <- c(0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2, 3.5)

cat("Compiling cogmod's Stan lpdf...\n")
lpdf <- cogmod_lba2_lpdf_expose()

rows <- list()
for (mn in names(modes)) for (dec in 0:1) for (t in ts) {
  p <- modes[[mn]]
  # Stan: full mixture, so set ndt tiny and poutlier ~0 to isolate the race
  stan <- lpdf(Y = t + 1e-3, mu = p$mu, driftone = p$driftone, sigmazero = 1, sigmaone = p$sigmaone,
               sigmabias = p$sigmabias, boundary = p$boundary, ndt = 1e-3, poutlier = 1e-300, dec = dec)
  r <- cogmod:::.lba2_ldens(t, dec, list(mu = p$mu, driftone = p$driftone, sigmazero = 1,
         sigmaone = p$sigmaone, sigmabias = p$sigmabias, boundary = p$boundary))
  if (dec == 0) ref <- ref_ldens(t, p$mu, 1, p$driftone, p$sigmaone, p$sigmabias, p$boundary)
  else ref <- ref_ldens(t, p$driftone, p$sigmaone, p$mu, 1, p$sigmabias, p$boundary)
  rtd <- log(rtdists::n1PDF(t, A = p$sigmabias, b = p$sigmabias + p$boundary, t0 = 0,
           mean_v = if (dec == 0) c(p$mu, p$driftone) else c(p$driftone, p$mu),
           sd_v = if (dec == 0) c(1, p$sigmaone) else c(p$sigmaone, 1), silent = TRUE))
  rows[[length(rows) + 1]] <- data.frame(mode = mn, dec = dec, t = t, stan = stan, R = r, ref = ref, rtdists = rtd)
}
res <- do.call(rbind, rows)
res$stan_minus_ref <- res$stan - res$ref
res$rtd_minus_ref <- res$rtdists - res$ref
print(res, digits = 6)
cat("\nmax |stan - ref| (finite rows):", max(abs(res$stan_minus_ref[is.finite(res$stan_minus_ref)])), "\n")

# Gradient smoothness: central differences at two step sizes along driftone and
# sigmaone for an error trial and a correct trial, over a sweep of driftone that
# covers the path a chain takes from the init (+3) down to mode B (-7) and past.
if (Sys.getenv("SKIP_GRAD") == "1") quit("no")
cat("\nFinite-difference gradients along driftone (h = 1e-4 vs 1e-6):\n")
g <- function(fun, x, h) (fun(x + h) - fun(x - h)) / (2 * h)
sweep <- c(3, 1, 0.01, -0.01, -1, -3, -5, -7, -10, -15, -20)
gr <- list()
for (dec in 0:1) for (t in c(0.15, 0.4, 1.0)) for (v1 in sweep) {
  s1 <- 2.67 * max(1, abs(v1) / 7)  # roughly along the |v| / s^2 ray
  fun <- function(v) lpdf(Y = t + 1e-3, mu = 2.6, driftone = v, sigmazero = 1, sigmaone = s1,
                          sigmabias = 0.5, boundary = 1, ndt = 1e-3, poutlier = 1e-300, dec = dec)
  gr[[length(gr) + 1]] <- data.frame(dec = dec, t = t, driftone = v1, sigmaone = s1, lp = fun(v1),
                                     g4 = g(fun, v1, 1e-4), g6 = g(fun, v1, 1e-6))
}
gr <- do.call(rbind, gr)
gr$rel_diff <- abs(gr$g4 - gr$g6) / pmax(abs(gr$g4), 1e-8)
print(gr, digits = 5)
saveRDS(list(res = res, gr = gr), "lba_modes/lba_numerics.rds")
