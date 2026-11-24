# Are the two modes of the 2026-10-13 refit of gam_lnr6_mullerlyer one ridge,
# or separated by low density? (AGENT.md 4.10 B, "Another sampler would not
# join them".) Reads only the draws, through posterior, so it is safe with the
# laptop's brms (no restructure()). Usage:
#   Rscript lnr6_ml_ridge.R path/to/gam_lnr6_mullerlyer_MullerLyer.rds
# The mode membership of the chains (A, B) is that fit's: re-derive it from the
# per-chain means of b_Intercept before reusing this on another fit.
suppressPackageStartupMessages(library(posterior))
m <- readRDS(commandArgs(trailingOnly = TRUE)[1])
d <- as_draws_array(m$fit)
vars <- variables(d)
pop <- vars[grepl("^(b_|bs_|sds_|sd_)", vars)]
x <- as_draws_array(subset_draws(d, variable = pop))
A <- c(1, 5, 6, 7); B <- c(2, 3, 4, 8)
arr <- unclass(x)  # iter x chain x var
# standardise every parameter by its pooled within-chain SD
wsd <- apply(arr, 3, \(v) sqrt(mean(apply(v, 2, var))))
z <- sweep(arr, 3, wsd, "/")
mA <- apply(z[, A, ], 3, mean); mB <- apply(z[, B, ], 3, mean)
dir <- (mB - mA); gap <- sqrt(sum(dir^2)); u <- dir / gap
cat("distance between mode centres, in pooled within-chain SDs:", round(gap, 1), "\n")
proj <- apply(z, c(1, 2), \(v) sum(v * u))
cat("projection on the between-mode axis, per chain (mean / sd):\n")
print(round(rbind(mean = colMeans(proj), sd = apply(proj, 2, sd)), 2))
cat("largest draw in mode A:", round(max(proj[, A]), 2), " smallest in mode B:", round(min(proj[, B]), 2), "\n")
# within-mode covariance: does its leading direction point along the between-mode axis?
for (nm in c("A", "B")) {
  ch <- get(nm)
  w <- do.call(rbind, lapply(ch, \(c) scale(z[, c, ], scale = FALSE)))
  e <- eigen(cov(w), symmetric = TRUE)
  cosv <- abs(as.vector(t(e$vectors[, 1:5]) %*% u))
  cat(sprintf("mode %s: within-mode eigenvalues (top 5) %s; |cos| of those axes with the between-mode axis %s; variance along the between-mode axis %.2f\n",
              nm, paste(round(e$values[1:5], 1), collapse = " / "), paste(round(cosv, 2), collapse = " / "),
              as.numeric(t(u) %*% cov(w) %*% u)))
}
cat("\nwithin-chain cor(b_Intercept, b_sigmabias_Intercept):",
    round(sapply(1:8, \(c) cor(arr[, c, "b_Intercept"], arr[, c, "b_sigmabias_Intercept"])), 2), "\n")
cat("top loadings of the between-mode axis:\n")
print(round(head(sort(abs(u), decreasing = TRUE), 10), 2))
