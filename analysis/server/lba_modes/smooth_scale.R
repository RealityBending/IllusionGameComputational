# What does a smooth-SD prior mean on the link scale, for this t2 basis?
# For coefficients zs ~ N(0, 1), the smooth term k contributes sds_k * (Zs_k %*% zs)
# to the linear predictor, whose variance at row i is sds_k^2 * sum_j Zs_k[i, j]^2.
# So the RMS over rows of sqrt(rowSums(Zs_k^2)) converts sds_k into "link units of
# wiggle". The null space (Xs, the unpenalised linear part) gets class-b priors.
# Then: the posterior smooths of gam_lba in those units, by mode.
options(width = 200)
suppressMessages(library(brms))
m <- readRDS("../models/gam_lba_MullerLyer.rds")
sd <- standata(m)
a <- as.array(m$fit)
nm <- dimnames(a)[[3]]
tag <- "t2Illusion_DifferenceZIllusion_StrengthZ_"
modeA <- c(1, 2, 4, 6, 8); modeB <- c(3, 5, 7)

cat("== basis scale: RMS link-unit SD per unit sds (per t2 penalty), and SD of the Xs columns\n")
k_scale <- sapply(1:3, function(k) sqrt(mean(rowSums(sd[[paste0("Zs_1_", k)]]^2))))
print(round(c(setNames(k_scale, paste0("Zs_", 1:3)), setNames(apply(sd$Xs, 2, sd), paste0("Xs_", 1:3))), 3))
cat("exponential(1) median 0.69 -> link-unit SD per penalty:", round(0.693 * k_scale, 2), "\n")
cat("student_t(3, 0, 2.5) half median 1.9 -> link-unit SD per penalty:", round(1.9 * k_scale, 2), "\n")

cat("\n== posterior: SD over rows of each smooth component, and range of the whole smooth (link units)\n")
rows <- list()
for (dp in c("mu", "driftone", "sigmaone", "sigmabias", "boundary", "ndt")) {
  sfx <- if (dp == "mu") "" else paste0("_", dp)
  for (md in c("A", "B")) {
    ch <- if (md == "A") modeA else modeB
    post <- function(p) colMeans(matrix(a[, ch, p], ncol = length(p)))
    tot <- as.vector(sd[[paste0("Xs", sfx)]] %*% post(paste0("bs", sfx, "_", tag, 1:3)))
    lin_sd <- sd(tot)
    comp <- sapply(1:3, function(k) {
      f <- as.vector(sd[[paste0("Zs", sfx, "_1_", k)]] %*% post(nm[startsWith(nm, paste0("s", sfx, "_", tag, k, "["))]))
      tot <<- tot + f
      sd(f)
    })
    sds <- post(paste0("sds", sfx, "_", tag, 1:3))
    rows[[length(rows) + 1]] <- data.frame(dpar = dp, mode = md, linear_sd = lin_sd,
      comp1_sd = comp[1], comp2_sd = comp[2], comp3_sd = comp[3],
      sds1 = sds[1], sds2 = sds[2], sds3 = sds[3],
      smooth_range = diff(range(tot)), smooth_sd = sd(tot))
  }
}
print(do.call(rbind, rows), digits = 3)
