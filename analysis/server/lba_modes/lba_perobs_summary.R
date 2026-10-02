options(width = 200)
scratch <- "lba_modes"
o <- readRDS(file.path(scratch, "lba_perobs.rds"))
cat("== validation (diff must be constant across chains)\n"); print(o$val, digits = 8)
r <- o$res
r$d <- r$llB - r$llA
r$side <- ifelse(r$Illusion_StrengthZ > 0, "conflicting", ifelse(r$Illusion_StrengthZ < 0, "facilitating", "none"))
r$diffbin <- cut(r$Illusion_DifferenceZ, quantile(r$Illusion_DifferenceZ, 0:4 / 4), include.lowest = TRUE, labels = paste0("D", 1:4))
cat("\n== total B - A:", sum(r$d), "\n")
cat("\n== by response\n"); print(aggregate(d ~ Error, r, function(x) c(sum = sum(x), per_obs = mean(x))))
cat("\n== by cell (sum of B - A, error rate)\n")
cell <- aggregate(cbind(d, Error) ~ side + diffbin, r, function(x) c(s = sum(x), m = mean(x)))
cell <- do.call(data.frame, cell); print(cell[order(cell$side, cell$diffbin), c("side", "diffbin", "d.s", "Error.m")], digits = 4)
cat("\n== by cell x response (sum of B - A)\n")
print(xtabs(d ~ interaction(side, diffbin) + Error, r), digits = 4)
r$rtbin <- cut(r$RT, c(0, .4, .5, .6, .8, 1.2, 2, Inf))
cat("\n== by RT bin x response\n"); print(xtabs(d ~ rtbin + Error, r), digits = 4)

pp <- aggregate(d ~ Participant, r, sum)
pp <- pp[order(-abs(pp$d)), ]
cat("\n== participants: share of |B-A| in top 1/5/10%:\n")
tot <- sum(abs(pp$d)); k <- nrow(pp)
print(round(c(top1 = sum(abs(pp$d[1:ceiling(k * .01)])), top5 = sum(abs(pp$d[1:ceiling(k * .05)])),
              top10 = sum(abs(pp$d[1:ceiling(k * .10)]))) / tot, 3))
cat("participants favouring B:", mean(pp$d > 0), " quantiles of per-ppt B-A:\n"); print(round(quantile(pp$d, c(0, .01, .05, .25, .5, .75, .95, .99, 1)), 2))
pe <- aggregate(cbind(Error, RT) ~ Participant, r, mean)
pp <- merge(pp, pe)
cat("cor(per-ppt B-A, ppt error rate):", cor(pp$d, pp$Error), " with mean RT:", cor(pp$d, pp$RT), "\n")

cat("\n== posterior-mean parameters by cell and mode\n")
v <- c("mu", "driftone", "sigmaone", "sigmabias", "boundary", "ndt")
pm <- aggregate(r[, c(paste0(v, "_A"), paste0(v, "_B"))], r[c("side", "diffbin")], mean)
print(pm[order(pm$side, pm$diffbin), ], digits = 3)
# realised mean drift of the error accumulator: mean of N(v, s) truncated at 0
tm <- function(v, s) v + s * dnorm(v / s) / pnorm(v / s)
r$realA <- tm(r$driftone_A, r$sigmaone_A); r$realB <- tm(r$driftone_B, r$sigmaone_B)
r$rayA <- abs(r$driftone_A) / r$sigmaone_A^2; r$rayB <- abs(r$driftone_B) / r$sigmaone_B^2
cat("\n== error accumulator: v/s, realised mean drift, |v|/s^2 by cell\n")
print(aggregate(cbind(zA = driftone_A / sigmaone_A, zB = driftone_B / sigmaone_B, realA, realB, rayA, rayB) ~ side + diffbin, r, median), digits = 3)

cat("\n== posterior predictive by cell: error rate and RT quantiles (observed, A, B)\n")
s <- o$sims
s$side <- r$side[s$row]; s$diffbin <- r$diffbin[s$row]; s$obs_err <- r$Error[s$row]
summ <- function(rt, resp, side, diffbin) {
  dd <- data.frame(rt, resp, side, diffbin)
  do.call(rbind, lapply(split(dd, list(dd$side, dd$diffbin), drop = TRUE), function(x) {
    data.frame(side = x$side[1], diffbin = x$diffbin[1], err = mean(x$resp),
               c_q50 = median(x$rt[x$resp == 0]), c_q90 = quantile(x$rt[x$resp == 0], .9),
               e_q10 = quantile(x$rt[x$resp == 1], .1), e_q50 = median(x$rt[x$resp == 1]),
               e_q90 = quantile(x$rt[x$resp == 1], .9))
  }))
}
obs <- summ(r$RT, r$Error, r$side, r$diffbin); obs$who <- "obs"
sA <- summ(s$rt[s$mode == "A"], s$response[s$mode == "A"], s$side[s$mode == "A"], s$diffbin[s$mode == "A"]); sA$who <- "A"
sB <- summ(s$rt[s$mode == "B"], s$response[s$mode == "B"], s$side[s$mode == "B"], s$diffbin[s$mode == "B"]); sB$who <- "B"
ppc <- rbind(obs, sA, sB)
ppc <- ppc[order(ppc$side, ppc$diffbin, match(ppc$who, c("obs", "A", "B"))), c("side", "diffbin", "who", "err", "c_q50", "c_q90", "e_q10", "e_q50", "e_q90")]
rownames(ppc) <- NULL
print(ppc, digits = 3)
saveRDS(list(cell = cell, pm = pm, ppc = ppc, pp = pp), file.path(scratch, "lba_perobs_summary.rds"))
