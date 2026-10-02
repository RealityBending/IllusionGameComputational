# Split lp__ of each chain of the combined gam_lba fit into
#   log-likelihood + lprior + std_normal(z) + std_normal(zs) + log-Jacobians
# z and zs are not saved (save_pars(all = FALSE)), so they are rebuilt as
# r / sd and s / sds. Constants are dropped; only differences between chains
# matter.
m <- readRDS("../models/gam_lba_MullerLyer.rds")
a <- as.array(m$fit)
nm <- dimnames(a)[[3]]
nch <- dim(a)[2]

dpars <- c("", "driftone", "sigmaone", "sigmabias", "boundary", "ndt", "poutlier")
out <- list()
for (ch in seq_len(nch)) {
  x <- a[, ch, ]
  zterm <- 0; jac <- 0; zsterm <- 0
  for (dp in dpars) {
    sdn <- if (dp == "") "sd_Participant__Intercept" else paste0("sd_Participant__", dp, "_Intercept")
    rpref <- if (dp == "") "r_Participant[" else paste0("r_Participant__", dp, "[")
    rcols <- nm[startsWith(nm, rpref)]
    z <- x[, rcols] / x[, sdn]
    zterm <- zterm + rowSums(-0.5 * z^2)
    jac <- jac + log(x[, sdn])
    if (dp != "poutlier") {
      tag <- if (dp == "") "" else paste0("_", dp)
      for (k in 1:3) {
        sdsn <- paste0("sds", tag, "_t2Illusion_DifferenceZIllusion_StrengthZ_", k)
        spref <- paste0("s", tag, "_t2Illusion_DifferenceZIllusion_StrengthZ_", k, "[")
        scols <- nm[startsWith(nm, spref)]
        zs <- x[, scols, drop = FALSE] / x[, sdsn]
        zsterm <- zsterm + rowSums(-0.5 * zs^2)
        jac <- jac + log(x[, sdsn])
      }
    }
  }
  ll <- x[, "lp__"] - x[, "lprior"] - zterm - zsterm - jac
  out[[ch]] <- data.frame(chain = ch, lp = mean(x[, "lp__"]), lprior = mean(x[, "lprior"]),
                          z = mean(zterm), zs = mean(zsterm), jac = mean(jac), loglik = mean(ll),
                          loglik_sd = sd(ll))
}
res <- do.call(rbind, out)
res$mode <- ifelse(res$lp > 37500, "B", "A")
print(res, digits = 7)
cat("\nmode means:\n")
print(aggregate(cbind(lp, lprior, z, zs, jac, loglik) ~ mode, res, mean), digits = 7)
saveRDS(res, "lba_modes/lba_decompose.rds")
