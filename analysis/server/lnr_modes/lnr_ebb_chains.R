# Is gam_lnr_ebbinghaus's non-convergence the gam_lba pattern (chains sitting
# stably in distinct modes) or something else (chains not stationary, stuck)?
# Per chain: lp__ level and drift, sampler state, and the population-level
# parameters whose chain means disagree most. Runs on the shards, before any
# combine. Usage: Rscript lnr_ebb_chains.R <models_dir> [suffix]
# suffix "" = the 2026-10-11 refit (cogmod 0.3.4); "_20261002" = the old fit.
options(width = 220)
args <- commandArgs(TRUE)
dir <- args[1]
suffix <- if (length(args) > 1) args[2] else ""
files <- file.path(dir, sprintf("gam_lnr_ebbinghaus_Ebbinghaus_%d%s.rds", 1:4, suffix))

arrs <- list()
for (i in seq_along(files)) {
  m <- readRDS(files[i])
  a <- as.array(m$fit)
  nm <- dimnames(a)[[3]]
  keep <- nm[grepl("^(b_|bs_|sds_|sd_|Intercept|lp__|lprior)", nm)]
  a <- a[, , keep, drop = FALSE]
  dimnames(a)[[2]] <- paste0("s", i, "c", seq_len(dim(a)[2]))
  arrs[[i]] <- a
  rm(m); gc()
}
# bind along chains without abind: put chains last, concatenate, put them back
perm <- lapply(arrs, aperm, c(1, 3, 2))
a <- aperm(array(unlist(perm), dim = c(dim(perm[[1]])[1:2], sum(sapply(perm, function(x) dim(x)[3]))),
                 dimnames = list(NULL, dimnames(perm[[1]])[[2]], unlist(lapply(perm, function(x) dimnames(x)[[3]])))),
           c(1, 3, 2))
nch <- dim(a)[2]; nd <- dim(a)[1]
cat("suffix '", suffix, "': ", nch, " chains x ", nd, " draws, ", dim(a)[3], " population-level params\n\n", sep = "")

# lp__ per chain: level, spread, and drift (mean of 4 consecutive quarters)
q <- split(seq_len(nd), cut(seq_len(nd), 4, labels = FALSE))
lp <- a[, , "lp__"]
tab <- data.frame(chain = colnames(lp), lp_mean = colMeans(lp), lp_sd = apply(lp, 2, sd),
                  t(sapply(seq_len(nch), function(j) sapply(q, function(ix) mean(lp[ix, j])))),
                  lprior = colMeans(a[, , "lprior"]))
names(tab)[4:7] <- paste0("lp_Q", 1:4)
# within-chain split-Rhat: > 1.05 means the chain itself is not stationary
within <- sapply(seq_len(nch), function(j) {
  r <- apply(a[, j, , drop = FALSE], 3, function(x) posterior::rhat(matrix(x, ncol = 1)))
  c(max_within_rhat = max(r, na.rm = TRUE), n_within_gt105 = sum(r > 1.05, na.rm = TRUE))
})
tab <- cbind(tab, t(within))
print(tab, digits = 6, row.names = FALSE)

# Which parameters separate the chains: between-chain SD of chain means over
# the mean within-chain SD. Big ratio = chains in different places.
cm <- apply(a, c(2, 3), mean)
wsd <- apply(apply(a, c(2, 3), sd), 2, mean)
ratio <- apply(cm, 2, sd) / wsd
top <- names(sort(ratio, decreasing = TRUE))[1:25]
cat("\nTop 25 parameters by between/within chain spread (chain means):\n")
short <- function(x) sub("t2Illusion_DifferenceZIllusion_StrengthZ", "t2", x)
out <- t(round(cm[, top], 3)); rownames(out) <- paste0(short(top), "  [", round(ratio[top], 1), "]")
print(out)

# Key scalars for every chain, whatever their ratio
keyn <- dimnames(a)[[3]][grepl("^(Intercept|sd_Participant)", dimnames(a)[[3]])]
cat("\nIntercepts and participant SDs, chain means:\n")
out2 <- t(round(cm[, keyn], 3)); rownames(out2) <- keyn
print(out2)
