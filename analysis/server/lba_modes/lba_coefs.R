options(width = 200)
m <- readRDS("../models/gam_lba_MullerLyer.rds")
a <- as.array(m$fit)
nm <- dimnames(a)[[3]]
pick <- nm[startsWith(nm, "bs_") | startsWith(nm, "sds_") | startsWith(nm, "Intercept")]
tab <- t(apply(a[, , pick], c(2, 3), mean))
colnames(tab) <- paste0(c("A", "A", "B", "A", "B", "A", "B", "A"), 1:8)
rownames(tab) <- sub("t2Illusion_DifferenceZIllusion_StrengthZ_", "t2_", rownames(tab))
print(round(tab, 2))
# Rhat restricted to within-mode chains: does each mode mix on its own?
rh <- function(ch) {
  s <- posterior::summarise_draws(posterior::as_draws_array(a[, ch, c(pick, "lp__")]), "rhat")
  setNames(s$rhat, s$variable)
}
cat("\nmax Rhat within A chains:", round(max(rh(c(1, 2, 4, 6, 8)), na.rm = TRUE), 3),
    " within B chains:", round(max(rh(c(3, 5, 7)), na.rm = TRUE), 3), "\n")
rb <- rh(c(3, 5, 7)); print(round(sort(rb, decreasing = TRUE)[1:8], 2))
ra <- rh(c(1, 2, 4, 6, 8)); print(round(sort(ra, decreasing = TRUE)[1:8], 2))
