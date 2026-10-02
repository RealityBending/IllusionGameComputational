# The error rate implied by cogmod_lba2's starting values, as a function of the
# starting driftone, with everything else at the registry init
# (mu = 3, sigmazero = sigmaone = 1, sigmabias = 0.5, boundary = 0.5).
# P(response 1) = integral over t of the race's defective density for response 1.
suppressMessages(library(cogmod))
p_err <- function(v1, mu = 3, s0 = 1, s1 = 1, A = 0.5, B = 0.5) {
  p <- list(mu = mu, driftone = v1, sigmazero = s0, sigmaone = s1, sigmabias = A, boundary = B)
  integrate(function(t) exp(cogmod:::.lba2_ldens(t, rep(1, length(t)), p)), 0, Inf,
            rel.tol = 1e-8, subdivisions = 1000)$value
}
cat("check: P(0) + P(1) at driftone = 1:",
    p_err(1) + integrate(function(t) exp(cogmod:::.lba2_ldens(t, rep(0, length(t)),
      list(mu = 3, driftone = 1, sigmazero = 1, sigmaone = 1, sigmabias = .5, boundary = .5))), 0, Inf)$value, "\n")
grid <- c(3, 2, 1.5, 1, 0.5, 0, -0.5, -1, -2)
print(data.frame(driftone = grid, p_error = round(sapply(grid, p_err), 4)))
solve_for <- function(target) uniroot(function(v) p_err(v) - target, c(-10, 3))$root
targets <- c(0.01, 0.05, 0.10, 0.237, 0.40)
print(data.frame(target_error_rate = targets, driftone_start = round(sapply(targets, solve_for), 3)))
