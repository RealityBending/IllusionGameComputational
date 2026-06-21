# =========================================================================
# Fitting -----------------------------------------------------------------
# =========================================================================
# TODO: this needs to be run on the cluster

library(brms)
library(cogmod) # remotes::install_github("DominiqueMakowski/cogmod")

df <- rbind(
  read.csv("https://raw.githubusercontent.com/RealityBending/IllusionGameComputational/refs/heads/main/data/illusion_part1.csv"),
  read.csv("https://raw.githubusercontent.com/RealityBending/IllusionGameComputational/refs/heads/main/data/illusion_part2.csv")
)
df$Illusion_Difference <- abs(df$Illusion_Difference)

df <- df[1:150000, ]


# LNR ---------------------------------------------------------------------


f <- bf(
  RT | dec(Error) ~ poly(Illusion_Difference, 2) * poly(Illusion_Strength, 2),
  nuone ~ poly(Illusion_Difference, 2) * poly(Illusion_Strength, 2),
  sigmazero ~ poly(Illusion_Difference, 2) * poly(Illusion_Strength, 2),
  sigmaone ~ poly(Illusion_Difference, 2) * poly(Illusion_Strength, 2),
  tau ~ 1,
  minrt = min(df$RT),
  family = lnr()
)

# brms::get_prior(f, data = df)

priors <- brms::empty_prior()
for (p in c(
  "polyIllusion_Difference21", "polyIllusion_Strength21", "polyIllusion_Difference22", "polyIllusion_Strength22",
  "polyIllusion_Difference21:polyIllusion_Strength21", "polyIllusion_Difference21:polyIllusion_Strength22",
  "polyIllusion_Difference22:polyIllusion_Strength21", "polyIllusion_Difference22:polyIllusion_Strength22"
)) {
  for (dpar in c("", "nuone", "sigmazero", "sigmaone")) {
    if (grepl(":", p)) { # Interaction
      prior <- "normal(0, 0.1)"
    } else if (grepl("22", p)) {
      prior <- "normal(0, 1)"
    } else if (grepl("21", p)) {
      prior <- "normal(0, 10)"
    } else {
      prior <- "normal(0, 5)"
    }
    priors <- c(priors, brms::set_prior(prior, class = "b", dpar = dpar, coef = p))
  }
}
priors <- brms::validate_prior(priors, f, data = df)



m_lnr_muller <- brm(f,
  data = df[df$Illusion_Type == "MullerLyer", ],
  init = 0,
  prior = priors,
  family = lnr(),
  stanvars = lnr_stanvars(),
  iter = 4000,
  chains = 8, threads = threading(4),
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE),
  algorithm = "pathfinder", max_lbfgs_iters = 4000, draws = 4000
)

# TODO: WAIC instead of loo will probably be faster.
# m_lnr_muller <- brms::add_criterion(m_lnr_muller, "loo", moment_match = FALSE)
m_lnr_muller <- brms::add_criterion(m_lnr_muller, "waic")

saveRDS(m_lnr_muller, "models/m_lnr_muller.rds")
