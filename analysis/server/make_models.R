# =========================================================================
# Fitting -----------------------------------------------------------------
# =========================================================================
# TODO: this needs to be run on the cluster

library(brms)
library(cogmod) # remotes::install_github("DominiqueMakowski/cogmod")

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
total_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "2"))

chains_per_node <- 2
threads_per_chain <- total_cores / chains_per_node  # 16 cores / 2 chains = 8 threads per chain

warmup <- 1000
iter <- warmup + 200


df <- rbind(
  read.csv("https://raw.githubusercontent.com/RealityBending/IllusionGameComputational/refs/heads/main/data/illusion_part1.csv"),
  read.csv("https://raw.githubusercontent.com/RealityBending/IllusionGameComputational/refs/heads/main/data/illusion_part2.csv"),
  read.csv("https://raw.githubusercontent.com/RealityBending/IllusionGameComputational/refs/heads/main/data/illusion_part3.csv")
)
df$Illusion_Difference <- abs(df$Illusion_Difference)

df <- df[1:10001, ] # TODO: remove this when running full models on the cluster.


# LNR ---------------------------------------------------------------------

# Formula
f <- bf(
  RT | dec(Error) ~ poly(Illusion_Difference, 2) * poly(Illusion_Strength, 2),
  nuone ~ poly(Illusion_Difference, 2) * poly(Illusion_Strength, 2),
  sigmazero ~ poly(Illusion_Difference, 2) * poly(Illusion_Strength, 2),
  sigmaone ~ poly(Illusion_Difference, 2) * poly(Illusion_Strength, 2),
  tau ~ Illusion_Difference + Illusion_Strength,
  minrt = min(df$RT),
  family = lnr()
)


# Weakly informative priors
# brms::get_prior(f, data = df)
priors <- c(
  brms::set_prior("normal(0.3, 1)", class = "Intercept", dpar = ""),
  brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "nuone"),
  brms::set_prior("normal(-0.5, 0.5)", class = "Intercept", dpar = "sigmazero"),
  brms::set_prior("normal(0, 0.5)", class = "Intercept", dpar = "sigmaone"),
  brms::set_prior("normal(3, 10)", class = "Intercept", dpar = "tau")
)
for (p in c(
  "polyIllusion_Difference21", "polyIllusion_Strength21",
  "polyIllusion_Difference22", "polyIllusion_Strength22",
  "polyIllusion_Difference21:polyIllusion_Strength21", "polyIllusion_Difference21:polyIllusion_Strength22",
  "polyIllusion_Difference22:polyIllusion_Strength21", "polyIllusion_Difference22:polyIllusion_Strength22"
)) {
  for (dpar in c("", "nuone", "sigmazero", "sigmaone")) {
    if (grepl(":", p)) { # Interaction
      prior <- "normal(0, 1)"
    } else if (grepl("22", p)) {
      prior <- "normal(0, 5)"
    } else if (grepl("21", p)) {
      prior <- "normal(0, 10)"
    } else {
      prior <- "normal(0, 5)"
    }
    priors <- c(priors, brms::set_prior(prior, class = "b", dpar = dpar, coef = p))
  }
}
priors <- c(priors, brms::set_prior("normal(0, 1)", class = "b", dpar = "tau", coef = c("Illusion_Difference", "Illusion_Strength")))
priors <- brms::validate_prior(priors, f, data = df)


# Fit
m_lnr_muller <- brm(f,
  data = df[df$Illusion_Type == "MullerLyer", ], # muller-lyer only? maybe add 1-2 other types? ebbinghaus and vertical-horizontal
  init = 0,
  prior = priors,
  family = lnr(),
  stanvars = lnr_stanvars(),
  warmup = warmup,
  iter = iter,
  chains = chains_per_node,
  cores = chains_per_node,
  threads = threading(threads_per_chain),
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE),
  algorithm = "sampling" #algorithm = "pathfinder", max_lbfgs_iters = 6000, draws = 4000
)

# TODO: WAIC instead of loo will probably be faster.
# m_lnr_muller <- brms::add_criterion(m_lnr_muller, "loo", moment_match = FALSE)
# m_lnr_muller <- brms::add_criterion(m_lnr_muller, "waic")

dir.create("models", showWarnings = FALSE)

saveRDS(m_lnr_muller, sprintf("models/m_lnr_muller_%d.rds", task_id))



# BayesFlow-Assisted Process ----------------------------------------------


# # Formula
# f <- bf(
#   RT | dec(Error) ~ Illusion_Difference,
#   nuone ~ Illusion_Difference,
#   sigmazero ~ Illusion_Difference,
#   sigmaone ~ Illusion_Difference,
#   tau ~ Illusion_Difference,
#   minrt = min(df$RT),
#   family = lnr()
# )
#
# # Fit
# m_lnr_muller <- brm(f,
#   data = df[df$Illusion_Type == "MullerLyer", ],
#   init = 0,
#   family = lnr(),
#   stanvars = lnr_stanvars(),
#   iter = 4000,
#   chains = 2, threads = threading(2),
#   backend = "cmdstanr",
#   save_pars = save_pars(all = TRUE),
#   algorithm = "pathfinder", max_lbfgs_iters = 500, draws = 100
# )
#
# insight::find_parameters(m_lnr_muller, effects = "all", component = "all")
# brms::get_prior(f, data = df)
