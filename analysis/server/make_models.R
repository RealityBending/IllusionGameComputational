# =========================================================================
# Fitting -----------------------------------------------------------------
# =========================================================================
# TODO: this needs to be run on the cluster

library(brms)
library(cogmod) # remotes::install_github("DominiqueMakowski/cogmod")
library(dplyr)

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
total_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "2"))

chains_per_node <- 2
threads_per_chain <- total_cores / chains_per_node # 16 cores / 2 chains = 8 threads per chain

warmup <- 1000
iter <- warmup + 200


df <- rbind(
  read.csv("https://raw.githubusercontent.com/RealityBending/IllusionGameComputational/refs/heads/main/data/illusion_part1.csv"),
  read.csv("https://raw.githubusercontent.com/RealityBending/IllusionGameComputational/refs/heads/main/data/illusion_part2.csv"),
  read.csv("https://raw.githubusercontent.com/RealityBending/IllusionGameComputational/refs/heads/main/data/illusion_part3.csv")
)
df$Illusion_Difference <- abs(df$Illusion_Difference)
df$Illusion_Effect <- factor(ifelse(df$Illusion_Strength >= 0, "Conflicting", "Facilitating"), levels = c("Conflicting", "Facilitating"))

# Normalize predictors
# summarize(df,
#   MinDiff = min(Illusion_Difference),
#   MaxDiff = max(Illusion_Difference),
#   MinStrength = min(Illusion_Strength),
#   MinStrength_abs = min(abs(Illusion_Strength)),
#   MaxStrength = max(abs(Illusion_Strength)),
#   .by = c("Illusion_Effect", "Illusion_Type")
# )
df <- mutate(df,
  Illusion_DifferenceZ = as.numeric(datawizard::normalize(Illusion_Difference)),
  Illusion_StrengthZ = sign(Illusion_Strength) * as.numeric(datawizard::normalize(abs(Illusion_Strength))),
  .by = "Illusion_Type"
)
# plot(df$Illusion_Difference, df$Illusion_DifferenceZ, col = as.numeric(as.factor(df$Illusion_Type)), pch = 19, cex = 0.5)
# plot(df$Illusion_Strength, df$Illusion_StrengthZ, col = as.numeric(as.factor(df$Illusion_Type)), pch = 19, cex = 0.5)


# TODO: data subset (to be removed in the final version)
df <- df[
  df$Participant %in% unique(df$Participant)[1:30],
]

# =========================================================================
# GAMs --------------------------------------------------------------------
# =========================================================================

# LNR ---------------------------------------------------------------------

# Formula
f <- bf(
  RT | dec(Error) ~ t2(Illusion_DifferenceZ, Illusion_StrengthZ,
    k = c(5, 5),
    bs = c("cr", "cr")
  ),
  nuone ~ t2(Illusion_DifferenceZ, Illusion_StrengthZ,
    k = c(5, 5),
    bs = c("cr", "cr")
  ),
  sigmazero ~ t2(Illusion_DifferenceZ, Illusion_StrengthZ,
    k = c(5, 5),
    bs = c("cr", "cr")
  ),
  sigmaone ~ t2(Illusion_DifferenceZ, Illusion_StrengthZ,
    k = c(5, 5),
    bs = c("cr", "cr")
  ),
  ndt ~ t2(Illusion_DifferenceZ, Illusion_StrengthZ,
    k = c(5, 5),
    bs = c("cr", "cr")
  ),
  poutlier ~ 1,
  family = cogmod_lnr()
)


# Informative priors
# brms::get_prior(f, data = df)
priors <- c(
  cogmod_priors(f, df[df$Illusion_Type == "MullerLyer", ]),
  brms::prior("normal(0, 1)", class = "b", dpar = "", coef = ""),
  brms::prior("normal(0, 1)", class = "b", dpar = "nuone", coef = ""),
  brms::prior("normal(0, 1)", class = "b", dpar = "sigmazero", coef = ""),
  brms::prior("normal(0, 1)", class = "b", dpar = "sigmaone", coef = ""),
  brms::prior("normal(0, 1)", class = "b", dpar = "ndt", coef = ""),
  replace = TRUE
)


# Fit
gam_lnr_muller <- brm(f,
  data = df[df$Illusion_Type == "MullerLyer", ],
  prior = priors,
  init = cogmod_inits(f, df[df$Illusion_Type == "MullerLyer", ]),
  stanvars = cogmod_stanvars(f),
  backend = "cmdstanr",
  warmup = warmup,
  iter = iter,
  algorithm = "pathfinder", chains = 16, single_path_draws = 4000, max_lbfgs_iters = 8000, threads = 8
  # chains = chains_per_node,
  # cores = chains_per_node,
  # threads = threading(threads_per_chain),
  # save_pars = save_pars(all = TRUE),
  # algorithm = "sampling"
)

# TODO: WAIC instead of loo will probably be faster.
# m_lnr_muller <- brms::add_criterion(m_lnr_muller, "loo", moment_match = FALSE)
# m_lnr_muller <- brms::add_criterion(m_lnr_muller, "waic")


saveRDS(gam_lnr_muller, sprintf("../models/gam_lnr_muller_%d.rds", task_id))

# =========================================================================
# LINEAR ------------------------------------------------------------------
# =========================================================================

# Formula
f <- bf(
  RT | dec(Error) ~ Illusion_Effect * Illusion_DifferenceZ * abs(Illusion_StrengthZ),
  nuone ~ Illusion_Effect * Illusion_DifferenceZ * abs(Illusion_StrengthZ),
  sigmazero ~ Illusion_Effect * Illusion_DifferenceZ * abs(Illusion_StrengthZ),
  sigmaone ~ Illusion_Effect * Illusion_DifferenceZ * abs(Illusion_StrengthZ),
  ndt ~ Illusion_Effect * Illusion_DifferenceZ * abs(Illusion_StrengthZ),
  poutlier ~ 1,
  family = cogmod_lnr()
)


# Informative priors
# brms::get_prior(f, data = df)
priors <- c(
  cogmod_priors(f, df[df$Illusion_Type == "MullerLyer", ]),
  brms::prior("normal(0, 1)", class = "b", dpar = "", coef = ""),
  brms::prior("normal(0, 1)", class = "b", dpar = "nuone", coef = ""),
  brms::prior("normal(0, 1)", class = "b", dpar = "sigmazero", coef = ""),
  brms::prior("normal(0, 1)", class = "b", dpar = "sigmaone", coef = ""),
  brms::prior("normal(0, 1)", class = "b", dpar = "ndt", coef = ""),
  replace = TRUE
)


# Fit
m1_lnr_muller <- brm(f,
  data = df[df$Illusion_Type == "MullerLyer", ],
  prior = priors,
  init = cogmod_inits(f, df[df$Illusion_Type == "MullerLyer", ]),
  stanvars = cogmod_stanvars(f),
  backend = "cmdstanr",
  warmup = warmup,
  iter = iter,
  algorithm = "pathfinder", chains = 16, single_path_draws = 4000, max_lbfgs_iters = 8000, threads = 8
  # chains = chains_per_node,
  # cores = chains_per_node,
  # threads = threading(threads_per_chain),
  # save_pars = save_pars(all = TRUE),
  # algorithm = "sampling"
)

# TODO: WAIC instead of loo will probably be faster.
# m_lnr_muller <- brms::add_criterion(m_lnr_muller, "loo", moment_match = FALSE)
# m_lnr_muller <- brms::add_criterion(m_lnr_muller, "waic")


saveRDS(m1_lnr_muller, sprintf("../models/m1_lnr_muller_%d.rds", task_id))
