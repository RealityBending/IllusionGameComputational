# =========================================================================
# Fitting -----------------------------------------------------------------
# =========================================================================
# TODO: this needs to be run on the cluster

library(brms)
library(cogmod) # remotes::install_github("DominiqueMakowski/cogmod")

df <- rbind(
  read.csv("https://raw.githubusercontent.com/RealityBending/IllusionGameComputational/refs/heads/main/data/illusion_part1.csv"),
  read.csv("https://raw.githubusercontent.com/RealityBending/IllusionGameComputational/refs/heads/main/data/illusion_part2.csv")
) |>
  mutate(Illusion_Difference = abs(Illusion_Difference))

df <- df[1:100000, ]


# LNR ---------------------------------------------------------------------


f <- bf(
  RT | dec(Error) ~ Illusion_Difference * Illusion_Strength,
  nuone ~ Illusion_Difference * Illusion_Strength,
  sigmazero ~ Illusion_Difference,
  sigmaone ~ Illusion_Difference,
  tau ~ 1,
  minrt = min(df$RT),
  family = lnr()
)

m_lnr_muller <- brm(f,
  data = df[df$Illusion_Type == "MullerLyer", ],
  init = 0,
  family = lnr(),
  stanvars = lnr_stanvars(),
  chains = 4, iter = 4000, backend = "cmdstanr",
  save_pars = save_pars(all = TRUE),
  algorithm = "pathfinder"
)

m_lnr_muller <- brms::add_criterion(m_lnr_muller, "loo", moment_match = FALSE)

saveRDS(m_lnr_muller, "models/m_lnr_muller.rds")
