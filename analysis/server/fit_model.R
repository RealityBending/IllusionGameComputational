# =========================================================================
# Fit ONE model, as one shard of an array job
# =========================================================================
# Which model is IGC_MODEL (see models.R); which shard is SLURM_ARRAY_TASK_ID.
# Each shard writes <model>_<illusion>_<task>.rds into IGC_MODELS_DIR and
# combine_model.R merges them. Submitted by `./hpc fit <model>`.
#
# Production fits are cold starts: cogmod_inits() plus Stan's own adaptation.
# The warm-start machinery (inits / inverse metric / step size carried over
# from a pilot fit) was removed on 2026-09-18 -- the decision is to run cold,
# and cogmod 0.3.3 fixed the two things that made cold starts expensive.

library(brms)
library(cogmod) # remotes::install_github("DominiqueMakowski/cogmod@dev")
library(dplyr)

# cogmod 0.3.3 fixed the LNR/LogNormal tail gradient (chains rejecting their
# initial value, divergences) and the jitter that started smooths far from
# their targets. Fitting with an older one wastes days, so fail in the first
# second rather than the twentieth hour. See cogmod_inits_issue.md.
if (utils::packageVersion("cogmod") < "0.3.3") {
  stop("cogmod ", utils::packageVersion("cogmod"), " is too old; need >= 0.3.3. ",
       "Run  ./hpc install cogmod", call. = FALSE)
}

source("models.R")

spec <- igc_model(Sys.getenv("IGC_MODEL", unset = ""))

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
total_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "2"))

# Chains per array task. Each chain is a separate OS process with its own copy
# of the data and AD stack, so this is the main driver of per-task memory.
# Two is the working compromise: threading scales nearly linearly, so more
# chains per task buys no throughput, while a single chain would put the whole
# task on one chain's survival.
chains_per_task <- as.integer(Sys.getenv("IGC_CHAINS", unset = "2"))
threads_per_chain <- total_cores / chains_per_task # 16 cores / 2 chains = 8 threads

# Output directory (set by ./hpc; falls back to ./models locally)
models_dir <- Sys.getenv("IGC_MODELS_DIR", unset = "models")
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

# "never" makes a resubmitted array skip shards whose .rds already exists, so a
# job killed at the wall resumes instead of restarting. The catch is that a
# stale shard is reused silently even if the formula or data changed -- so
# combine_model.R deletes shards once they are combined, and you can force a
# clean refit with:  IGC_FILE_REFIT=always ./hpc fit <model>
file_refit <- Sys.getenv("IGC_FILE_REFIT", unset = "never")

warmup <- as.integer(Sys.getenv("IGC_WARMUP", unset = "1000"))
iter <- warmup + as.integer(Sys.getenv("IGC_SAMPLES", unset = "500"))

cat(sprintf(
  "config: model=%s illusion=%s shard=%d warmup=%d samples=%d chains=%d threads/chain=%g refit=%s\n",
  spec$name, spec$illusion, task_id, warmup, iter - warmup,
  chains_per_task, threads_per_chain, file_refit
))


# Data --------------------------------------------------------------------
# Read straight from GitHub: the compute nodes have outbound internet, so
# nothing but code needs pushing.

df <- rbind(
  read.csv("https://raw.githubusercontent.com/RealityBending/IllusionGameComputational/refs/heads/main/data/illusion_part1.csv"),
  read.csv("https://raw.githubusercontent.com/RealityBending/IllusionGameComputational/refs/heads/main/data/illusion_part2.csv"),
  read.csv("https://raw.githubusercontent.com/RealityBending/IllusionGameComputational/refs/heads/main/data/illusion_part3.csv")
)
df$Illusion_Difference <- abs(df$Illusion_Difference)
df$Illusion_Effect <- factor(
  ifelse(df$Illusion_Strength >= 0, "Conflicting", "Facilitating"),
  levels = c("Conflicting", "Facilitating")
)

# Normalize the predictors within illusion, keeping the sign of the strength
# (which side of the illusion the trial is on).
df <- mutate(df,
  Illusion_DifferenceZ = as.numeric(datawizard::normalize(Illusion_Difference)),
  Illusion_StrengthZ = sign(Illusion_Strength) * as.numeric(datawizard::normalize(abs(Illusion_Strength))),
  .by = "Illusion_Type"
)

# Subset size. Production is "all"; set IGC_NPARTICIPANTS=30 for a smoke test.
n_participants <- Sys.getenv("IGC_NPARTICIPANTS", unset = "all")
if (!identical(tolower(n_participants), "all")) {
  keep <- unique(df$Participant)[seq_len(as.integer(n_participants))]
  df <- df[df$Participant %in% keep, ]
}

data <- df[df$Illusion_Type == spec$illusion, ]
cat("participants:", length(unique(data$Participant)), " rows:", nrow(data), "\n")


# Fit ---------------------------------------------------------------------

f <- spec$formula()

# cogmod's family-aware priors, plus a standard normal on the slopes brms would
# otherwise leave flat. poutlier keeps cogmod's own prior.
#
# Only the FLAT ones. cogmod already names a slope prior for the dpars it knows
# to be hard to identify, and it is tighter than normal(0, 1) on purpose:
# normal(0, 0.2) for `ndt` and normal(0, 0.5) for `sigmadrift` / `sigmabias` /
# `sigmandt`. Until 2026-09-18 this loop overwrote those with normal(0, 1)
# unconditionally, i.e. it widened by 5x and 2x exactly the priors cogmod had
# tightened. Measured over gam_ddm7: `mu`, `bias` and `boundary` arrive flat and
# want filling; the other four arrive set and must be left alone.
priors <- cogmod_priors(f, data)
for (par in c("", setdiff(unique(priors$dpar), c("poutlier", "")))) {
  # The blanket `b` row for this dpar: class "b" with no coef, no group.
  blanket <- priors$class == "b" & priors$dpar == par &
    !nzchar(priors$coef) & !nzchar(priors$group)
  if (any(blanket & nzchar(priors$prior))) next # cogmod set it; keep it
  priors <- c(priors, brms::prior_string("normal(0, 1)", class = "b", dpar = par),
              replace = TRUE)
}

# One-line summary per fit so runs can be compared from the .out logs alone.
report_fit <- function(m, name, wall_min) {
  np <- brms::nuts_params(m) # post-warmup draws only
  stat <- vapply(split(np$Value, np$Parameter), mean, numeric(1))
  rh <- brms::rhat(m)
  ne <- brms::neff_ratio(m)
  cat(sprintf(
    "REPORT %s | wall %.1f min | n_leapfrog %.0f | treedepth %.2f | stepsize %.3g | divergent %.3f | accept %.2f | max Rhat %.3f | min neff_ratio %.3f | n params %d\n",
    name, wall_min, stat[["n_leapfrog__"]], stat[["treedepth__"]], stat[["stepsize__"]],
    stat[["divergent__"]], stat[["accept_stat__"]], max(rh, na.rm = TRUE),
    min(ne, na.rm = TRUE), length(rh)
  ))
  md <- attr(m$fit, "metadata")
  if (!is.null(md$time)) {
    cat("REPORT time per chain (s):\n")
    print(md$time)
  }
  if (!is.null(md$step_size)) {
    cat("REPORT adapted step sizes:", format(unlist(md$step_size), digits = 3), "\n")
  }
  invisible(stat)
}

t0 <- Sys.time()
m <- brm(f,
  data = data,
  prior = priors,
  init = cogmod_inits(f, data),
  stanvars = cogmod_stanvars(f),
  backend = "cmdstanr",
  warmup = warmup,
  iter = iter,
  algorithm = "sampling",
  chains = chains_per_task,
  cores = chains_per_task,
  threads = threading(threads_per_chain),
  # save_pars(all = TRUE) dropped: it keeps every Stan parameter, including
  # the latents behind ~6 x n_participants group-level coefficients, which
  # dominates memory at full scale. Only moment-matched loo needs it; we use
  # waic. Restore it if you ever want loo(moment_match = TRUE).
  stan_model_args = list(
    stanc_options = list("O1"),
    cpp_options = list(STAN_CPP_OPTIMS = TRUE, STAN_NO_RANGE_CHECKS = TRUE)
  ),
  file = file.path(models_dir, sprintf("%s_%s_%d.rds", spec$name, spec$illusion, task_id)),
  file_refit = file_refit
)
wall_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

cat(spec$name, "-", spec$illusion, "shard", task_id, ": SUCCESSFUL.\n")
tryCatch(report_fit(m, paste0(spec$name, "_", spec$illusion), wall_min),
  error = function(e) cat("REPORT failed:", conditionMessage(e), "\n")
)
