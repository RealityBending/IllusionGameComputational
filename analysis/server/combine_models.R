library(brms)
library(cmdstanr)
library(rstan)
library(loo)


# Get the number of cores and task ID from the environment variables
options(mc.cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))


setwd("/mnt/lustre/users/psych/oc236/IGComputational/models/combined/")

# List of models
model_names <- c("gam_lnr_muller")

# Select the model name based on the task ID
name <- model_names[task_id]

combine_and_save <- function(name) {
  print(paste0(" **", name, ": ", Sys.time()))
  files <- list.files(".", pattern = paste0(name, "_.*rds$"), full.names = TRUE)
  print(paste0("** Found ", length(files), " files for model: ", name))

  # Full
  m <- brms::combine_models(mlist = lapply(files, readRDS))
  m <- brms::add_criterion(m, "waic", ndraws = 1500, file = name) # waic is faster than loo
  saveRDS(m, paste0(name, ".rds"))

  # Mini
  mini <- brms::combine_models(mlist = lapply(files[1:2], readRDS))
  mini <- brms::add_criterion(mini, "waic", ndraws = 500, file = name) # waic is faster than loo
  saveRDS(mini, paste0(name, "_mini.rds"))

  print(paste0("** Finished: ", name, " at ", Sys.time()))
}

combine_and_save(name)
print(paste0("** Completed: ", Sys.time()))
