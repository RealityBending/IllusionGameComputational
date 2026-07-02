library(brms)
library(cmdstanr)
library(rstan)
library(loo)


# Get the number of cores and task ID from the environment variables
options(mc.cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))


setwd("/mnt/lustre/users/psych/oc236/IGComputational/models/combined/")

# List of models
model_names <- c("m_lnr_muller")

# Select the model name based on the task ID
name <- model_names[task_id]

combine_and_save <- function(name) {
  print(paste0(" **", name, ": ", Sys.time()))
  files <- list.files(".", pattern = paste0(name, "_.*rds$"), full.names = TRUE)
  print(paste0("** Found ", length(files), " files for model: ", name))
  m <- brms::combine_models(mlist = lapply(files, readRDS))
  
  # Add criterion and save
  # m <- brms::add_criterion(m, "waic", ndraws = 1500, file = name)
  
  # Save combined model
  saveRDS(m, paste0(name, ".rds"))
  
  print(paste0("** Finished: ", name, " at ", Sys.time()))
}

combine_and_save(name)
print(paste0("** Completed: ", Sys.time()))