# =========================================================================
# Posterior predictions of ONE combined model, for 1_modelcomparison.qmd
# =========================================================================
# Reads IGC_MODELS_DIR/combined/<model>_<illusion>.rds (written by
# combine_model.R) and writes IGC_MODELS_DIR/predictions/<model>_<illusion>.rds:
# ~20 MB of tidy results -- efficiency, loo, predictive check, parameter
# curves and heatmaps -- that the qmd relabels and plots without ever loading
# the fit. What is computed, and how, is predictions.R; this only runs it.
# Which model is IGC_MODEL. Submitted by `./hpc predict <model>`.
#
#   IGC_PRED_<SETTING>   overrides one of igc_prediction_settings(), by its
#                        name in capitals, e.g. IGC_PRED_PPC_NDRAWS=1000 or
#                        IGC_PRED_CURVES_ITERATIONS=2000. Numbers only.
#
# Re-run it after anything in predictions.R changes the output (the grid, the
# settings), and after re-combining the model: the qmd warns when the fit a
# prediction file was made from is not the one it has locally, but only if it
# has that fit locally.

library(brms)
# cogmod provides the families' posterior_predict / posterior_epred methods
library(cogmod)

source("models.R")
source("predictions.R")

spec <- igc_model(Sys.getenv("IGC_MODEL", unset = ""))
name <- paste0(spec$name, "_", spec$illusion)

models_dir <- Sys.getenv("IGC_MODELS_DIR", unset = "models")
in_file <- file.path(models_dir, "combined", paste0(name, ".rds"))
out_dir <- file.path(models_dir, "predictions")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, paste0(name, ".rds"))

cat("**", name, ":", format(Sys.time()), "\n")
if (!file.exists(in_file)) stop("no combined fit at ", in_file, " -- run ./hpc combine first")

# Settings: the defaults, the cores SLURM gave, and any IGC_PRED_* override
settings <- list(cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
for (k in names(igc_prediction_settings())) {
  v <- Sys.getenv(paste0("IGC_PRED_", toupper(k)), unset = "")
  if (nzchar(v)) settings[[k]] <- as.numeric(v)
}
str(do.call(igc_prediction_settings, settings))

m <- readRDS(in_file)
attr(m, "illusion") <- spec$illusion
cat("** read", in_file, "with", brms::ndraws(m), "draws\n")

t0 <- Sys.time()
res <- run_predictions(m, name, fit_file = in_file, settings = settings)
saveRDS(res, out)
cat(sprintf("** wrote %s (%.1f MB) in %.1f min\n", out, file.size(out) / 1e6,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("** finished:", name, "at", format(Sys.time()), "\n")
