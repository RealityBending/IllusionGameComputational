# =========================================================================
# Priors against posteriors for every combined fit (priors.R)
# =========================================================================
# Reads IGC_MODELS_DIR/combined/<name>.rds and writes
# IGC_MODELS_DIR/priors/<name>.rds (a few hundred kB each), for 1_models.qmd.
# Fits with a date suffix (_YYYYMMDD, superseded) are left out. A fit whose
# file has not changed size since its prior file was written is skipped;
# IGC_PRIORS_FORCE=1 redoes them all. IGC_MODEL, if set and not 'all', limits
# it to the fits whose name starts with it. Submitted by `./hpc priors`.

library(brms)
source("priors.R")

models_dir <- Sys.getenv("IGC_MODELS_DIR", unset = "models")
want <- Sys.getenv("IGC_MODEL", unset = "all")
force <- Sys.getenv("IGC_PRIORS_FORCE", unset = "") == "1"
out_dir <- file.path(models_dir, "priors")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(file.path(models_dir, "combined"), pattern = "\\.rds$", full.names = TRUE)
files <- files[!grepl("_[0-9]{8}\\.rds$", files)]
if (nzchar(want) && want != "all") files <- files[startsWith(basename(files), paste0(want, "_"))]

for (f in files) {
  name <- sub("\\.rds$", "", basename(f))
  out <- file.path(out_dir, basename(f))
  if (!force && file.exists(out)) {
    old <- readRDS(out)$meta
    if (identical(old$version, igc_priors_version) && identical(old$fit_size, file.size(f))) {
      cat("** up to date, skipped:", name, "\n")
      next
    }
  }
  t0 <- Sys.time()
  res <- tryCatch({
    m <- readRDS(f)
    r <- igc_prior_check(m, name, fit_file = f)
    rm(m); gc()
    r
  }, error = function(e) { cat("!! ", name, ": ", conditionMessage(e), "\n", sep = ""); NULL })
  if (is.null(res)) next
  saveRDS(res, out)
  tab <- res$table
  cat(sprintf("** %s (cogmod %s): %d params, lprior max |diff| %.2g, %d conflict, %d tail, %d weak data, %.1f min\n",
              name, res$meta$cogmod, nrow(tab), res$meta$lprior_maxdiff, sum(tab$Flag == "conflict"),
              sum(tab$Flag == "tail"), sum(tab$Flag == "weak data"),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat("** finished at", format(Sys.time()), "\n")
