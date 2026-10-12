# Participants whose own parameters do not mix
#
# In the models that converge (DDM-4, DDM-5, LNR, RDM), the population terms
# agree between chains, but a few dozen of the ~2,200 participants each sit in
# one of two local solutions -- typically a lower ndt against a higher
# boundary or poutlier -- and no chain switches between them. This lists them
# for 2_modelcomparison.qmd (Convergence). They are not excluded: those that
# still do not mix after the longer refits are a result (AGENT.md 4.10).
#
# Only the draws are read, so the laptop's brms is safe here (AGENT.md 3.9).
#
#   Rscript server/mixing_participants.R              # every combined fit in models/
#   Rscript server/mixing_participants.R a.rds b.rds  # these fits
#
# Run from analysis/. Writes server/poorly_mixing_participants.csv.

suppressPackageStartupMessages(library(posterior))

threshold <- 1.1
out <- "server/poorly_mixing_participants.csv"

files <- commandArgs(TRUE)
if (length(files) == 0) {
  # Combined fits only: shards and dated backups carry a further suffix
  files <- list.files("models", full.names = TRUE,
                      pattern = "^gam_[a-z0-9]+(_[a-z]+)?_(MullerLyer|VerticalHorizontal|Ebbinghaus)[.]rds$")
}

rows <- lapply(files, function(f) {
  name <- sub("[.]rds$", "", basename(f))
  m <- readRDS(f)
  d <- subset_draws(as_draws_array(m$fit), variable = "^r_Participant", regex = TRUE)
  rm(m)
  rhat <- vapply(variables(d), \(v) rhat(extract_variable_matrix(d, v)), numeric(1))

  # A per-participant list means something only when the fit as a whole has
  # converged. Where the chains are in different modes (LNR-6, RDM-5, LBA:
  # AGENT.md 4.10) or frozen, nearly every participant is flagged.
  if (mean(rhat > threshold, na.rm = TRUE) > 0.05) {
    message(sprintf("%s: %.0f%% of participant terms above %.2f; not converged as a whole, skipped",
                    name, 100 * mean(rhat > threshold, na.rm = TRUE), threshold))
    return(NULL)
  }

  bad <- rhat[!is.na(rhat) & rhat > threshold]
  if (length(bad) == 0) return(NULL)
  ppt <- sub("^[^[]*[[]([^,]+),.*$", "\\1", names(bad))
  dpar <- sub("^r_Participant_*([^[]*)[[].*$", "\\1", names(bad))
  dpar[dpar == ""] <- "mu"
  res <- do.call(rbind, lapply(split(seq_along(bad), ppt), function(i) {
    i <- i[order(-bad[i])]
    data.frame(Participant = ppt[i[1]], Rhat = max(bad[i]),
               Parameters = paste(dpar[i], collapse = ";"))
  }))
  res$Model <- sub("^gam_([a-z0-9]+).*$", "\\1", name)
  res$Illusion <- sub("^.*_", "", name)
  message(sprintf("%s: %d participants above %.2f", name, nrow(res), threshold))
  res
})

res <- do.call(rbind, rows)
res$Dataset <- sub("_.*$", "", res$Participant)
res <- res[order(-res$Rhat), c("Participant", "Dataset", "Illusion", "Model", "Rhat", "Parameters")]
res$Rhat <- round(res$Rhat, 3)
write.csv(res, out, row.names = FALSE)
message("wrote ", out, ": ", nrow(res), " rows, ", length(unique(res$Participant)), " participants")
