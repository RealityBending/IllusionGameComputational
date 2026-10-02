# The priors gam_lba was actually fitted with, as brms resolved them per parameter.
options(width = 200)
suppressMessages(library(brms))
m <- readRDS("../models/gam_lba_MullerLyer.rds")
ps <- prior_summary(m, all = TRUE)
ps <- ps[ps$class %in% c("Intercept", "b", "sd", "sds"), c("prior", "class", "coef", "group", "dpar", "source")]
ps$coef <- sub("t2Illusion_DifferenceZIllusion_StrengthZ", "t2", ps$coef)
print(as.data.frame(ps), row.names = FALSE)
cat("\ninits: cogmod_inits(f, data), default jitter (0.25 population, 0.05 hierarchical)\n")
str(m$stan_args$control)
