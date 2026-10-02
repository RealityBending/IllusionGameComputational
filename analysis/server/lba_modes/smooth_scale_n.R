# Does the scale of brms's t2 random-effect basis (Zs) depend on the number of
# rows? If it does, a fixed prior or init on `sds` means a different amount of
# wiggle at different data sizes. Same formula, same covariates, nested subsets.
suppressMessages({library(brms); library(dplyr)})
df <- do.call(rbind, lapply(1:3, function(i) read.csv(sprintf("../../data/illusion_part%d.csv", i))))
df$Illusion_Difference <- abs(df$Illusion_Difference)
df <- mutate(df,
  Illusion_DifferenceZ = 2 * as.numeric(datawizard::normalize(Illusion_Difference)) - 1,
  Illusion_StrengthZ = sign(Illusion_Strength) * as.numeric(datawizard::normalize(abs(Illusion_Strength))),
  .by = "Illusion_Type")
df <- df[df$Illusion_Type == "MullerLyer", ]
f <- bf(RT ~ t2(Illusion_DifferenceZ, Illusion_StrengthZ, k = c(5, 5), bs = c("cr", "cr")))
ppl <- unique(df$Participant)
out <- t(sapply(c(30, 120, 480, length(ppl)), function(n) {
  d <- df[df$Participant %in% ppl[seq_len(n)], ]
  s <- make_standata(f, data = d)
  sc <- sapply(1:3, function(k) sqrt(mean(rowSums(s[[paste0("Zs_1_", k)]]^2))))
  c(participants = n, rows = nrow(d), Zs1 = sc[1], Zs2 = sc[2], Zs3 = sc[3],
    Zs2_x_sqrtN = sc[2] * sqrt(nrow(d)))
}))
print(signif(out, 3))
