df <- do.call(rbind, lapply(1:3, function(i) read.csv(sprintf("../../data/illusion_part%d.csv", i))))
df <- df[df$Illusion_Type == "MullerLyer", ]
cat("rows", nrow(df), " participants", length(unique(df$Participant)), "\n")
cat("error rate", mean(df$Error), "\n")
q <- function(x) round(quantile(x, c(.1, .25, .5, .75, .9)), 3)
cat("correct RT quantiles:\n"); print(q(df$RT[df$Error == 0]))
cat("error RT quantiles:\n"); print(q(df$RT[df$Error == 1]))
df$side <- ifelse(df$Illusion_Strength > 0, "conflicting", ifelse(df$Illusion_Strength < 0, "facilitating", "none"))
df$diffbin <- cut(abs(df$Illusion_Difference), quantile(abs(df$Illusion_Difference), 0:4 / 4), include.lowest = TRUE, labels = paste0("D", 1:4))
tab <- aggregate(cbind(Error, RT) ~ side + diffbin, df, mean)
tab$n <- aggregate(Error ~ side + diffbin, df, length)$Error
tab$errRT_med <- aggregate(RT ~ side + diffbin, df[df$Error == 1, ], median)$RT[match(paste(tab$side, tab$diffbin), with(aggregate(RT ~ side + diffbin, df[df$Error == 1, ], median), paste(side, diffbin)))]
tab$corRT_med <- aggregate(RT ~ side + diffbin, df[df$Error == 0, ], median)$RT
print(tab[order(tab$side, tab$diffbin), ], digits = 3)
pe <- aggregate(Error ~ Participant, df, mean)$Error
cat("participant error-rate quantiles:\n"); print(round(quantile(pe, c(0, .05, .25, .5, .75, .95, 1)), 3))
