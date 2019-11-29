range(summary_algo[summary_algo[, "algo"] == "NUTS", "x_ess"])
range(summary_algo[summary_algo[, "algo"] == "eHMC", "x_ess"])

max(summary_algo$param)
df <- summary_algo[summary_algo$param == 1, ]
df

df_nuts <- summary_algo[summary_algo[, "algo"] == "NUTS", ]
df_ehmc <- summary_algo[summary_algo[, "algo"] == "eHMC", ]

i <- max(df_nuts$param)
i <- 1
j <- max(df_ehmc$param)
j <- 1
df_nuts[df_nuts$param == i, ]
df_ehmc[df_ehmc$param == j, ]

temp <- details_algo[details_algo == "NUTS", ]
temp <- temp[temp$delta == 0.6, ]
temp <- temp[temp$param == i, ]
df_ehmc[df_ehmc$param == j, ]

boxplot(ks_ehmc, ks_nuts)

boxplot(as.numeric(n_leapfrog_ehmc[1, -1]), as.numeric(n_leapfrog_nuts[1, -1]))

plot(dist_pi_ehmc[,1])
lines(dist_pi_ground, col = "red")
