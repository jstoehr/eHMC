# --- For manipulating data
library(tidyverse)
# --- For Rstan
library(rstan)
# --- For hand made parallel computation
library(doParallel)
# --- For the monotone positive sequence
library(mcmc)

source("/home/users/stoehr/eHMC/function/ess.R")

# --- Loading arguments
# !/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# --- Model
model_name <- args[2]
setwd(paste("/home/users/stoehr/eHMC/", model_name, "/", sep = ""))
# ------ Model expression for Stan
model_stan <- paste(model_name, "stan", sep = ".")
# ------ Data specification
source(paste(model_name, "R", sep = "."))

# --- Running NUTS on the cluster
options(mc.cores = as.numeric(args[1]))
rstan_options(auto_write = TRUE)

# ------ Parameters
M_type <- paste(args[3], "e", sep = "_")
n_chain <- as.numeric(args[4])
iter <- as.numeric(args[5])
warmup <- as.numeric(args[6])
delta <- as.numeric(args[7])

# ------ File name
exp_name <- paste(model_name, M_type, as.character(100 * delta), sep = "_")
algo_name <- paste("NUTS", exp_name, ".RData", sep = "_")


fit <- stan(model_stan,
  data = data,
  chains = n_chain,
  iter = iter,
  warmup = warmup,
  save_dso = FALSE,
  verbose = FALSE,
  algorithm = "NUTS",
  control = list(
    adapt_engaged = TRUE, stepsize = 0.01, adapt_delta = delta,
    metric = M_type, max_treedepth = 14
  )
)

# --- Saving output
save(fit, file = algo_name, sep = "_")

# --- Computing the number of leapfrog for each chains
sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
n_leapfrog_per_chain_nuts <- sapply(
  sampler_params,
  function(x) sum(x[, "n_leapfrog__"])
)
n_leapfrog_x_chain <- sum(n_leapfrog_per_chain_nuts)

# --- Computing the number of divergence for each chains
div_per_chain_nuts <- sapply(
  sampler_params,
  function(x) sum(x[, "divergent__"])
)

# --- Loading data into an array for computation
data <- as.array(fit)
n_iter <- dim(data)[1]
n_param <- dim(data)[3]
param_names <- c(colnames(data[, 1, -n_param]), "log_prob")

# --- Open Cluster
n_cores <- as.numeric(args[1])
cl <- makeCluster(n_cores, type = "FORK")
clusterEvalQ(cl, .libPaths("~/local/R_libs/"))
clusterExport(cl, c(
  "diagnosis_single_chain", "diagnosis_x_chain",
  "par_diagnosis_chain"
))
registerDoParallel(cl)

# --- Computing the ess
stat_nuts <- parSapply(cl, seq_len(n_chain), par_diagnosis_chain, data[, , 1])
temp <- row.names(stat_nuts)
stat_nuts <- data.frame(t(matrix(unlist(stat_nuts), nrow = nrow(stat_nuts))))
colnames(stat_nuts) <- temp

ess_per_leap <- unlist(
  stat_nuts %>% select(starts_with("ess"))
) / n_leapfrog_per_chain_nuts

x_stat_nuts <- diagnosis_x_chain(stat_nuts, n_iter, n_leapfrog_x_chain)

stat_nuts <- stat_nuts %>% select(
  -starts_with("lag"),
  -starts_with("auto"),
  -starts_with("gam")
)
stat_nuts <- add_column(stat_nuts,
  chain = seq_len(n_chain),
  param = param_names[1], i_param = 1,
  .before = "mean"
)
stat_nuts <- add_column(stat_nuts,
  ess_per_leap = ess_per_leap,
  .after = "ess"
)

for (k in 2:n_param) {
  ans <- parSapply(cl, seq_len(n_chain), par_diagnosis_chain, data[, , k])
  temp <- row.names(ans)
  ans <- data.frame(t(matrix(unlist(ans), nrow = nrow(ans))))
  colnames(ans) <- temp

  ess_per_leap <- unlist(
    ans %>% select(starts_with("ess"))
  ) / n_leapfrog_per_chain_nuts

  x_ans <- diagnosis_x_chain(ans, n_iter, n_leapfrog_x_chain)

  ans <- ans %>% select(-starts_with("lag"), -starts_with("gam"))
  ans <- add_column(ans,
    chain = seq_len(n_chain),
    param = param_names[k], i_param = k,
    .before = "mean"
  )
  ans <- add_column(ans,
    ess_per_leap = ess_per_leap,
    .after = "ess"
  )

  stat_nuts <- rbind(stat_nuts, ans)
  x_stat_nuts <- rbind(x_stat_nuts, x_ans)
}

stat_nuts <- add_column(stat_nuts,
  algo = "NUTS", delta = delta,
  .before = "chain"
)
x_stat_nuts <- add_column(x_stat_nuts,
  algo = "NUTS", delta = delta,
  param = param_names, i_param = seq_len(n_param),
  .before = "mean"
)

# --- Computing statistics on the chains
summary_nuts <- data.frame(
  algo = "NUTS", delta = delta, chain = seq_len(n_chain),
  a_rate = NA, n_leapfrog = n_leapfrog_per_chain_nuts,
  esjd = NA, ks = NA, div = div_per_chain_nuts
)
n_leapfrog_nuts <- matrix(0, n_iter, n_chain)
energy_nuts <- matrix(0, n_iter, n_chain)

for (j in seq_len(n_chain)) {
  ans <- data[, j, -n_param]
  summary_nuts[j, "a_rate"] <- length(unique(ans[, 1])) / n_iter
  summary_nuts[j, "esjd"] <- mean(rowSums((ans[-1, ] - ans[-n_iter, ])^2))
  n_leapfrog_nuts[, j] <- sampler_params[[j]][, "n_leapfrog__"]
  energy_nuts[, j] <- sampler_params[[j]][, "energy__"]
}

log_prob_nuts <- data[, , n_param]

save(log_prob_nuts, n_leapfrog_nuts, energy_nuts,
  file = paste("log_prop", algo_name, ".RData", sep = "_")
)
save(stat_nuts, x_stat_nuts, summary_nuts,
  file = paste("result", algo_name, ".RData", sep = "_")
)
# --- Close Cluster
stopCluster(cl)
