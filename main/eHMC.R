# --- For manipulating data
library(tidyverse)
# --- For Rstan
library(rstan)
# --- For hand made parallel computation
library(doParallel)
# --- For the monotone positive sequence
library(mcmc)

source("/home/users/stoehr/eHMC/function/leapfrog.R")
source("/home/users/stoehr/eHMC/function/empirical_dist_L.R")
source("/home/users/stoehr/eHMC/function/ess.R")
source("/home/users/stoehr/eHMC/function/eHMC.R")

# --- Definition of the function for parallel computation
par_eHMC <- function(i, start, iter, U, grad_U, eps, L0,
                     inv_M, data, algo_name,
                     model_stan, do_fit = TRUE) {
  if (do_fit) {
    fit <- stan(model_stan, data = data, chains = 0)
    data <- list_modify(data, fit = fit)
  }

  M_i <- diag(1 / as.numeric(inv_M[, i]))
  inv_M_i <- diag(as.numeric(inv_M[, i]))
  chol_M_i <- sqrt(M_i)
  emp_L_i <- L0[, i]

  # --- Old Version
  # burn_in <- learn_emp_dist(start[i, ], warmup,
  #                           U, grad_U,
  #                           as.numeric(eps[i]),
  #                           L0[i],
  #                           inv_M_i,
  #                           chol_M_i,
  #                           data = data
  # )
  # time_burn <- proc.time() - time_elapsed
  #
  # emp_L <- burn_in[, "l"]
  # emp_L <- emp_L[which(!burn_in[, "div"])]

  time_elapsed <- proc.time()
  ehmc_sample <- eHMC(start[i, ], iter,
    U, grad_U,
    as.numeric(eps[i]),
    emp_L_i,
    inv_M_i,
    chol_M_i,
    data = data
  )
  time_elapsed <- proc.time() - time_elapsed

  save(ehmc_sample, time_elapsed,
    file = paste(algo_name, i, ".RData", sep = "_")
  )

  log_prob_ehmc <- ehmc_sample %>% select(
    starts_with("log"), starts_with("L"), starts_with("energy")
  )

  # --- Computing statistics on the chain
  # ------ Average accepted moves
  acceptance_rate <- mean(ehmc_sample$move)
  # ------ Total number of leapfrog
  n_leapfrog <- sum(ehmc_sample$L)
  # ------ Total number of divergence
  div <- sum(ehmc_sample$div)
  # ------ Computing ess and esjd
  df <- ehmc_sample %>% select(starts_with("theta"), starts_with("log_p"))
  # --------- esjd
  esjd <- mean(rowSums((df[-1, -ncol(df)] - df[-iter, -ncol(df)])^2))
  # --------- ess
  ans <- sapply(seq_len(ncol(df)), par_diagnosis_chain, df)
  temp <- row.names(ans)
  ans <- data.frame(t(matrix(unlist(ans), nrow = nrow(ans))))
  colnames(ans) <- temp
  temp <- unlist(ans %>% select(starts_with("ess"))) / n_leapfrog
  ans <- add_column(ans, ess_per_leap = temp, .after = "ess")

  param_names <- c(colnames(data$named_param), "log_prob")

  return(list(
    stat_ehmc = data.frame(
      chain = i, param = param_names,
      i_param = seq_len(ncol(df)), ans
    ),
    summary_ehmc = data.frame(
      chain = i, a_rate = acceptance_rate,
      n_leapfrog = n_leapfrog, esjd = esjd, ks = NA, div = div
    ),
    log_prob_ehmc = data.frame(chain = i, log_prob_ehmc)
  ))
}

x_chain_ess <- function(k, df, n_iter_per_chain, n_leapfrog_x_chain) {
  ans <- df[df[, "i_param"] == k, ]
  return(data.frame(
    i_param = k,
    diagnosis_x_chain(ans, n_iter_per_chain, n_leapfrog_x_chain)
  ))
}



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
source("U_grad_U.R")

# --- Loading eHMC settings
options(mc.cores = as.numeric(args[1]))
rstan_options(auto_write = TRUE)

M_type <- paste(args[3], "e", sep = "_")
first_chain <- as.numeric(args[4])
last_chain <- as.numeric(args[5])
n_chain <- last_chain - first_chain + 1
iter <- as.numeric(args[6])
delta <- as.numeric(args[7])

# --- File name
exp_name <- paste(model_name, M_type, as.character(100 * delta), sep = "_")
algo_name <- paste("eHMC", exp_name, sep = "_")

# --- Fit from stan required
do_fit <- TRUE
if (model_name %in% c("BLP", "Funel", "MVN")) {
  do_fit <- FALSE
}

# --- Loading samplers parameters from NUTS warmup
load(paste("adapt_info_NUTS", exp_name, ".RData", sep = "_"))

n_param <- ncol(start)
named_param <- data.frame(matrix(rep(NA, n_param), 1))
colnames(named_param) <- colnames(start)
data <- list_modify(data, named_param = named_param)

# --- Open Cluster
n_cores <- as.numeric(args[1])
cl <- makeCluster(n_cores, type = "FORK")
clusterSetRNGStream(cl, as.integer(args[8]))
clusterEvalQ(cl, .libPaths("~/local/R_libs/"))
clusterExport(cl, c(
  "U", "grad_U", "leapfrog", "single_leapfrog",
  "diagnosis_single_chain", "par_diagnosis_chain",
  "eHMC"
))
registerDoParallel(cl)

# --- Running eHMC on the different cores
result_ehmc <- parLapplyLB(
  cl, first_chain:last_chain, par_eHMC, start,
  iter, U, grad_U, eps, L0, inv_M,
  data, algo_name, model_stan, do_fit
)

# --- Combining the results
log_prob_ehmc <- matrix(0, iter, n_chain)
n_leapfrog_ehmc <- matrix(0, iter, n_chain)
energy_ehmc <- matrix(0, iter, n_chain)

stat_ehmc <- result_ehmc[[1]]$stat_ehmc
summary_ehmc <- result_ehmc[[1]]$summary_ehmc
temp <- result_ehmc[[1]]$log_prob_ehmc
log_prob_ehmc[, 1] <- temp$log_prob
n_leapfrog_ehmc[, 1] <- temp$L
energy_ehmc[, 1] <- temp$energy


for (i in 2:length(result_ehmc)) {
  stat_ehmc <- rbind(
    stat_ehmc,
    result_ehmc[[i]]$stat_ehmc
  )
  
  summary_ehmc <- rbind(
    summary_ehmc,
    result_ehmc[[i]]$summary_ehmc
  )

  temp <- result_ehmc[[i]]$log_prob_ehmc
  log_prob_ehmc[, i] <- temp$log_prob
  n_leapfrog_ehmc[, i] <- temp$L
  energy_ehmc[, i] <- temp$energy
}

# --- Computing the x chain results
n_leapfrog_x_chain <- sum(summary_ehmc$n_leapfrog)

x_result_ehmc <- parLapplyLB(
  cl, seq_len(n_param + 1), x_chain_ess,
  stat_ehmc, iter, n_leapfrog_x_chain
)

x_stat_ehmc <- x_result_ehmc[[1]]
for (i in (seq_len(n_param) + 1)) {
  x_stat_ehmc <- rbind(x_stat_ehmc, x_result_ehmc[[i]])
}

# --- Saving the output
# ------ Output with the stat for each chain and each parameter
stat_ehmc <- stat_ehmc %>% select(
  -starts_with("lag"),
  -starts_with("auto"),
  -starts_with("gam")
)

stat_ehmc <- add_column(stat_ehmc,
  algo = "eHMC", delta = delta, 
  .before = "chain"
)

# ------ Output with the stat for each chain
summary_ehmc <- add_column(summary_ehmc,
  algo = "eHMC", delta = delta, 
  .before = "chain"
)

# ------ Output with the stat for each parameter accross the chains
x_stat_ehmc <- add_column(x_stat_ehmc,
  algo = "eHMC", delta = delta,
  param = c(colnames(start), "log_prob"),
  .before = "i_param"
)

save(log_prob_ehmc, n_leapfrog_ehmc, energy_ehmc,
     file = paste("log_prob", algo_name, ".RData", sep = "_"))

save(stat_ehmc, x_stat_ehmc, summary_ehmc, 
  file = paste("result", algo_name, ".RData", sep = "_")
)

# --- Close Cluster
stopCluster(cl)
