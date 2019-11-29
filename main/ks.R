# --- For Rstan
library(rstan)
# --- For hand made parallel computation
library(doParallel)
# ---
library(mcmc)
# ---
library(tidyverse)
library(dplyr)
library(tidyselect)

# --- Main function for computing the ecdf in parallel
dist_pi_hat <- function(i, data, x) {
  f <- ecdf(data[, i])
  return(f(x))
}


# --- Loading arguments
# !/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
options(mc.cores = as.numeric(args[1]))

# --- Model
model_name <- args[2]
setwd(paste("/home/users/stoehr/eHMC/", model_name, "/", sep = ""))
# ------ Model expression for Stan
model_stan <- paste(model_name, "stan", sep = ".")
# ------ Data specification
source(paste(model_name, "R", sep = "."))

# --- Running eHMC on the cluster
M_type <- paste(args[3], "e", sep = "_")
n_chain <- as.numeric(args[4])

# --- File name
exp_name <- paste(model_name, M_type, 80, sep = "_")

# --- Empirical distribution function for the ground truth
load(paste("NUTS", exp_name, "long_run_.RData", sep = "_"))
data <- as.array(fit)
bound <- range(data[, , "lp__"])
f_hat_ground <- ecdf(data[, , "lp__"])
rm(data, fit)

# --- Empirical distribution function
ks_nuts <- numeric(n_chain)
ks_ehmc <- numeric(n_chain)
# --- Data for NUTS
load(paste("NUTS", exp_name, ".RData", sep = "_"))
data_nuts <- as.array(fit)
data_nuts <- data_nuts[, , "lp__"]
rm(fit)
# --- Data for eHMC
load_data <- function(i, exp_name) {
  load(paste("eHMC", exp_name, as.character(i), ".RData", sep = "_"))
  return(ehmc_sample[, "log_prob"])
}

# --- Open Cluster
n_cores <- as.numeric(args[1])
cl <- makeCluster(n_cores, type = "FORK")
clusterEvalQ(cl, .libPaths("~/local/R_libs/"))
clusterExport(cl, c("load_data", "dist_pi_hat"))
registerDoParallel(cl)

data_ehmc <- parSapplyLB(cl, seq_len(n_chain), load_data, exp_name)

# --- Computing the values of pi for the chain
bound <- range(c(bound, range(data_nuts), range(data_ehmc)))
x <- seq(bound[1], bound[2], 0.01)
dist_pi_ground <- f_hat_ground(x)
dist_pi_nuts <- parSapplyLB(cl, seq_len(n_chain), dist_pi_hat, data_nuts, x)
dist_pi_ehmc <- parSapplyLB(cl, seq_len(n_chain), dist_pi_hat, data_ehmc, x)

# --- Computing the ks divergence
ks_nuts <- apply(abs(dist_pi_nuts - dist_pi_ground), 2, max)
ks_ehmc <- apply(abs(dist_pi_ehmc - dist_pi_ground), 2, max)

# --- Close Cluster
stopCluster(cl)

# --- Empirical distribution function for ehmc
save(x, dist_pi_ground, dist_pi_nuts, dist_pi_ehmc,
  ks_nuts, ks_ehmc,
  file = paste("ks", exp_name, ".RData", sep = "_")
)
