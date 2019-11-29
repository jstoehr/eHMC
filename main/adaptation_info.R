library(rstan)

# --- Main function
extract_inv_M <- function(x) {
  ans <- strsplit(x, "\n")[[1]][4]
  ans <- strsplit(ans, "#")[[1]][2]
  return(as.vector(as.numeric(strsplit(ans, ",")[[1]])))
}

# --- Loading arguments
# !/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# --- Parameters
model_name <- args[1]
M_type <- paste(args[2], "e", sep = "_")
delta <- as.numeric(args[3])
n_calib <- as.numeric(args[4])
# --- File name
exp_name <- paste(model_name, M_type, as.character(100 * delta), sep = "_")

setwd(paste("/home/users/stoehr/eHMC/", model_name, "/", sep = ""))

# --- Loading output of NUTS
load(paste("NUTS", exp_name, ".RData", sep = "_"))

# --- Getting the calibrated stepsize
sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
eps <- sapply(sampler_params, function(x) x[1, "stepsize__"])

# --- Getting L0 for eHMC
# --- Old version
# L0 <- sapply(sampler_params, function(x) median(x[1:warmup, "n_leapfrog__"]) %/% 2)
L0 <- sapply(sampler_params, function(x) x[1:n_calib, "n_leapfrog__"])

# --- Getting the calibrated inverse matrix
inv_M <- sapply(get_adaptation_info(fit), extract_inv_M)

# --- Strating points for eHMC
start <- as.array(fit)[1, , seq_len(nrow(inv_M))]

save(eps, L0, inv_M, start,
  file = paste("adapt_info_NUTS", exp_name, ".RData", sep = "_")
)
