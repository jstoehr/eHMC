run_nuts <- function(opt, data) {
  stanmodels <- ehmcexamples::get_stanmodels()
  
  run <- rstan::sampling(
    stanmodels[[opt$model]],
    data = data,
    chains = opt$chains,
    iter = opt$iter + opt$warmup,
    warmup = opt$warmup,
    seed = opt$seed,
    algorithm = "NUTS",
    control = list(
      adapt_delta = opt$delta, 
      metric = paste0(opt$metric, "_e")
    )
  )
  
  return(run)
}

# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------

run_ehmc <- function(
    opt,
    log_pdf,
    grad_log_pdf,
    warmup_samples,
    sampling_weights,
    tuned_pars = NULL
) {
  control <- list(
    metric = opt$metric,
    epoch_adapt = opt$epoch,
    adapt_delta = opt$delta,
    stepsize = switch(
      opt$model,
      "banana" = 0.1,
      "BLP" = 0.1,
      "MVNorm" = 0.065,
      0.1
    )
  )
  
  if (opt$algo == "prhmc") {
    control$refresh_rate <- opt$refresh
    control$rand_leapfrog <- opt$rand
  }
  
  run <- ehmc::ehmc(
    opt$iter,
    log_pdf,
    grad_log_pdf,
    warmup_samples,
    sampling_weights = sampling_weights,
    chains = opt$chains,
    control = control,
    tuned_pars = tuned_pars
  )
  
  return(run)
}