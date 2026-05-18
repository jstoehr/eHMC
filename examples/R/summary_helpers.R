load_run_output <- function(opt) {
  model_dir <- file.path(opt$repo_root, "examples", opt$model)
  algo_dir <- file.path(model_dir, opt$algo)
  
  filename <- file.path(
    algo_dir,
    paste0(get_algoname(opt), get_basename(opt), ".RData")
  )
  
  if (!file.exists(filename)) {
    stop("Run file not found: ", filename)
  }
  
  env <- new.env(parent = emptyenv())
  load(filename, envir = env)
  
  if (!exists("output", envir = env)) {
    stop("Expected object `output` in file: ", filename)
  }
  
  get("output", envir = env)
}

# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------

extract_full_array <- function(output, opt, pars) {
  if (opt$algo == "nuts") {
    return(rstan::extract(output, permuted = FALSE, inc_warmup = FALSE))
  }
  
  out <- array(
    NA_real_,
    dim = c(opt$iter + 1, opt$chains, length(pars) + 1)
  )
  
  for (i in seq_len(opt$chains)) {
    out[, i, ] <- cbind(
      output$sim$samples[[i]],
      output$sim$lp__[[i]]
    )
  }
  
  return(out)
}

# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------

extract_parameter_array <- function(output, opt, pars) {
  if (opt$algo == "nuts") {
    return(rstan::extract(
      output,
      pars = pars,
      permuted = FALSE,
      inc_warmup = FALSE
    ))
  }
  
  out <- array(
    NA_real_,
    dim = c(opt$iter + 1, opt$chains, length(pars))
  )
  
  for (i in seq_len(opt$chains)) {
    out[, i, ] <- output$sim$samples[[i]]
  }
  
  return(out)
}

# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------

extract_run_diagnostics <- function(output, opt) {
  if (opt$algo == "nuts") {
    rstan::get_sampler_params(output, inc_warmup = FALSE)
  } else {
    output$sim$diagnosis
  }
}

# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------

get_run_metric <- function(run_data, key) {
  sapply(run_data, function(x) x[, key])
}

# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------

add_run_metadata <- function(df, opt, total_leapfrog = NA_real_) {
  df %>%
    dplyr::mutate(
      model = opt$model,
      algo = opt$algo,
      seed = opt$seed,
      chains = opt$chains,
      warmup = opt$warmup,
      iter = opt$iter,
      delta = opt$delta,
      metric = opt$metric,
      resampling = ifelse(opt$algo == "nuts", NA_real_, opt$resampling),
      n_train = ifelse(opt$algo == "nuts", NA_real_, opt$n_train),
      ess_train = ifelse(opt$algo == "nuts", NA_real_, opt$ess_train),
      ess_target = ifelse(opt$algo == "nuts", NA_real_, opt$ess_target),
      refresh_rate = ifelse(opt$algo == "prhmc", opt$refresh, NA_real_),
      rand_leap = ifelse(opt$algo == "prhmc", opt$rand, NA)
    )
}

# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------

get_summary <- function(opt, output = NULL, data_dir = NULL) {
  
  if (is.null(output)) {
    output <- load_run_output(opt)
  }
  
  data <- get_data(data_dir)
  pars <- get_pars_name(data)
  n_pars <- length(pars)
  
  # --------------------------------------------------
  # --- COMPUTING PARAMETER SUMMARY STATISTICS
  # --------------------------------------------------
  
  # --- COMPUTING ESS FOR THE FIRST MOMENT
  df <- extract_full_array(output, opt, pars)
  
  ans <- rstan::monitor(
    df,
    warmup = 0,
    probs = c(0.5),
    print = FALSE
  ) %>%
    tibble::as_tibble(rownames = "par")
  
  if (opt$algo != "nuts") {
    ans$par <- c(pars, "lp__")
  }
  
  # --- COMPUTING ESS FOR THE SECOND MOMENT
  df_pars <- extract_parameter_array(output, opt, pars)
  
  temp <- rstan::monitor(
    df_pars^2,
    warmup = 0,
    probs = c(0.5),
    print = FALSE
  )
  
  n_eff_var <- c(temp$n_eff, NA_real_)
  
  # --- COMPUTING ESJD PER PARAMETER
  df_diff <- array(NA_real_, dim = dim(df) - c(1, 0, 0))
  
  for (i in seq_len(opt$chains)) {
    df_diff[, i, ] <- diff(df[, i, ])
  }
  
  df_diff_mat <- matrix(df_diff, ncol = dim(df_diff)[3])
  
  esjd <- colMeans(df_diff_mat^2, na.rm = TRUE)
  diff_lp_mean <- mean(df_diff_mat[, n_pars + 1], na.rm = TRUE)
  
  # --- EXTRACTING DIAGNOSTIC METRICS
  run_data <- extract_run_diagnostics(output, opt)
  
  mh_prob_per_chain <- get_run_metric(run_data, "accept_stat__")
  stepsize_per_chain <- get_run_metric(run_data, "stepsize__")
  n_leapfrog_per_chain <- get_run_metric(run_data, "n_leapfrog__")
  divergences_per_chain <- sapply(run_data, function(x) sum(x[, "divergent__"]))
  
  v_lf <- as.numeric(n_leapfrog_per_chain)
  total_leapfrog <- sum(v_lf)
  
  # ------------------------------ #
  # --- STATISTICS PER PARAMETER
  # ------------------------------ #
  
  param_stats <- ans %>%
    dplyr::mutate(
      n_eff_norm = n_eff / total_leapfrog,
      n_eff_var = n_eff_var,
      n_eff_var_norm = n_eff_var / total_leapfrog,
      ESJD = esjd
    ) %>%
    add_run_metadata(opt)
  
  # ------------------------------ #
  # --- STATISTICS ACCROSS PARAMETERS
  # ------------------------------ #
  
  diagnosis <- tibble::tibble(
    mh_prob_stat = mean(mh_prob_per_chain, na.rm = TRUE),
    total_leapfrog = total_leapfrog,
    ESJD = sum(esjd[seq_len(n_pars)], na.rm = TRUE),
    ESJD_norm = sum(esjd[seq_len(n_pars)], na.rm = TRUE) / total_leapfrog,
    diff_lp__ = diff_lp_mean,
    divergent = sum(divergences_per_chain),
    adapt_accept_stat = ifelse(
      opt$algo == "nuts",
      NA_real_,
      mean(output$tuned_pars$accept_stat__, na.rm = TRUE)
    ),
    freq_accept_stat = ifelse(
      opt$algo == "nuts",
      NA_real_,
      mean(output$sim$freq_accept, na.rm = TRUE)
    )
  ) %>%
    add_run_metadata(opt)

  return(list(param_stats = param_stats, diagnosis = diagnosis))
}
