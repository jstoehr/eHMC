diagnosis_single_chain <- function(f_chain) {
  n <- length(f_chain)
  # --- Moments of the chain
  mean_chain <- mean(f_chain)
  var_chain <- var(f_chain)
  # --- Estimating the ESS using the initial convex sequence estimator of Geyer
  temp <- initseq(f_chain)
  lag_cut <- length(temp$Gamma.con)
  gam <- rep(0, n %/% 2 + 1)
  gam[1:lag_cut] <- temp$Gamma.con

  return(data.frame(
    mean = mean_chain,
    var = var_chain,
    ess = n * temp$gamma0 / temp$var.con,
    lag_cut = length(temp$Gamma.con),
    auto_cov_0 = temp$gamma0,
    gam = matrix(gam, nrow = 1)
  ))
}

par_diagnosis_chain <- function(i, f_chain) {
  # --- Useful for multi_D parameter
  return(diagnosis_single_chain(f_chain[, i]))
}

diagnosis_x_chain <- function(df, n_iter_per_chain, n_leapfrog_x_chain) {
  n_chain <- nrow(df)
  # ---Between chain variance
  b <- n_iter_per_chain * var(df$mean)
  # --- Within chain variance
  w <- mean(df$var)
  # --- Variance estimator
  var_hat <- ((n_iter_per_chain - 1) * w + b) / n_iter_per_chain
  # --- Potential scale reduction
  r_hat <- sqrt(var_hat / w)

  # --- Estimating the auto correlation for multiple chains
  cutoff <- max(df$lag_cut)
  gam <- (df %>% select(starts_with("gam")))[, 1:cutoff]
  gam <- 2 * (1 - w / var_hat) + colMeans(gam) / var_hat
  # --- Initial monotone sequence estimator
  test <- which(sign(diff(gam)) == 1)
  if (length(test) != 0) {
    gam <- gam[1:min(test)]
  }
  ess <- n_chain * n_iter_per_chain / (2 * sum(gam) - 1)

  return(data.frame(
    mean = mean(df$mean),
    var_hat = var_hat,
    r_hat = r_hat,
    x_ess = ess,
    x_ess_per_leap = ess / n_leapfrog_x_chain,
    b = b,
    w = w
  ))
}
