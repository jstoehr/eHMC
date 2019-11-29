U <- function(theta, data) {
  mu <- theta[1]
  alpha <- theta[2]
  gamma <- theta[3]
  h_std <- theta[-(1:3)]
  phi <- (exp(alpha) - 1) / (exp(alpha) + 1)
  sigma <- exp(gamma)
  h <- h_std * sigma
  h[1] <- h[1] / sqrt(1 - phi * phi)
  h <- h + mu
  for (t in 2:data$T) {
    h[t] <- h[t] + phi * (h[t - 1] - mu)
  }

  if (alpha > 709 || gamma > 709 || abs(phi) >= 1 ||
    sum(h > 1418) > 0 || sum(is.infinite(h)) > 0 ||
    sum(is.nan(theta) + is.na(theta)) > 0) {
    return(NA)
  } else {
    # save(theta, file = data$filename_U)
    return(-log_prob(data$fit, theta))
  }
}

grad_U <- function(theta, data) {
  mu <- theta[1]
  alpha <- theta[2]
  gamma <- theta[3]
  h_std <- theta[-(1:3)]
  phi <- (exp(alpha) - 1) / (exp(alpha) + 1)
  sigma <- exp(gamma)
  h <- h_std * sigma
  h[1] <- h[1] / sqrt(1 - phi * phi)
  h <- h + mu
  for (t in 2:data$T) {
    h[t] <- h[t] + phi * (h[t - 1] - mu)
  }
  
  if (alpha > 709 || gamma > 709 || abs(phi) >= 1 ||
      sum(h > 1418) > 0 || sum(is.infinite(h)) > 0 ||
      sum(is.nan(theta) + is.na(theta)) > 0) {
    return(NA)
  } else {
    # save(theta, file = data$filename_grad_U)
    return(-grad_log_prob(data$fit, theta))
  }
}
