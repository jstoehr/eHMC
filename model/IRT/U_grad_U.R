U <- function(theta, data) {
  if(sum(is.nan(theta) + is.na(theta)) > 0) {
    return(NA)
  } else {
    return(-log_prob(data$fit, theta))
  }
}

grad_U <- function(theta, data) {
  if(sum(is.nan(theta) + is.na(theta)) > 0) {
    return(rep(NA, length(theta)))
  } else {
    return(-grad_log_prob(data$fit, theta))
  }
}
