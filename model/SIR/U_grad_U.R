U <- function(theta, data) {
  return(-log_prob(data$fit, theta))
}

grad_U <- function(theta, data) {
  return(-grad_log_prob(data$fit, theta))
}
