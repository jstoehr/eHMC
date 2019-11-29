U <- function(theta, data, sd = 10) {
  temp <- data$x %*% theta
  return(as.numeric(sum(log(1 + exp(temp))) - data$y %*% temp + (theta %*% theta) / (2 * sd * sd)))
}

grad_U <- function(theta, data, sd = 10) {
  temp <- data$x %*% theta
  temp <- as.numeric(1/(1 + exp(-temp)) - data$y)
  return(as.numeric(theta / (sd * sd) + temp %*% data$x))
}