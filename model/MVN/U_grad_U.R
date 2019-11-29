U <- function(theta, data) {
  return(as.numeric(0.5 * theta %*% (data$inv_A %*% theta)))
}

grad_U <- function(theta, data) {
  return(as.numeric(data$inv_A %*% theta))
}