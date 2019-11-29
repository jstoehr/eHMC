leapfrog <- function(theta, v, grad_U, eps, L, inv_M, data = NULL) {
  new_theta <- theta
  new_v <- v - 0.5 * eps * grad_U(new_theta, data)
  for (l in 1:(L - 1)) {
    new_theta <- new_theta + eps * as.numeric(inv_M %*% new_v)
    new_v <- new_v - eps * grad_U(new_theta, data)
  }
  new_theta <- new_theta + eps * as.numeric(inv_M %*% new_v)
  return(list(
    theta = new_theta,
    v = new_v - 0.5 * eps * grad_U(new_theta, data)
  ))
}


single_leapfrog <- function(theta, v, grad_U, eps, inv_M, data) {
  new_v <- v - 0.5 * eps * grad_U(theta, data)
  new_theta <- theta + eps * as.numeric(inv_M %*% new_v)
  new_v <- new_v - 0.5 * eps * grad_U(new_theta, data)
  return(list(theta = new_theta, v = new_v))
}
