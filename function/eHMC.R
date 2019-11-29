eHMC <- function(theta, n_iter, U, grad_U, eps, emp_L,
                 inv_M, chol_M, data = NULL) {
  current_U <- U(theta, data)
  dim_theta <- length(theta)
  theta <- matrix(theta, n_iter + 1, dim_theta, byrow = TRUE)
  ans <- data.frame(
    log_prob = rep(0, n_iter), rho = rep(0, n_iter),
    L = rep(1, n_iter), div = rep(F, n_iter),
    energy = rep(0, n_iter), move = rep(0, n_iter)
  )

  for (i in 1:n_iter) {
    current_v <- as.numeric(chol_M %*% rnorm(dim_theta))
    current_K <- current_v %*% (inv_M %*% current_v)
    current_energy <- current_U + 0.5 * current_K

    L <- sample(emp_L, 1)
    prop <- leapfrog(theta[i, ], current_v, grad_U, eps, L, inv_M, data)
    new_U <- U(prop$theta, data)
    new_K <- prop$v %*% (inv_M %*% prop$v)
    new_energy <- new_U + 0.5 * new_K
    rho <- current_energy - new_energy

    if (is.na(rho) || rho < log(runif(1))) {
      theta[i + 1, ] <- theta[i, ]
      move <- 0
    } else {
      theta[i + 1, ] <- prop$theta
      current_U <- new_U
      move <- 1
    }

    ans[i, ] <- c(
      -current_U, rho, L, is.na(rho),
      current_energy, move
    )
  }
  return(data.frame(theta = theta[-1, ], ans))
}
