longest_batch <- function(theta, v, grad_U, eps, L,
                          inv_M, data = NULL, max_length = 2^12) {
  prop <- list(theta = theta, v = v)
  first_turn <- FALSE
  cond <- 0
  l <- L

  for (i in 1:L) {
    prop <- single_leapfrog(prop$theta, prop$v, grad_U, eps, inv_M, data)
    cond <- as.numeric((prop$theta - theta) %*% (inv_M %*% prop$v))
    
    if (is.na(cond)) {
      return(list(
        theta = theta, v = v, l = l,
        first_turn = first_turn, div = TRUE
      ))
    } else if (!first_turn && (cond < 0)) {
      l <- i
      first_turn <- TRUE
    }
  }

  if (!first_turn) {
    while (!first_turn && (l < max_length)) {
      l <- l + 1
      prop <- single_leapfrog(prop$theta, prop$v, grad_U, eps, inv_M, data)
      cond <- as.numeric((prop$theta - theta) %*% (inv_M %*% prop$v))
      
      if (is.na(cond)) {
        return(list(
          theta = theta, v = v, l = l,
          first_turn = first_turn, div = TRUE
        ))
      }
      first_turn <- (cond < 0)
    }
  }

  return(list(
    theta = prop$theta, v = prop$v, l = l,
    first_turn = first_turn, div = FALSE
  ))
}


learn_emp_dist <- function(theta, n_iter, U, grad_U, eps, L0,
                           inv_M, chol_M, data = NULL, max_length = 2^12) {
  current_U <- U(theta, data)
  dim_theta <- length(theta)
  theta <- matrix(theta, n_iter + 1, dim_theta, byrow = TRUE)
  ans <- data.frame(
    rho = rep(0, n_iter), L = rep(L0, n_iter), l = rep(L0, n_iter),
    first_turn = rep(F, n_iter), div = rep(F, n_iter), energy = rep(0, n_iter)
  )

  for (i in 1:n_iter) {
    current_v <- as.numeric(chol_M %*% rnorm(dim_theta))
    current_K <- current_v %*% (inv_M %*% current_v)

    L <- sample(1:L0, 1)
    prop <- longest_batch(
      theta[i, ], current_v, grad_U, eps, L, inv_M,
      data, max_length
    )

    new_U <- U(prop$theta, data)
    new_K <- prop$v %*% (inv_M %*% prop$v)
    rho <- current_U - new_U + 0.5 * (current_K - new_K)

    ans[i, ] <- c(rho, L, prop$l, prop$first_turn, prop$div, 
                  current_U + 0.5 * current_K)
    
    if (is.na(rho) || rho < log(runif(1))) {
      theta[i + 1, ] <- theta[i, ]
    } else {
      theta[i + 1, ] <- prop$theta
      current_U <- new_U
    }
  }
  return(data.frame(theta = theta[-1, ], log_prob = - current_U, ans))
}
