# -------------------- 
# --- MVNORM MODEL
# --------------------

get_data <- function(data_dir = NULL) {
  d <- 100
  rho <- 0.99
  A <- matrix(0, d, d)
  for (i in 1:d) {
    for (j in i:d) {
      A[i , j] <- rho^(j - i)
    }
  }
  A <- A + t(A) - diag(1, d)
  A_inv <- solve(A)
  
  return(list(d = d, A = A, A_inv = A_inv))
}

# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------

get_pars_name <- function(data = get_data()) {
  return(paste0("x[", 1:data$d, "]"))
}