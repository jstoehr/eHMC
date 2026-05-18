# ----------------- 
# --- BLP MODEL
# -----------------

get_data <- function(data_dir = NULL) {
  ext_data <- file.path(data_dir, "german_data_numeric.csv")
  
  df <- as.matrix(read.table(ext_data))
  N <- nrow(df)
  K <- ncol(df)
  X <- cbind(rep(1, N), df[, 1:(K - 1)])

  return(list(N = N, K = K, X = X, y = 2 - df[, K], sig = 1))
}

# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------

get_pars_name <- function(data = get_data()) {
  return(paste0("theta[", 1:data$K, "]"))
}