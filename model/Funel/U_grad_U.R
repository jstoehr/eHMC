U <- function(theta, data) {
  data$named_param[1, ] <- theta
  y <- as.numeric(data$named_param %>% select(starts_with("y_raw")))
  x <- as.numeric(data$named_param %>% select(starts_with("x_raw")))
  return(0.5 * y * y + 0.5 * sum(x * x))
}

grad_U <- function(theta, data) {
  return(theta)
}


