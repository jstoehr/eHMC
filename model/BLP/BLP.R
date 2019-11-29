load("German_Credit.RData")
d <- ncol(x)
N <- nrow(x)
data <- list(d = d, N = N, x = x, y = y)