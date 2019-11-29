d <- 100
mu <- rep(0, d)
A <- 0.99^abs(matrix(1:d, d, d) - matrix(1:d, d, d, byrow = T))
data <- list(d = d, mu = mu, A = A, inv_A = solve(A))