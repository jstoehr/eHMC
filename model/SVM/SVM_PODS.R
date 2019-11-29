set.seed(678)
T.0 <- 500
mu <- -1.02
phi <- 0.95
sigma <- 0.25
x <- rnorm(T.0, mu, sigma)
x[1] <- x[1] / sqrt(1 - phi * phi)
for(t in 2:T.0) {
  x[t] <- x[t] + phi * (x[t - 1] - mu)
}
y <- rnorm(T.0, sd = exp(x / 2))

save.image(file = "SVM_PODS.RData")
