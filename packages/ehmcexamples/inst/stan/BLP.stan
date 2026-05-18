data {
  int<lower=0> N; // Number of data points
  int<lower=0> K; // number of predictors with intercept
  matrix[N, K] X; // predictor matrix with intercept
  int<lower=0, upper=1> y[N]; // outcome vector
  real<lower=0> sig; // prior standard deviation
}
parameters {
  vector[K] theta;
}
model {
  theta ~ normal(0, sig);
  y ~ bernoulli_logit(X * theta);
}
