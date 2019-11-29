data {
  int<lower=1> d;
  int<lower=0> N;
  row_vector[d] x[N];
  int<lower=0,upper=1> y[N];
}
parameters {
  vector[d] theta;
}
model {
  theta ~ normal(0, 10);
  for (n in 1:N)
    y[n] ~ bernoulli_logit((x[n] * theta));
}
