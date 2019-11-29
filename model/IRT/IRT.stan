 // Based on stan-dev's code
// https://github.com/stan-dev/stat_comp_benchmarks/tree/master/benchmarks/irt_2pl
data {
  int<lower=0> I;
  int<lower=0> J;
  int<lower=0, upper=1> y[I, J];
}
parameters {
  real log_sigma_theta;
  vector[J] theta;

  real log_sigma_a;
  vector[I] log_a;
  
  real mu_b;
  real log_sigma_b;
  vector[I] b;
}
transformed parameters {
  real<lower=0> sigma_theta = exp(log_sigma_theta);
  real<lower=0> sigma_a = exp(log_sigma_a);
  real<lower=0> sigma_b = exp(log_sigma_b);
  vector[I] a = exp(log_a);
}
model {
  sigma_theta ~ cauchy(0, 2);
  theta ~ normal(0, sigma_theta); 

  sigma_a ~ cauchy(0, 2);
  log_a ~ normal(0, sigma_a);

  mu_b ~ normal(0, 5);
  sigma_b ~ cauchy(0, 2);
  b ~ normal(mu_b, sigma_b);

  for (i in 1:I)
    y[i] ~ bernoulli_logit(a[i] * (theta - b[i]));
    
  target += log_sigma_theta + log_sigma_a + log_sigma_b;
}

