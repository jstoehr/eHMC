data {
  int<lower=0> T;
  vector[T] y;
}
parameters {
  real mu;
  real alpha; // Parameter corresponding to the change of parametrisation for phi
  real gamma; // Parameter corresponding to the change of parametrisation for sigma
  vector[T] h_std;
}
transformed parameters {
  real<lower=-1,upper=1> phi = (exp(alpha) - 1) / (exp(alpha) + 1);
  real<lower=0> sigma = exp(gamma);
  vector[T] h = h_std * sigma;
  h[1] /= sqrt(1 - phi * phi);
  h += mu;
  for (t in 2:T)
    h[t] += phi * (h[t-1] - mu);
}
model {
  phi ~ uniform(-1, 1);
  sigma ~ cauchy(0, 5);
  mu ~ cauchy(0, 10);
  h_std ~ std_normal();
  y ~ normal(0, exp(h / 2));
  target += alpha - 2 * log(exp(alpha) + 1); // Jocobian term due to the reparametrisation
  target += gamma;
}

