data {
  real sigma1;
  real sigma2;
  real lambda;
  real lag;
}
parameters {
  real theta1;
  real theta2;
}
model {
  theta1 ~ normal(0, sigma1);
  theta2 ~ normal(lambda * (theta1^2 - lag), sigma2);
}

