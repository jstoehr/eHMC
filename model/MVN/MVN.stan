data {
  int<lower=0> d;
  matrix[d, d] A;
  vector[d] mu;
}
parameters {
  vector[d] theta;
}
model {
  target += multi_normal_lpdf(theta | mu, A);
}
