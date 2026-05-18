data {
  int d;
  matrix[d, d] A_inv;
}
parameters {
  vector[d] x;
}
model {
  //x ~ normal(rep_vector(0, d), A);
  target += -0.5 * quad_form(A_inv, x);
}
