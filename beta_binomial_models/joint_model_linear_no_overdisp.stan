data {
  int<lower=0> N; 
  int<lower=0> P; 
  int<lower=0> K; 
  real library_size[N]; 
  matrix[N,P] x_1; 
  matrix[N,P] x_2; 
  int<lower=0> ys[N,K];
  int<lower=0> ns[N,K];
  int<lower=0> gene_counts[N]; 
}
parameters {
  vector[P] beta;
}
model {
  vector[N] xb_1; 
  vector[N] xb_2; 
  xb_1 <- x_1 * beta;
  xb_2 <- x_2 * beta;
  for (n in 1:N) {
    real te; 
    real p; 
    te <- xb_1[n] + xb_2[n]; 
    p <- xb_1[n] / te; 
    for (k in 1:K) {
      ys[n,k] ~ binomial(ns[n,k], p);
    }
    gene_counts[n] ~ poisson( te * library_size[n] ); 
  }
}
