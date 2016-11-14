# EAGLE alternative implementation in stan by modeling the variance of allelic expression. 
# One SNP version. 
data {
  int<lower=0> N; 
  int<lower=0> P;
  matrix[N,P] x; 
  int<lower=0> ys[N];
  int<lower=0> ns[N];
  real<lower=0> concShape; 
  real<lower=0> concRate;  
}
parameters {
  real<lower=0> conc; 
  vector[P] beta;
}
model {
  vector[N] xb; 
  xb <- x * beta;
  for (n in 1:N) {
    real p; 
    p <- .5 * exp(xb[n]); 
    ys[n] ~ beta_binomial(ns[n], p , p ); 
  }
  conc ~ gamma(concShape, concRate);
}
