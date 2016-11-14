data {
  int<lower=0> N; 
  int<lower=0> P; 
  int<lower=0> K; 
  matrix[N,P] x_1; 
  matrix[N,P] x_2; 
  int<lower=0> ys[N,K];
  int<lower=0> ns[N,K];
  real<lower=0> concShape; 
  real<lower=0> concRate;  
}
parameters {
  real<lower=0> conc[K]; 
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
      ys[n,k] ~ beta_binomial(ns[n,k], conc[k]*p, conc[k]*(1.0-p));
    }
  }
  conc ~ gamma(concShape, concRate);
}
