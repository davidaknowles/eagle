data {
  int<lower=0> N; // number of samples
  int<lower=0> K; // number of groups
  int<lower=0> P; // # fixed effects
  matrix[N,P] x; // fixed effect covariates
  int z[N]; // group indices
  int<lower=0> ys[N]; 
  int<lower=0> ns[N];
  real<lower=0> concShape; 
  real<lower=0> concRate;  
}
parameters {
  real<lower=0> conc; 
  real logrev; // log random effect variance
  vector[P] beta;
  real u[K]; // random effects
}
model {
  vector[N] xb; 
  real a[N];
  real b[N];
  real p[N]; 
  real rev;
  rev <- exp(logrev);
  u ~ normal(0, rev);
  xb <- x * beta;
  for (n in 1:N) {
    p[n] <- inv_logit( xb[n] + u[z[n] ]); 
    a[n] <- conc*p[n];
    b[n] <- conc*(1.0-p[n]);
  }
  // beta ~ normal(0,5);
  conc ~ gamma(concShape, concRate);
  ys ~ beta_binomial(ns, a, b);
}
