# EAGLE alternative implementation in stan using flip variables (which are marginalized out)
# BB(n,a,sigma(phi*beta*x)) is the model, where phi in (-1,1) corresponding to the phase. 
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
    real a;
    real b;
    real p; 
    p <- inv_logit(xb[n]); 
    a <- conc*p;
    b <- conc*(1.0-p);
    increment_log_prob( log(.5) + beta_binomial_log(ys[n], ns[n], a, b), log(.5) + beta_binomial_log(ys[n], ns[n], b, a) );
  }
  conc ~ gamma(concShape, concRate);
}
