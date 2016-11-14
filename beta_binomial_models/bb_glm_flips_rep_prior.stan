# EAGLE extension for repeated measurements of the same individual, in T=2 conditions. 
# BB(n,a,sigma(phi*beta*x)) is the model, where phi in (-1,1) corresponding to the phase. 
# One SNP version. 
# No constraints on cofficients.
functions {
  # would be more efficient to pre-calc p
  real beta_binomial_reparam_log(int y, int n, real g, real conc) {
    real p; 
    p <- inv_logit(g);
    return beta_binomial_log(y, n, conc*p, conc*(1.0-p));
  }
}
data {
  int<lower=0> N; 
  int<lower=0> P;
  matrix[N,P] x[2]; 
  int<lower=0> ys[N,2];
  int<lower=0> ns[N,2];
  real<lower=0> concShape; 
  real<lower=0> concRate;  
}
parameters {
  real<lower=0> conc; 
  vector[P] beta;
  real<lower=0,upper=1> p; 
}
model {
  vector[N] xb[2]; 
  for (k in 1:2)
    xb[k] <- x[k] * beta;
  for (n in 1:N) {
    real log_prob_unflipped;
    real log_prob_flipped; 
    log_prob_unflipped <- beta_binomial_reparam_log(ys[n,1], ns[n,1], xb[1][n], conc) + beta_binomial_reparam_log(ys[n,2], ns[n,2], xb[2][n], conc);
    log_prob_flipped <- beta_binomial_reparam_log(ys[n,1], ns[n,1], -xb[1][n], conc) + beta_binomial_reparam_log(ys[n,2], ns[n,2], -xb[2][n], conc);
    increment_log_prob( log_sum_exp( log(p) + log_prob_unflipped, log(1.0-p) + log_prob_flipped ) );
  }
  conc ~ gamma(concShape, concRate);
}
