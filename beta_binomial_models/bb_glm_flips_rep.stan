# EAGLE extension for repeated measurements of the same individual, in T=2 conditions. 
# BB(n,a,sigma(phi*beta*x)) is the model, where phi in (-1,1) corresponding to the phase. 
# One SNP version. 
# No constraints on cofficients.
# No learnt prior on mixture weights (symmetric [0.5,0.5] on phi={-1,+1})
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
}
model {
  vector[N] xb[2]; 
  for (k in 1:2)
    xb[k] <- x[k] * beta;
  for (n in 1:N) {
    real p1;
    real p2; 
    p1 <- beta_binomial_reparam_log(ys[n,1], ns[n,1], xb[1][n], conc) + beta_binomial_reparam_log(ys[n,2], ns[n,2], xb[2][n], conc);
    p2 <- beta_binomial_reparam_log(ys[n,1], ns[n,1], -xb[1][n], conc) + beta_binomial_reparam_log(ys[n,2], ns[n,2], -xb[2][n], conc);
    increment_log_prob( log_sum_exp( log(.5) + p1, log(.5) + p2 ) );
  }
  conc ~ gamma(concShape, concRate);
}
