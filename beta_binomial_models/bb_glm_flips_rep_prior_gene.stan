# EAGLE extension for repeated measurements of the same individual, in T conditions. 
# BB(n,a,sigma(phi*beta*x)) is the model, where phi in (-1,1) corresponding to the phase. 
# This version allows for K SNPs
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
  int<lower=0> N; # individuals
  int<lower=0> P; # covariates
  int<lower=0> K; # SNPs
  int<lower=0> T; # time points
  matrix[N,P] x[T]; 
  int<lower=0> ys[N,T,K];
  int<lower=0> ns[N,T,K];
  real<lower=0> concShape; 
  real<lower=0> concRate;
}
parameters {
  real<lower=0> conc[K]; 
  vector[P] beta;
  real<lower=0,upper=1> p[K]; 
}
model {
  vector[N] xb[T]; 
  for (t in 1:T)
    xb[t] <- x[t] * beta;
  for (n in 1:N) {
    for (k in 1:K) {
      vector[T] log_prob_unflipped;
      vector[T] log_prob_flipped; 
      for (t in 1:T) {
        log_prob_unflipped[t] <- beta_binomial_reparam_log(ys[n,t,k], ns[n,t,k], xb[t][n], conc[k]); 
        log_prob_flipped[t] <- beta_binomial_reparam_log(ys[n,t,k], ns[n,t,k], -xb[t][n], conc[k]);
      }
      increment_log_prob( log_sum_exp( log(p[k]) + sum(log_prob_unflipped), log(1.0-p[k]) + sum(log_prob_flipped) ) );
    }
  }
  conc ~ gamma(concShape, concRate);
}
