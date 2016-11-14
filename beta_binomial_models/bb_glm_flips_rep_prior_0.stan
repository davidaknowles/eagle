# EAGLE extension for repeated measurements of the same individual, in T=2 conditions in this case. 
# BB(n,a,sigma(phi*beta*x)) is the model, where phi in (-1,0,1) corresponding to the phase with the 
# unobserved cis-SNP (0 if homozygous)
functions {
  # would be more efficient to pre-calc p
  real beta_binomial_reparam_log(int y, int n, real g, real conc) {
    real p; 
    p <- inv_logit(g);
    return beta_binomial_log(y, n, conc*p, conc*(1.0-p));
  }
}
data {
  int<lower=0> N; # Number of individuals
  int<lower=0> P; # Number of covariates
  matrix[N,P] x[2]; # Covariates for allele 1 and 2 (in condition A and B)
  int<lower=0> ys[N,2]; # Alternative allele counts
  int<lower=0> ns[N,2]; # Total allelic counts
  real<lower=0> concShape; # Hyperparameters for concentration
  real<lower=0> concRate;  
}
parameters {
  real<lower=0> conc; # Concentration parameter
  vector[P] beta; # Cofficients
  simplex[3] p; # Mixture weights (for phi=1,-1,0 respectively)
}
model {
  vector[N] xb[2]; # X times beta
  beta ~ normal(0,.01); # Weak regularization on cofficients
  for (k in 1:2)
    xb[k] <- x[k] * beta;
  for (n in 1:N) { # Foreach individual
    vector[3] temp; # Likelihood for each mixture component
    temp[1] <- beta_binomial_reparam_log(ys[n,1], ns[n,1], xb[1][n], conc) + beta_binomial_reparam_log(ys[n,2], ns[n,2], xb[2][n], conc);
    temp[2] <- beta_binomial_reparam_log(ys[n,1], ns[n,1], -xb[1][n], conc) + beta_binomial_reparam_log(ys[n,2], ns[n,2], -xb[2][n], conc);
    temp[3] <- beta_binomial_reparam_log(ys[n,1], ns[n,1], 0, conc) + beta_binomial_reparam_log(ys[n,2], ns[n,2], 0, conc);
    increment_log_prob( log_sum_exp( log(p) + temp ) );
  }
  conc ~ gamma(concShape, concRate);
}
