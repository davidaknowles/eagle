# EAGLE alternative implementation in stan using flip variables (which are marginalized out)
# BB(n,a,sigma(phi*beta*x)) is the model, where phi in (-1,1) corresponding to the phase. 
# Data can additionally come from a outlier model to account for incorrect genotyping/imprinting. 
# One SNP version. 
data {
  int<lower=0> N; 
  int<lower=0> P;
  matrix[N,P] x; 
  int<lower=0> ys[N];
  int<lower=0> ns[N];
  real<lower=0> concShape; 
  real<lower=0> concRate;
  vector<lower=0>[3] alpha; # e.g. 1,100,1
  real<lower=0> eps_a; # e.g. 1
  real<lower=0> eps_b; # e.g. 1000
}
parameters {
  real<lower=0> conc; 
  vector[P] beta;
  simplex[3] pi; 
  real<lower=0,upper=1> eps; # e.g. .001
}
model {
  vector[N] xb; 
  pi ~ dirichlet(alpha); 
  eps ~ beta(eps_a,eps_b);
  xb <- x * beta;
  for (n in 1:N) {
    real ps[3];
    real p; 
    p <- inv_logit(xb[n]); 
    ps[1] <- log(pi[1]) + binomial_log(ys[n], ns[n], eps); 
    ps[2] <- log(pi[2]) + beta_binomial_log(ys[n], ns[n], conc*p, conc*(1.0-p)); 
    ps[3] <- log(pi[3]) + binomial_log(ys[n], ns[n], 1.0-eps); 
    increment_log_prob(log_sum_exp(ps)); 
  }
  // beta ~ normal(0,5);
  conc ~ gamma(concShape, concRate);
  
}
