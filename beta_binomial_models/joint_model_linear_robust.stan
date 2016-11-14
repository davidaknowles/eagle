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
  real<lower=0> concShape; 
  real<lower=0> concRate;  
  vector<lower=0>[3] alpha; # e.g. 1,100,1
  real<lower=0> eps_a; # e.g. 1
  real<lower=0> eps_b; # e.g. 1000
}
parameters {
  real<lower=0> conc[K]; 
  real<lower=0> nb_conc; 
  vector[P] beta;
  simplex[3] pi[K]; 
  real<lower=0,upper=1> eps; # e.g. .001
}
model {
  vector[N] xb_1; 
  vector[N] xb_2; 
  xb_1 <- x_1 * beta;
  xb_2 <- x_2 * beta;
  for (k in 1:K) 
    pi[k] ~ dirichlet(alpha); 
  eps ~ beta(eps_a,eps_b);
  for (n in 1:N) {
    real te; 
    real p; 
    te <- xb_1[n] + xb_2[n]; 
    p <- xb_1[n] / te; 
    for (k in 1:K) {
      real ps[3];
      ps[1] <- log(pi[k][1]) + binomial_log(ys[n,k], ns[n,k], eps); 
      ps[2] <- log(pi[k][2]) + beta_binomial_log(ys[n,k], ns[n,k], conc[k]*p, conc[k]*(1.0-p)); 
      ps[3] <- log(pi[k][3]) + binomial_log(ys[n,k], ns[n,k], 1.0-eps); 
      increment_log_prob(log_sum_exp(ps)); 
    }
    gene_counts[n] ~ neg_binomial_2( te * library_size[n], nb_conc ); 
  }
  conc ~ gamma(concShape, concRate);
  nb_conc ~ gamma(concShape, concRate);
}
