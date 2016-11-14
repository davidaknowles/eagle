data {
  int<lower=0> N; 
  int<lower=0> P; 
  int<lower=0> K; 
  int<lower=0> T;
  real library_size[N,T]; # set to 0 for missing sample
  matrix[N,P] x_1[T]; 
  matrix[N,P] x_2[T]; 
  int<lower=0> ys[N,T,K];
  int<lower=0> ns[N,T,K];
  int<lower=0> gene_counts[N,T]; 
  real<lower=0> concShape; 
  real<lower=0> concRate;  
}
parameters {
  real<lower=0> conc; 
  vector[P] beta;
  real<lower=0,upper=1> prior[K]; 
}
model {
  vector[N] xb_1[T]; 
  vector[N] xb_2[T]; 
  for (t in 1:T) {
    xb_1[t] <-  x_1[t] * beta ;
    xb_2[t] <-  x_2[t] * beta ;
  }
  for (n in 1:N) {
    vector[T] p; 
    for (t in 1:T) {
      real te; 
      real allele1; 
      real allele2; 
      allele1 = 1e-3 + log1p_exp( xb_1[t][n] ); 
      allele2 = 1e-3 + log1p_exp( xb_2[t][n] ); 
      te <- allele1 + allele2;
      p[t] <- allele1 / te; # cache for next loop
      if (library_size[n,t]>0.0) {
        gene_counts[n,t] ~ neg_binomial_2( te * library_size[n,t], (1.0+gene_counts[n,t]) * conc ); 
      }
    }
    for (k in 1:K) {
      vector[T] log_prob_unflipped;
      vector[T] log_prob_flipped; 
      for (t in 1:T) {
        real conc_here; 
        if (ns[n,t,k] > 0) {
          conc_here <- conc * ns[n,t,k]; 
          log_prob_unflipped[t] <-  beta_binomial_log(ys[n,t,k], ns[n,t,k], conc_here*p[t], conc_here*(1.0-p[t]));
          log_prob_flipped[t] <-  beta_binomial_log(ys[n,t,k], ns[n,t,k], conc_here*(1.0-p[t]), conc_here*p[t]);
        } else {
          log_prob_unflipped[t] <- 0.0;
          log_prob_flipped[t] <- 0.0;
        }
      }
      increment_log_prob( log_sum_exp( log(prior[k]) + sum(log_prob_unflipped), log(1.0-prior[k]) + sum(log_prob_flipped) ) );
    }
  }
  conc ~ gamma(concShape, concRate);
}

