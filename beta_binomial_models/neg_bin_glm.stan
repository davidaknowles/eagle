data {
  int<lower=0> N; 
  int<lower=0> P; 
  real library_size[N]; 
  matrix[N,P] x_1; 
  matrix[N,P] x_2; 
  int<lower=0> gene_counts[N]; 
  real<lower=0> concShape; 
  real<lower=0> concRate;  
}
parameters {
  real<lower=0> nb_conc; 
  vector[P] beta;
}
model {
  vector[N] xb_1; 
  vector[N] xb_2; 
  xb_1 <- x_1 * beta;
  xb_2 <- x_2 * beta;
  for (n in 1:N) {
    real te; 
    te <- xb_1[n] + xb_2[n]; 
    gene_counts[n] ~ neg_binomial_2( te * library_size[n], nb_conc ); 
  }
  nb_conc ~ gamma(concShape, concRate);
}
