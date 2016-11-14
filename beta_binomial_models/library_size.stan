data {
  int<lower=0> N; # num individuals
  int<lower=0> P; # num genes
  int<lower=0> gene_counts[N,P]; 
  real<lower=0> concShape; 
  real<lower=0> concRate;  
}
parameters {
  real<lower=0> conc[P]; 
  real<lower=0> library_size[N]; 
  simplex[P] gene_means;
}
model {
  for (n in 1:N) {
    for (p in 1:N) {
      gene_counts[n,p] ~ neg_binomial_2( gene_means[p] * library_size[n], conc[p] ); 
    }
  }
  conc ~ gamma(concShape, concRate);
}
