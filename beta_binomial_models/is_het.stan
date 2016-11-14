# Straightforward reimplementation of EAGLE using a beta-binomial likelihood. 
data {
  int<lower=0> N; 
  int<lower=0> T; 
  int<lower=0> ys[N,T];
  int<lower=0> ns[N,T];
  real<lower=0> concShape; 
  real<lower=0> concRate;  
  real<lower=0> errorRate; 
}
parameters {
  real<lower=0> conc; 
}
transformed parameters {
  vector[3] probs[N]; 
  for (n in 1:N) {
    probs[n] <- rep_vector(0, 3);
    for (t in 1:T) {
        # likelihood of being het
        probs[n][1] <- probs[n][1] + beta_binomial_log(ys[n,t], ns[n,t], conc * .5, conc * .5);
        # likelihood of being homozygous ref
        probs[n][2] <- probs[n][2] + binomial_log(ys[n,t], ns[n,t], errorRate);
        # likelihood of being homozygous alt
        probs[n][3] <- probs[n][3] + binomial_log(ys[n,t], ns[n,t], 1.0-errorRate);
    }
  }
}
model {
  for (n in 1:N) {
    increment_log_prob(log_sum_exp(probs[n]));
  }
  conc ~ gamma(concShape, concRate);
}
