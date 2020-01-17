// STAN model
data {
  int<lower=0> N;    // number of sample sites 
  vector[N] g;       // process
  matrix[N, N] dist; // distance matrix
}
parameters {
  real mu; 
  real<lower=1.0e-6, upper=2.0> lambda;
  real<lower=0.1, upper=sqrt(100)> eta;
}
transformed parameters {
}
model {
  matrix[N, N] Sigma;
  vector[N] mu_vec;
  
  mu       ~ normal(0,5);
  lambda   ~ uniform(1.0e-6,1.0);
  eta      ~ uniform(0.1,sqrt(100));

  Sigma = eta^2 * exp(-dist/lambda);
  
  for (i in 1:N){
    mu_vec[i] = mu;
  }
  
  g ~ multi_normal(mu_vec, Sigma);
}