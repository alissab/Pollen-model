// STAN model
data {
  int<lower=0> N;       // number of sample sites 
  int<lower=0> N_taxa;  // number of taxa 
  int y[N,N_taxa];      // process
  matrix[N, N] dist;    // distance matrix
}
parameters {
  vector[N_taxa] mu; 
  vector<lower=1.0e-6, upper=2.0>[N_taxa] lambda;
  vector<lower=0.1, upper=sqrt(100)>[N_taxa] eta;
  vector[N] g[N_taxa];
}
transformed parameters {
}
model {
  matrix[N, N] Sigma[N_taxa];
  vector[N] mu_vec[N_taxa];
  vector[N] sum_exp_g;
  vector[N_taxa] r[N];
  
  for (n in 1:N_taxa){
    mu[n]       ~ normal(0,5);
    lambda[n]   ~ uniform(1.0e-6,2.0);
    eta[n]      ~ uniform(0.1,sqrt(100));
    
    Sigma[n] = eta[n]^2 * exp(-dist/lambda[n]);
    
    for (i in 1:N){
      mu_vec[n][i] = mu[n];
    }
    
    g[n] ~ multi_normal(mu_vec[n], Sigma[n]);
  }
  
  // sum process vals for each i
  for (i in 1:N) {
    sum_exp_g[i] = 0.0;
    for (n in 1:N_taxa)
      sum_exp_g[i] = sum_exp_g[i] + exp(g[n,i]);
  }

  for (i in 1:N) {
    for (n in 1:N_taxa){
      r[i,n] = exp(g[n,i])/sum_exp_g[i];
    }
  }

  // link composition vector to count data through multinomial
  for (i in 1:N) 
    y[i] ~ multinomial(r[i]);
}
