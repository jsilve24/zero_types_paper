data {
  int<lower=1> N; 
  int<lower=0> N_person;
  int<lower=1> N_batch;
  int<lower=0> y[N];
  int<lower=1, upper=N_person> z[N]; 
  int<lower=1, upper=N_batch> x[N];
  
  // prior hyperparameters
  vector[2] lambda_prior;
}
parameters {
  real<lower=0> eta_id[N_batch-1];
  real<lower=0> lambda[N_person];
  real lambda_mean; 
}
transformed parameters{
  real<lower=0> eta[N_batch];
  eta[1] = 1;
  eta[2:N_batch] = eta_id;
}
model {
  // priors
  lambda_mean ~ normal(-1, lambda_prior[1]);
  lambda ~ lognormal(lambda_mean, lambda_prior[2]);
  eta_id ~ lognormal(0, 2);
  
  // likelihood
  for (i in 1:N)
    y[i] ~ poisson(lambda[z[i]]*eta[x[i]]);
} 
generated quantities {
  real<lower=0> lambda0[N_person];
  for (i in 1:N_person){
    lambda0[i] = lognormal_rng(lambda_prior[1], lambda_prior[2]);  
  }
}
