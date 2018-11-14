data {
  int<lower=1> N; 
  int<lower=0> N_person;
  int<lower=0> y[N];
  int<lower=1, upper=N_person> z[N]; 
  
  // prior hyperparameters
  vector[2] lambda_prior;
}
parameters {
  real<lower=0> lambda[N_person];
  real lambda_mean; 
}
model {
  // priors
  lambda_mean ~ normal(-1, lambda_prior[1]);
  lambda ~ lognormal(lambda_mean, lambda_prior[2]);
  
  // likelihood
  y ~ poisson(lambda[z]);
} 
generated quantities {
  real<lower=0> lambda0[N_person];
  for (i in 1:N_person){
    lambda0[i] = lognormal_rng(lambda_prior[1], lambda_prior[2]);  
  }
}
