data {
  int<lower=1> N; 
  int<lower=1> N_batch;
  int<lower=1> N_person;
  int<lower=0> y[N];
  int<lower=1, upper=N_person> z[N];
  int<lower=1, upper=N_batch> x[N];
  
  // prior hyperparameters
  vector[2] lambda_prior;
  vector<lower=0>[2] theta_prior;
}
parameters {
  real<lower=0> lambda[N_person];
  real lambda_mean;
  real<lower=0,upper=1> theta[N_batch];
}
model {
  // priors
  lambda_mean ~ normal(-1, lambda_prior[1]);
  lambda ~ lognormal(lambda_mean, lambda_prior[2]);
  theta ~ beta(theta_prior[1],theta_prior[2]);
  
  // likelihood
  for (i in 1:N){
    if (y[i] == 0){
      target += log_sum_exp(log(theta[x[i]]), log1m(theta[x[i]])+poisson_lpmf(0|lambda[z[i]]));
    } else {
      target += log1m(theta[x[i]])+poisson_lpmf(y[i]|lambda[z[i]]);
    }
  }
}
generated quantities {
  real<lower=0> lambda0[N_person];
  for (i in 1:N_person){
    lambda0[i] = lognormal_rng(lambda_prior[1], lambda_prior[2]);  
  }
}
