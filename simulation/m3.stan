data {
  int<lower=1> N; 
  int<lower=0> N_person;
  int<lower=0> y[N];
  int<lower=1, upper=N_person> z[N]; 
  
  // prior hyperparameters
  vector[2] lambda_prior;
  vector<lower=0>[2] gamma_prior;
}
parameters {
  real<lower=0> lambda[N_person];
  real lambda_mean;
  real<lower=0,upper=1> gamma[N_person];
}
model {
  // priors
  lambda_mean ~ normal(-1, lambda_prior[1]);
  gamma ~ beta(gamma_prior[1],gamma_prior[2]);
  for (i in 1:N_person){
    target += log_sum_exp(log(gamma[i])+normal_lpdf(lambda[i]|0,0.0001), 
                          log1m(gamma[i])+lognormal_lpdf(lambda[i]|lambda_mean,lambda_prior[2]));
  }
  
  // likelihood
  y ~ poisson(lambda[z]);
}
generated quantities {
  real<lower=0> lambda0[N_person];
  { int ind;  // local variable
    int i = 0;                  // local variable
  for (n in 1:N_person){
    ind = bernoulli_rng(gamma[n]);
    if (ind == 0){
      lambda0[n] = lognormal_rng(lambda_mean, lambda_prior[2]);
    } else { 
      // need rejection sampler for truncated distributions in 
      // generated quantities block in stan (ugly but works)
      lambda0[n] = -1;
      while (lambda0[n] < 0 && i < 1000){
        lambda0[n] = normal_rng(0, 0.001);
        i = i + 1;
      } // lambda0 will be negative if can't draw in 1000 tries. 
    }} 
  }
}
