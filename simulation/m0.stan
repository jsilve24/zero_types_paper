data {
  int<lower=1> N; 
  int<lower=0> N_person;
  int<lower=0> y[N];
  int<lower=1, upper=N_person> z[N]; 
  real<lower=0> k; // pseudocount
  
  // prior hyperparameters
  vector[2] lambda_prior;
}
transformed data {
  vector<lower=0>[N] yadj;
  for (i in 1:N)
    yadj[i]= y[i]+k;
}
parameters {
  real<lower=0> lambda[N_person];
}
model {
  // priors
  lambda ~ normal(-1, lambda_prior[1]);
  
  // likelihood
  yadj ~ lognormal(lambda[z],lambda_prior[2]);
} 
generated quantities {
  real<lower=0> lambda0[N_person];
  for (i in 1:N_person){
    lambda0[i] = lognormal_rng(-1, lambda_prior[2]);
  }
}
