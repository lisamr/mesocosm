data{
  int<lower=1> N; //number of observations
  int<lower=0> K; //number of predictor variables
  matrix[N,K] X; //the model matrix 
  int I[N];//the response variable
  int n[N]; 
  int<lower=0> Nsims;
  matrix[Nsims,K] Xsim;
}
parameters{
  real a0; 
  vector[K] beta; //the regression parameters
}
transformed parameters{ 
  vector[N] p; //mean of the linear model
  p = inv_logit(a0 + X*beta);
}
model{
  a0 ~ normal(0, 1);
  beta ~ normal(0, .5);
  target += binomial_lpmf(I | n, p);
}
generated quantities{
  vector[N] log_lik;
  int<lower=0> y_rep[N];
  vector[Nsims] psim;
  for(i in 1:N) {
    y_rep[i] = binomial_rng(n[i], p[i]);
    log_lik[i] = binomial_lpmf(I[i] | n[i], p[i]);
  } 
  for(i in 1:Nsims){
    psim[i] = inv_logit(a0 + Xsim[i]*beta);
  }
}
