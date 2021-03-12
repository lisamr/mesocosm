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
  real<lower=0> phi;
}
transformed parameters{ 
  vector[N] p; //mean of the linear model
  vector<lower=0>[N] Alpha; //shape parameters for the beta distribtuion
  vector<lower=0>[N] Beta;
  real theta;
  p = inv_logit(a0 + X*beta);
  theta = phi + 2; //want the mass near 2 to get an even spread between 0 and 1
  Alpha = p*theta;
  Beta = (1-p)*theta;
}
model{
  a0 ~ normal(0, 1);
  beta ~ normal(0, 1);
  phi ~ exponential(.5);
  //likelihood
  target += beta_binomial_lpmf(I | n, Alpha, Beta);
}
generated quantities{
  vector[N] log_lik;
  int<lower=0> y_rep[N];
  vector[Nsims] psim;
  for(i in 1:N) {
    y_rep[i] = beta_binomial_rng(n[i], Alpha[i], Beta[i]);  
    log_lik[i] = beta_binomial_lpmf(I[i] | n[i], Alpha[i], Beta[i]);
  } 
  for(i in 1:Nsims){
    psim[i] = inv_logit(a0 + Xsim[i]*beta);
  }
}
