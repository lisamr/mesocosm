data{
  int<lower=1> N; //number of observations
  int<lower=1> S; //number of species
  int<lower=1> J; //number of trays
  int<lower=1> K; //number of predictors
  matrix[N,K] X; //tray-level predictors + intercept
  int<lower=0, upper=1> I[N];//infected or not
  int<lower=1, upper=S> spID[N]; //indices
  int<lower=1, upper=J> trayID[N];
  
  //for plotting predictions
  int<lower=1> Nsims; //number of predicted observations
  matrix[Nsims,K] Xsims; //covariates to predict to 
  int<lower=1> spIDsims[Nsims]; //simulated species
}

parameters{
  real<lower=0> sdj;
  vector[J] zj;
  vector<lower=0>[K] sigma_pars;
  vector[K] B[S];
  vector[K] Bbar;
  real a0;
  corr_matrix[K] rho;
}

transformed parameters{ 
  vector[J] aj;
  vector[N] p;
  //linear models
   aj = sdj*zj; //NC plot intercept
   for(i in 1:N){ //mean model, effects vary by species
     p[i] = X[i]*B[spID[i]] + aj[trayID[i]];
   }
}
model{
  //priors
  zj ~ normal(0,1);
  sdj ~ normal(0, 1);
  Bbar ~ normal(0, 1);
  sigma_pars ~ exponential(1);
  rho ~ lkj_corr(2);
  B ~ multi_normal(Bbar, quad_form_diag(rho, sigma_pars));
  
  //likelihood
  I ~ bernoulli_logit(p);
}
generated quantities{
  int<lower=0, upper=1> y_rep[N];
  vector[N] log_lik;
  vector[Nsims] psims;
  for(i in 1:N){
    y_rep[i] = bernoulli_logit_rng(p[i]);
    log_lik[i] = bernoulli_logit_lpmf(I[i]| p[i]);
  }
  for(i in 1:Nsims){
      psims[i] = inv_logit(Xsims[i]*B[spIDsims[i]]);
  }

}
