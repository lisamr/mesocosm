data{
  int<lower=1> N;
  int<lower=1> T;
  int<lower=0> I[N];
  int<lower=0> n[N];
  vector[N] richness;
  int<lower=1, upper=T> trt[N];
  
  //for simulating predicted probability
  int<lower=1> Nsims;
  int<lower=1, upper=T> trt_sim[Nsims];
  vector[Nsims] richness_sim;
}

parameters{
  vector<lower=0>[2] sigma_pars;
  vector[T] a0;
  vector[T] br;
  real bbar;
  real abar;
  corr_matrix[2] rho;
  real<lower=0> phi;
}

transformed parameters{ 
  vector[N] p;
  vector<lower=0>[N] Alpha; //shape parameters for the beta distribtuion
  vector<lower=0>[N] Beta;
  real theta;
  
  //linear models
   for(i in 1:N){
        p[i] = inv_logit(a0[trt[i]] + br[trt[i]]*richness[i]);
   }
  theta = phi + 2; //want the mass near 2 to get an even spread between 0 and 1
  Alpha = p*theta;
  Beta = (1-p)*theta;
}
model{
  //priors
  //multinormal prior for varying int & slope
  abar ~ normal(0, 1);
  bbar ~ normal(0, 1);
  sigma_pars ~ exponential(1);
  rho ~ lkj_corr(2);
  {
    vector[2] YY[T];
    vector[2] MU;
    MU = [abar, bbar]';
    for ( t in 1:T ) YY[t] = [a0[t] , br[t]]';
    YY ~ multi_normal(MU, quad_form_diag(rho, sigma_pars) );
  }
  phi ~ exponential(.5);

  //likelihood
  target += beta_binomial_lpmf(I | n, Alpha, Beta);
  //I ~ binomial_logit(n, p);
}
generated quantities{
  int<lower=0> y_rep[N];
  //vector[N] log_lik;
  vector[Nsims] p_sim;
  
  for(i in 1:N){
    y_rep[i] = beta_binomial_rng(n[i], Alpha[i], Beta[i]);
    //y_rep[i] = binomial_rng(n[i], inv_logit(p[i]));
    //log_lik[i] = binomial_logit_lpmf(I[i]| n[i], p[i]);
  }
  
  for(i in 1:Nsims){//predictions at the 'average' tray
    p_sim[i] = inv_logit(a0[trt_sim[i]] + br[trt_sim[i]]*richness_sim[i]);
  }
}
