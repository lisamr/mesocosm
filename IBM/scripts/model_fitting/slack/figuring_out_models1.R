library(rethinking)
library(tidyverse)
rm(list = ls())

#load data----

#try to recover parameters 
#beta = .01
#alpha = .1

#aggregated data
outAgg <- readRDS("IBM/outputs/sampledata/nonreplicated.RDS")
#not aggregated data. also, 2 replicates in an array.
#dim = [inds, time, replicate]
outInd <- readRDS('IBM/outputs/sampledata/replicated.RDS')

#define some stuff 
nind = dim(outInd)[1]
nreps = dim(outInd)[3]
tfinal = dim(outInd)[2]




#model non-replicated data with Binomial likelihood----

# format data for STAN
dat_Agg <- outAgg %>% 
  pivot_wider(names_from = state, values_from = n, values_fill = list(n=0)) %>% 
  dplyr::select(-time) %>% 
  as.list %>% 
  c(t = tfinal) %>% list_modify(N=168)
str(dat_Agg)

binom = 
  "
data { 
  //declare variables
  int<lower=0> N; //number of individuals in each tray
  int<lower=0> t; //number of time points
  //states
  int<lower = 0, upper =N> C[t];
  int<lower = 0, upper =N> S[t]; 
  int<lower = 0, upper =N> I[t];
}
transformed data { 
  //Calculate transitions, based on change in states
  int<lower=0, upper = N> z_ci[t-1]; 
  int<lower=0, upper = N> z_si[t-1];
  
  //I think it should be for(i in 2:t) z_ci[2] = C[1] - C[2]?
    for(i in 1:(t-1)){
    z_ci[i] = C[i] - C[i+1]; //new infections of C
    z_si[i] = S[i] - S[i+1]; //new infections of S
    }
}
parameters {
  real<lower=0, upper = 1> alpha; //I think it's also between 0 and 1.
  real<lower=0, upper = 1> beta;
}

model {
  real p_ci; //prob. of C -> I
  real p_si; //prob. of S -> I
  
  //priors
  alpha ~ beta(1,5); //previously had at gamma(1,1)
  beta ~ beta(1,5);

 //generate number of infecteds from C and S at each time
 //C, S, and I are the number of individuals in that state at time t
 for(i in 1:(t-1)){
    p_ci = 1 - exp(-(alpha + beta*I[i]));
    z_ci[i] ~ binomial(C[i], p_ci); 
  //only define z_si when there are infections. distribution is degenerate and STAN has trouble
    if(I[i] > 0){
      p_si = 1 - exp(-(beta*I[i]));
      z_si[i] ~ binomial(S[i], p_si); 
    }
 }
}
"
binom_mod = stan_model(model_code = binom)
fit_binom = sampling(binom_mod,dat_Agg, iter = 2000, chains = 3)

#check it out
precis(fit_binom) %>% plot #covers the parameters. Interval of alpha is big. 
print(fit_binom)
traceplot(fit_binom)





#model non-replicated data with Bernoulli likelihood----

# Bernoulli distribution (useful for when I want to add covariates and contact matrix)

#need to index data by time t AND individual N

# format data for STAN
#just use the first replicate
dat_Bern1 <- list(
  N = nrow(outInd), 
  t = tfinal,
  C = matrix(as.integer(outInd[,,1]=="C"), ncol = tfinal),
  S = matrix(as.integer(outInd[,,1]=="S"), ncol = tfinal),
  I = matrix(as.integer(outInd[,,1]=="I"), ncol = tfinal)
)
str(dat_Bern1) #make sure integers are classified as such


bern_NR = 
"
data { 
  //declare variables
  int<lower=0> N; //number of individuals
  int<lower=0> t; //number of time points
  int<lower = 0, upper =N> C[N,t]; //states
  int<lower = 0, upper =N> S[N,t]; 
  int<lower = 0, upper =N> I[N,t];
}
transformed data { 
  //Calculate transitions, based on change in states
  int<lower=0, upper = N> z_ci[N,t-1]; 
  int<lower=0, upper = N> z_si[N,t-1];
  for(T in 1:(t-1)){
    for(i in 1:N){
      z_ci[i,T] = C[i,T] - C[i,T+1]; //new infections of C
      z_si[i,T] = S[i,T] - S[i,T+1]; //new infections of S
    }
  }
}
parameters {
  real<lower=0, upper = 1> alpha;
  real<lower=0, upper = 1> beta;
}

model {
  real p_ci;
  real p_si;
  
  //priors
  alpha ~ beta(1,5); //might actually be a gamma dist
  beta ~ beta(1,5); //might actually be a gamma dist

 //generate number of infecteds from C and S at each time
 //C, S, and I are the number of individuals in that state at time t
 for(T in 1:(t-1)){
  for(i in 1:N){
    if(C[i,T]){ //transition for challenged
      p_ci = 1 - exp(-(alpha + beta*sum(I[,T])));
      z_ci[i,T] ~ bernoulli(p_ci); 
    }
    if(sum(I[,T]) > 0 && S[i,T]){ //transition for susceptible
      p_si = 1 - exp(-beta*sum(I[,T]));
      z_si[i,T] ~ bernoulli(p_si); 
      }
    }
 }
}
"
#run model
mod_bern_NR = stan_model(model_code = bern_NR)
fit_bern_NR = sampling(mod_bern_NR,dat_Bern1, iter = 1000, chains = 2, cores = 2)

#check it out
precis(fit_bern_NR) 
par(mfrow=c(1,2)) 
precis(fit_bern_NR) %>% plot; precis(fit_binom) %>% plot #binomial and bernoulli the same. dataset slightly different, accounting for differences. 
par(mfrow=c(1,1)) 
traceplot(fit_bern_NR)




#model replicated data with Bernoulli likelihood----

#model already takes long without any covariates!

# format data for STAN
dat_Bern2 <- list(
  N = nind,
  J = nreps,
  t = tfinal,
  ID = as.integer(1:nreps),
  C = abind(lapply(1:nreps, function(i) matrix(as.integer(outInd[,,i] == "C"), ncol = tfinal)), along = 3),
  S = abind(lapply(1:nreps, function(i) matrix(as.integer(outInd[,,i] == "S"), ncol = tfinal)), along = 3),
  I = abind(lapply(1:nreps, function(i) matrix(as.integer(outInd[,,i] == "I"), ncol = tfinal)), along = 3)
)
str(dat_Bern2) #make sure integers are classified as such


#try to fit where random effects are embedded within alpha and beta. individual estimates of alpha and beta look good, but the estimated mean (mu_a and mu_b) too high and doesn't make sense.
mod_bern1_2 <- stan_model(file = "IBM/scripts/model_fitting/sir_bernoulli1-2.stan")
start_time <- proc.time()
fit_bern1_2 <- sampling(mod_bern1_2, dat_Bern2, iter = 500, chains = 1, cores = 1)
total_time <- proc.time() - start_time

#check out model
precis(fit_bern1_2, depth = 2)
precis(fit_bern1_2, depth = 2, pars = c('alpha', 'beta')) %>% plot
traceplot(fit_bern1_2)


