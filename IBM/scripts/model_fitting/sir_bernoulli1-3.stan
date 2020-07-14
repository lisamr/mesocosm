data{ 
  //declare variables
  int<lower=0> N; //number of individuals
  int<lower=0> t; //number of time points
  int<lower=0> J; //number of replicated trays
  int ID[J]; //tray replicate ID
  int<lower = 0, upper =1> C[N,t,J]; //states in array
  int<lower = 0, upper =1> S[N,t,J]; 
  int<lower = 0, upper =1> I[N,t,J];
}

transformed data{ 
  //track when each individual makes a state change (new infection)
  int<lower=0, upper = 1> z_ci[N,t-1,J]; 
  int<lower=0, upper = 1> z_si[N,t-1,J];
  for(j in 1:J){
    for(T in 1:(t-1)){
      for(i in 1:N){
        z_ci[i,T,j] = C[i,T,j] - C[i,T+1,j]; 
        z_si[i,T,j] = S[i,T,j] - S[i,T+1,j]; 
      }
    }
  }
}

parameters{
  real<lower=0, upper=1> alpha; 
  real<lower=0, upper=1> beta;
  real<lower=0> za[J];
  real<lower=0> zb[J];
  real<lower=0> scale_a; 
  real<lower=0> scale_b; 
}

model{
  real p_ci;
  real p_si;
  
  //priors
  alpha ~ beta(1, 5); //I think needs to be 0-1
  beta ~ beta(1, 5);  
  za ~ gamma(1, scale_a); //force varying int to center around 1
  scale_a ~ exponential(2);
  zb ~ gamma(1, scale_b); 
  scale_b ~ exponential(2);
  //generate number of infecteds from C and S at each time
 //C, S, and I are the number of individuals in that state at time t
 for(j in 1:J){
   for(T in 1:(t-1)){
     for(i in 1:N){
       if(C[i,T,j]){ //transition for challenged
       p_ci = 1 - exp(-(alpha*za[ID[j]] + beta*zb[ID[j]]*sum(I[,T,j]) ));
       z_ci[i,T,j] ~ bernoulli(p_ci); 
       }
       
       if(sum(I[,T,j]) > 0 && S[i,T,j]){ //for susceptible
       p_si = 1 - exp(-(beta*zb[ID[j]]*sum(I[,T,j]) ));
       z_si[i,T,j] ~ bernoulli(p_si); 
      }
     }
   }
 }
} 
