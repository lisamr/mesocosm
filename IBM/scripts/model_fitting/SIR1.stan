data { 
  int<lower=1> N; //total number of individuals
  int<lower=1> J; //total number of trays
  int n[J]; //number of individuals in each tray
  int<lower=0> t; //number of time points
  int<lower=1> rep[N]; //tray rep
  int<lower=0> C[N]; //number of challenged
  int<lower=0> S[N]; //number of susceptible
  int<lower=0> I[N]; //number of infected
  real<lower=0> beta_max; //maximum value of beta (sometimes useful to help STAN converge)
}
transformed data { 
  //Calculate transitions, based on change in states
    int z_ci[J,t-1]; //new infections from C for each time and rep
    int z_si[J,t-1]; //new infections from S for each time and rep
    int c[J, t]; //matrix of challenged, turns vector into matrix
    int s[J, t]; //matrix of challenged, turns vector into matrix
    int i[J, t]; //matrix of challenged, turns vector into matrix

  for(x in 1:J){ //for every tray
    //declare variables of # of c, s, i in tray x
    int c[n[x]]; 
    int s[n[x]];
    int i[n[x]];
    
    //seperate individuals by tray
    c = segment(C, rep[x][1], )
    c = C[rep[x]];
    s = S[rep==tray];
    i = I[rep==tray];
    
    for(T in 1:(t-1)){
    z_ci[tray] = c[T] - c[T+1]; //new infections of C
    z_si[tray,] = s[T] - s[T+1]; //new infections of S
    }
  }
}
parameters {
  real<lower=0> alpha;
  real<lower=0, upper = beta_max> beta;
  vector[J] gamma;
  real<lower=0> sigma_g;
}

model {
  real p_ci;
  real p_si;
  
  //priors
  alpha ~ gamma(1,1);
  beta ~ uniform(0,beta_max);
  gamma ~ normal(0, sigma_g);
  sigma_g ~ exponential(1);

 //generate number of infecteds from C and S at each time
 //C, S, and I are the number of individuals in that state at time t
 for(i in 1:(t-1)){
    p_ci = 1 - exp(-(alpha + beta*I[i] + gamma[tray]));
    z_ci[i] ~ binomial(C[i], p_ci); 
  //only define z_si when there are infections. distribution is degenerate and STAN has trouble
    if(I[i] > 0){
      p_si = 1 - exp(-(beta*I[i] + gamma[tray]));
      z_si[i] ~ binomial(S[i], p_si); 
    }
 }
}