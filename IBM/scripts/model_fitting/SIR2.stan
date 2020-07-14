
functions{

    //convert a distance matrix into a matrix containing P(contact) values
  matrix f_kernel(int N, matrix dist_mat, real Sigma){
    matrix[N,N] contact_rt;
    matrix[N,N] dist_mat2;
    
    dist_mat2 = dist_mat .* dist_mat; //elementwise product of matrix x and y
    contact_rt = exp(-Sigma * dist_mat2);
    for(k in 1:N){ //reassign diagonals as 0
      contact_rt[k,k] = 0;
    }
    return(contact_rt);
  }
}

data { 
  //declare variables
  int<lower=0> N; //number of individuals
  int<lower=0> t; //number of time points
  int<lower = 0, upper =N> C[N,t]; //states
  int<lower = 0, upper =N> S[N,t]; 
  int<lower = 0, upper =N> I[N,t];
  matrix[N,N] dist_mat; //distance matrix. used to calculate probability of contact.
  real<lower = 0, upper=1> sigma_upper;
}

transformed data { 
  //Calculate transitions, based on change in states
  int<lower=0, upper = N> z_ci[N,t-1]; 
  int<lower=0, upper = N> z_si[N,t-1];
  
  //calculate infected matrix for each individual. I is treated as an integer array, not matrix
  matrix[N,t] I_mat;
  I_mat = to_matrix(I); //will be multipled by contact matrix in model
  
  //calculate newly infected C and S individuals
  for(T in 1:(t-1)){
    for(i in 1:N){
      z_ci[i,T] = C[i,T] - C[i,T+1]; //new infections of C
      z_si[i,T] = S[i,T] - S[i,T+1]; //new infections of S
    }
  }
}

parameters {
  real<lower=0, upper=1> amp_a; //.6
  //real<lower=0> rate_a; //had problems estimating. going to fix it at .5
  real<lower=0, upper=1> amp_b; //.1
  real<lower=0> rate_b; //1
  //real<lower=0> tq; //10
  real<lower=0, upper=sigma_upper> sigma; //.005

}

model {
  //declare variables
  real alpha;
  real beta;
  real p_ci;
  real p_si;
  matrix[N,N] contact_mat;
  matrix[N,t] I_mat_contact;
  
  //priors. simulated values in comments
  amp_a ~ beta(1,2); 
  //rate_a ~ gamma(1,.5); 
  amp_b ~ beta(1,2); 
  rate_b ~ gamma(1,.5); 
  //tq ~ gamma(3,.3); 
  sigma ~ beta(1,100); //most values under .1

  //calculate sum of infecteds weighted by contact matrix 
  contact_mat = f_kernel(N, dist_mat, sigma);
  //contact_mat = f_kernel(N, dist_mat, .005);
  I_mat_contact = contact_mat * I_mat;
  
 for(T in 1:(t-1)){
  //calculate alpha and beta at each time point
  alpha = amp_a*exp(-.5*T);
  beta = amp_b*exp(-rate_b*(log(T/10.0))^2);
  //alpha = .6*exp(-.5*T);
  //beta = .1*exp(-1*(log(T/10.0))^2);
 
  for(i in 1:N){
    if(C[i,T]){ //transition for challenged
      p_ci = 1 - exp(-(alpha + beta*I_mat_contact[i,T]));
      z_ci[i,T] ~ bernoulli(p_ci); 
    }
    if(I_mat_contact[i,T] > 0 && S[i,T]){ //transition for susceptible
      p_si = 1 - exp(-beta*I_mat_contact[i,T]);
      z_si[i,T] ~ bernoulli(p_si); 
      }
    }
 }
}

