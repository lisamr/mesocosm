#write simple simulation model to be evaluated by ABC fitting
#single species, fixed transmission coefficients, and nearest neigbor, 1 replicate

rm(list=ls())
source('ABC/ABC_functions.R')
require(spatstat)
require(sp)
library(tidyverse)



#setup community grid----
r <- 1.7 #interplanting distance in cm
pinoc <- .1 #percent inoculated
s <- 1 #number of spp
spp <- c("sp1") #species name
tfinal <- 20 #time steps

#simulate communities in a 9.5in square tray in hexagonal grid
grid <- sample_community(spp, pinoc, r) 
agents <- as.data.frame(grid) %>% mutate(ID = row_number())
N <- nrow(agents)

#distance matrix
d <- as.matrix(dist(cbind(agents$x, agents$y), method = 'euclidean', diag = TRUE, upper = TRUE))

#Nearest neighbor matrix
dnn <- round(d, 8)
dnn[dnn==0] <- NA
NN <- function(v){
  nn <- which(v == min(v,na.rm = T))
  v[nn] <- 1
  v[-nn] <- 0
  return(v)
}
dnn <- sapply(1:nrow(dnn), function(x) NN(dnn[x,]))

ggplot(agents, aes(x,y)) +
  geom_point(aes(color=state0)) +
  coord_equal()


#run simulation
f_sim1 <- function(alpha, beta){
  
  #keep track of all individuals' states
  states_matrix <- matrix(NA, nrow=N, ncol = tfinal)
  states_matrix[,1] <- as.character(agents$state0)
  
  #vectors denoting specific states
  challenged <- as.numeric(agents$state0=='C')
  susceptible <- 1-challenged
  infected <- rep(0, N)
  
  
  #Run the simulation
  for(t in 1:(tfinal-1)){ 
    
    #calculate probabilities of infections
    P_ci <- (1 - exp(-(alpha + (beta*dnn)%*%infected)))*challenged
    P_si <- (1 - exp(- (beta*dnn)%*%infected))*susceptible
    
    #generate new infections
    CI <- rbinom(N, 1, P_ci)
    SI <- rbinom(N, 1, P_si)
    new_infections <- CI + SI
    
    #update states
    challenged <- challenged - CI
    susceptible <- susceptible - SI
    infected <- infected + new_infections
    states_matrix[,t+1][challenged==1] <- 'C'
    states_matrix[,t+1][susceptible==1] <- 'S'
    states_matrix[,t+1][infected==1] <- 'I'
  }
  
  return(states_matrix)
}


#test generate fake data----
D <- f_sim1(1, .75) #"True" data
Dstar <- f_sim1(.5, .1) #simulated data
plot_spread_map(grid, D, animate = F)


#compare "data" with simulation
tmp <- data.frame(S = colSums(D=="S"),
                  Sstar = colSums(Dstar=="S"),
                  C = colSums(D=="C"),
                  Cstar = colSums(Dstar=="C"),
                  I = colSums(D=="I"),
                  Istar = colSums(Dstar=="I"),
                  t = 1:tfinal)

par(mfrow = c(1,3))
plot(tmp$t, tmp$S, type = 'o', col = 'blue', main = 'S->I')
lines(tmp$t, tmp$Sstar, type = 'o', col = 'red')
plot(tmp$t, tmp$C, type = 'o', col = 'blue', pch=16, main = 'C->I')
lines(tmp$t, tmp$Cstar, type = 'o', col = 'red', pch=16)
plot(tmp$t, tmp$I, type = 'o', col = 'blue', pch=16, main = 'I')
lines(tmp$t, tmp$Istar, type = 'o', col = 'red', pch=16)





#ABC parameters----

#simulate 'True' data
D <- f_sim1(alpha = 1, beta = .75) 

#set up parameters for ABC-SMC MNN
M <- 50 # Number of neighbours for covariance matrix calculations
P <- 1000 #number of particles
n <- 1 #number of simulations per parameter set

#thresholds automatically determined by median distance
epsilonC <- c(8, rep(NA, 6))
epsilonS <- c(200, rep(NA, 6))
G <- length(epsilonC) #number of generations
z <- rep(1, G) #keep track of n.simulations run
y <- rep(0, G) #keep track of accepted simulations

#upper and lower bounds of parameters: alpha, beta, sigmak
#drastically narrow sigmak bounds
lower <- c(0,0)
upper <- c(2, 5)


# Empty matrices to store results (2 mean differences + 3 model parameters)
res.old<-matrix(ncol=4,nrow=P)
res.new<-matrix(ncol=4,nrow=P)

# Empty vectors to store weights
w.old<-rep(NA, P)
w.new<-rep(NA, P)


#calculate summary statistic of data (C and S counts)
nC <- colSums(D=='C')
nS <- colSums(D=='S')






#ABC fitting----
g = 1

for(g in 1:G){
  
  i = 1 #initiate particle counter 
  if(g != 1){
    #assign new epsilons
    epsilonC[g] = median(res.old[,3]) 
    epsilonS[g] = median(res.old[,4]) 
  } 
  
  while(i <= P){
    
    #sample parameters
    if(g == 1){ 
      alpha <- runif(1, lower[1], upper[1])
      beta <- runif(1, lower[2], upper[2])
    }else{
      p<-sample(seq(1,P),1,prob=w.old) #index of theta*	
      sigma <- Sigma[[p]]
      #pull a draw of a, b from truncated MVN dist. mean = theta*
      par<- rtmvnorm(1, as.numeric(res.old[p,1:2]), sigma, lower, upper) 
      alpha <- par[1]
      beta <- par[2]
    }
    theta_star <- c(alpha, beta)
    
    #make sure parameters>0
    if(prior.non.zero(theta_star)){
      m = 0 #restart number of accepted simulations
      
      distanceC <- rep(NA, n)
      distanceS <- rep(NA, n)
      for(j in 1:n){
        D_star <- f_sim1(alpha, beta) #run model
        z[g]=z[g]+1 #track n.simulations/generation
        
        #calc difference, accept if all thresholds met
        nCstar <- colSums(D_star=='C')
        nSstar <- colSums(D_star=='S')

        distanceC[j] <- SSE(nC, nCstar)
        distanceS[j] <- SSE(nS, nSstar)
        
        if((distanceC[j] <= epsilonC[g]) & (distanceS[j] <= epsilonS[g])) m = m+1
        
      }
      y[g] = y[g]+m #track n.accepted simulations
      
      if(m>0){
        #store parameters
        res.new[i,] <- c(theta_star, mean(distanceC), mean(distanceS))
        
        
        #caculate weights
        w1<-prod(sapply(1:2, function(b) dunif(res.new[i,b], min=lower[b], max=upper[b]))) 
        if(g==1){w2 <- 1
        }else{
          w2 <- sum(sapply(1:P, function(a) w.old[a]* dtmvnorm(res.new[i,1:2], mean=res.old[a,1:2], sigma=sigma, lower=lower, upper=upper)))
        }
        w.new[i] <- m/n*w1/w2
        
        #update counter and print progress
        print(paste0('Generation: ', g, ', particle: ', i, ', acceptance rt: ', round(y[g]/z[g], 2) ))
        i = i+1
      }
    }
  }
  
  #after obtaining P particles, finish the generation 
  Sigma <- vector('list', length = P)
  for(p in 1:P) Sigma[[p]]<- getSigmaNeighbours(M, res.new[p,1:2], res.new[,1:2]) 
  res.old <- res.new
  w.old <- w.new/sum(w.new)
  
  write.csv(res.new, file = paste0("ABC/results/ABC_NN_n1_gen_",g,".csv"), row.names=FALSE)
  
}

#total simulations, accepted simulations per generation
print(z);print(y) 
beepr::beep(5)
sum(z)



#check out posteriors-----
# alpha = 1, beta = .75, sigmak = .5
#3:18pm

plotdens <- function(gen){
  par(mfrow=c(2,2))
  dens(gen$V1, main = paste('alpha'))
  abline(v=1, lwd=2, col='red')
  abline(v=median(gen$V1), lwd=2, lty=2, col='blue')
  dens(gen$V2, main = 'beta')
  abline(v=.75, lwd=2, col='red')
  abline(v=median(gen$V2), lwd=2, lty=2, col='blue')
  dens(gen$V3, main = 'diff C')
  dens(gen$V4, main = 'diff S')
  par(mfrow=c(1,1))
}


gen1 <- read_csv('ABC/results/ABC_NN_n1_gen_1.csv')
gen2 <- read_csv('ABC/results/ABC_NN_n1_gen_2.csv')
gen3 <- read_csv('ABC/results/ABC_NN_n1_gen_3.csv')
gen4 <- read_csv('ABC/results/ABC_NN_n1_gen_4.csv')
gen5 <- read_csv('ABC/results/ABC_NN_n1_gen_5.csv')
gen6 <- read_csv('ABC/results/ABC_NN_n1_gen_6.csv')
gen7 <- read_csv('ABC/results/ABC_NN_n1_gen_7.csv')
gen8 <- read_csv('ABC/results/ABC_NN_n1_gen_8.csv')

plotdens(gen1)
plotdens(gen2)
plotdens(gen3)
plotdens(gen4)
plotdens(gen5)
plotdens(gen6)
plotdens(gen7)
plotdens(gen8)


results <- bind_rows(gen1, gen2, gen3, gen4, gen5, gen6, gen7) %>% 
  mutate(generation = rep(1:7, each = P),
         particle = rep(1:P, 7))

p1 <- ggplot(results, aes(V1, group=generation, fill=as.factor(generation))) +
  geom_density(alpha=.5) +
  scale_fill_brewer(palette = 'YlOrRd') +
  geom_vline(xintercept = 1, lty=2, color='blue') +
  labs(x='alpha', fill='gen')

p2 <- ggplot(results, aes(V2, group=generation, fill=as.factor(generation))) +
  geom_density(alpha=.5) +
  scale_fill_brewer(palette = 'YlOrRd') +
  geom_vline(xintercept = .75, lty=2, color='blue') +
  labs(x='beta', fill='gen')

cowplot::plot_grid(p1, p2)
mean(gen7$V2)
mean(gen7$V1)


#simulate infections with posterior
sim_post <- function(i, gen, color='blue'){
  tmp <- f_sim1(alpha = gen$V1[i], beta=gen$V2[i])
  colSums(tmp=='I') %>% points(type='l', pch=16, col=scales::alpha(color, .2))
}



#20 different parameter sets. need to narrow epsilon!
Nsims <- 50
colSums(D=='I') %>% plot(type='o', pch=16)
sapply(1:Nsims, function(i) sim_post(i, gen3))
sapply(1:Nsims, function(i) sim_post(i, gen7, 'red'))


#check out the demographic stochasticity
colSums(D=='I') %>% plot(type='o', pch=16)
replicate(20, sim_post(1, gen5)) #20 different parameter sets


#simulate loss function
D_star7 <- lapply(1:400, function(i) {
  dstar = f_sim1(alpha = gen7$V1[i], beta=gen7$V2[i])
  return(dstar)
})
D_sametheta <- replicate(400, f_sim1(1, .75))

simn <- function(d){
  nc <- colSums(d=='C') 
  ns <- colSums(d=='S') 
  return(data.frame(nc, ns))
}
D_star7n <- lapply(1:Nsims, function(i) simn(D_star7[[i]])) 
D_samethetan <- lapply(1:Nsims, function(i) simn(D_sametheta[,,i])) 





par(mfrow=c(1,2))
plot(nC, type = 'l', lwd = 2, xlab='time', ylab='# Challenged', main='fitted simulations')
for(i in 1:Nsims) lines(D_star7n[[i]]$nc, col=scales::alpha('red', .2))
plot(nC, type = 'l', lwd = 2, xlab='time', ylab='# Challenged', main='same theta')
for(i in 1:Nsims) lines(D_samethetan[[i]]$nc, col=scales::alpha('blue', .2))


plot(nS, type = 'l', lwd = 2, xlab='time', ylab='# Susceptibles', main='fitted simulations')
for(i in 1:Nsims) lines(D_star7n[[i]]$ns, col=scales::alpha('red', .2))
plot(nS, type = 'l', lwd = 2, xlab='time', ylab='# Susceptibles', main='same theta')
for(i in 1:Nsims) lines(D_samethetan[[i]]$ns, col=scales::alpha('blue', .2))



#check out other loss functions?

#1. nearest infected neighbor?

nin <- sapply(1:400, function(i) sqrt(colSums(nobs(D, uniq.dist(d)) - nobs(D_star7[[i]], uniq.dist(d)))^2))

#compare to data with known theta
D_sametheta <- replicate(400, f_sim1(1, .75))
ninsametheta <- sapply(1:400, function(i) sqrt(colSums(nobs(D, uniq.dist(d)) - nobs(D_sametheta[,,i], uniq.dist(d)))^2))

#distribution of errors
colSums(nin) %>% hist(main='fitted data') 
colSums(ninsametheta) %>% hist

colSums(nin) %>% median
colSums(ninsametheta) %>% median
        