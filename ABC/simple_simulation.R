#write simple simulation model to be evaluated by ABC fitting
#single species, fixed transmission coefficients, NN? (gaussian kernel might be ok), 1 replicate

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

plot_maps(grid)


#visualize distance kernel
plot(seq(0,5,length.out = 100), K(seq(0,5,length.out = 100), 2), type = 'l')

#run simulation
f_sim1 <- function(beta, alpha, sigmak){
  
  #keep track of all individuals' states
  states_matrix <- matrix(NA, nrow=N, ncol = tfinal)
  states_matrix[,1] <- as.character(agents$state0)
  
  #vectors denoting specific states
  challenged <- as.numeric(agents$state0=='C')
  susceptible <- 1-challenged
  infected <- rep(0, N)
  
  #define distance kernel
  K1 <- K(d, sigmak)
  diag(K1) <- 0
  
  #Run the simulation
  for(t in 1:(tfinal-1)){ 
    
    #calculate probabilities of infections
    P_ci <- (1 - exp(-(alpha + (beta*K1)%*%infected)))*challenged
    P_si <- (1 - exp(- (beta*K1)%*%infected))*susceptible
    
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
D <- f_sim1(.75, 1, .5) #"True" data
Dstar <- f_sim1(.1, .5, .75) #simulated data
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
D <- f_sim1(alpha = 1, beta = .75, sigmak = .5) 

#set up parameters for ABC-SMC MNN
M <- 50 # Number of neighbours for covariance matrix calculations
P <- 100 #number of particles
n <- 1 #number of simulations per parameter set

#thresholds automatically determined by median distance
epsilonC <- c(1, rep(NA, 4))
epsilonS <- c(4, rep(NA, 4))
G <- length(epsilonC) #number of generations
z <- rep(1, G) #keep track of n.simulations run
y <- rep(0, G) #keep track of accepted simulations

#upper and lower bounds of parameters: alpha, beta, sigmak
#drastically narrow sigmak bounds
lower <- c(0,0,0)
upper <- c(2, 2, 2)


# Empty matrices to store results (2 mean differences + 3 model parameters)
res.old<-matrix(ncol=5,nrow=P)
res.new<-matrix(ncol=5,nrow=P)

# Empty vectors to store weights
w.old<-rep(NA, P)
w.new<-rep(NA, P)


#create point pattern object for dataset
npp <- ppp(x=agents$x, y = agents$y, c(min(agents$x), max(agents$x)),c(min(agents$y), max(agents$y)))

#calculate summary statistic of data
peak <- sort(unique(round(c(d), 8)))[3] #first peak 
CO <- sapply(1:tfinal, function(t) f_Oring(npp, D, t, 'C', R = peak)) 
SO <- sapply(1:tfinal, function(t) f_Oring(npp, D, t, 'S', R = peak)) 









#ABC fitting----
g = 1

for(g in 1:G){
  
  i = 1 #initiate particle counter 
  if(g != 1){
    #assign new epsilons
    epsilonC[g] = median(res.old[,4]) 
    epsilonS[g] = median(res.old[,5]) 
  } 
  
  while(i <= P){
    
    #sample parameters
    if(g == 1){ 
      alpha <- runif(1, lower[1], upper[1])
      beta <- runif(1, lower[2], upper[2])
      sigmak <- runif(1, lower[3], upper[3])
    }else{
      p<-sample(seq(1,P),1,prob=w.old) #index of theta*	
      sigma <- Sigma[[p]]
      #pull a draw of a, b from truncated MVN dist. mean = theta*
      par<- rtmvnorm(1, as.numeric(res.old[p,1:3]), sigma, lower, upper) 
      alpha <- par[1]
      beta <- par[2]
      sigmak <- par[3]
    }
    theta_star <- c(alpha, beta, sigmak)
    
    #make sure parameters>0
    if(prior.non.zero(theta_star)){
      m = 0 #restart number of accepted simulations
      
      distanceC <- rep(NA, n)
      distanceS <- rep(NA, n)
      for(j in 1:n){
        D_star <- f_sim1(alpha = alpha, beta = beta, sigmak = sigmak) #run model
        z[g]=z[g]+1 #track n.simulations/generation
        
        #calc difference, accept if all thresholds met
        COstar <- sapply(1:tfinal, function(t) f_Oring(npp, D_star, t, 'C', R = peak)) 
        SOstar <- sapply(1:tfinal, function(t) f_Oring(npp, D_star, t, 'S', R = peak)) 
        distanceC[j] <- SSE(CO, COstar) 
        distanceS[j] <- SSE(SO, SOstar)
        
        if((distanceC[j] <= epsilonC[g]) & (distanceS[j] <= epsilonS[g])) m = m+1
        
      }
      y[g] = y[g]+m #track n.accepted simulations
      
      if(m>0){
        #store parameters
        res.new[i,] <- c(theta_star, mean(distanceC), mean(distanceS))
          
        
        #caculate weights
        w1<-prod(sapply(1:3, function(b) dunif(res.new[i,b], min=lower[b], max=upper[b]))) 
        if(g==1){w2 <- 1
        }else{
          w2 <- sum(sapply(1:P, function(a) w.old[a]* dtmvnorm(res.new[i,1:3], mean=res.old[a,1:3], sigma=sigma, lower=lower, upper=upper)))
          }
        w.new[i] <- m/n*w1/w2
          
        #update counter and print progress
        print(paste0('Generation: ', g, ', particle: ', i, ', acceptance rt: ', y[g]/z[g]))
        i = i+1
        }
    }
  }
  
  #after obtaining P particles, finish the generation 
  Sigma <- vector('list', length = P)
  for(p in 1:P) Sigma[[p]]<- getSigmaNeighbours(M, res.new[p,1:3], res.new[,1:3]) 
  res.old <- res.new
  w.old <- w.new/sum(w.new)
  
  write.csv(res.new, file = paste0("ABC/results/kernel_ABC_peak2_gen_",g,".csv"), row.names=FALSE)
  
}

#total simulations, accepted simulations per generation
print(z);print(y) 
beepr::beep(5)



#check out posteriors-----
# alpha = 1, beta = .75, sigmak = .5
#3:18pm

plotdens <- function(gen){
  par(mfrow=c(1,3))
  dens(gen$V1, main = 'alpha')
  abline(v=1, lwd=2, col='red')
  dens(gen$V2, main = 'beta')
  abline(v=.75, lwd=2, col='red')
  dens(gen$V3, main = 'sigma')
  abline(v=.5, lwd=2, col='red')
  par(mfrow=c(1,1))
}


gen1 <- read_csv('ABC/results/kernel_ABC_peak2_gen_1.csv')
gen2 <- read_csv('ABC/results/kernel_ABC_peak2_gen_2.csv')
gen3 <- read_csv('ABC/results/kernel_ABC_peak2_gen_3.csv')
gen4 <- read_csv('ABC/results/kernel_ABC_peak2_gen_4.csv')
gen5 <- read_csv('ABC/results/kernel_ABC_peak2_gen_5.csv')

plotdens(gen1)
plotdens(gen2)
plotdens(gen3)
plotdens(gen4)
plotdens(gen5)

gen5$V3 %>% hist



#simulate infections with posterior
sim_post <- function(i, gen){
  tmp <- f_sim1(alpha = gen$V1[i], beta=gen$V2[i], sigmak = gen$V3[i])
  colSums(tmp=='I') %>% points(type='l', pch=16, col=scales::alpha('blue', .2))
}



#20 different parameter sets. need to narrow epsilon!
Nsims <- 50
colSums(D=='I') %>% plot(type='o', pch=16)
sapply(1:Nsims, function(i) sim_post(i, gen5))

#check out the demographic stochasticity
colSums(D=='I') %>% plot(type='o', pch=16)
replicate(20, sim_post(1, gen5)) #20 different parameter sets


sim_Oring <- function(gen, i){
  D_star <- f_sim1(alpha = gen$V1[i], beta=gen$V2[i], sigmak = gen$V3[i])
  COstar <- sapply(1:tfinal, function(t) f_Oring(npp, D_star, t, 'C', R = peak)) 
  SOstar <- sapply(1:tfinal, function(t) f_Oring(npp, D_star, t, 'S', R = peak)) 
  return(data.frame(COstar, SOstar))
}

test <- lapply(1:Nsims, function(i) sim_Oring(gen5, i))

par(mfrow=c(1,3))
plot(CO, type='o', pch=16, ylim=c(0,1))
sapply(1:Nsims, function(i) lines(test[[i]]$COstar, col=scales::alpha('blue', .2) ))

plot(SO, type='o', pch=16, ylim=c(0,2))
sapply(1:Nsims, function(i) lines(test[[i]]$SOstar, col=scales::alpha('blue', .2) ))

colSums(D=='I') %>% plot(type='o', pch=16)
sapply(1:Nsims, function(i) sim_post(i, gen5))
par(mfrow=c(1,1))

#check out how the summary statistic changes with same parameter set. quite a lot :(
test2 <- replicate(10, sim_Oring(gen5, 1)$SOstar)
colrank <- rank(sapply(1:10, function(i) SSE(SO, test2[,i])))
palz <- viridis::viridis(10)
plot(SO, type='o', pch=16, ylim=c(0,2))
sapply(1:10, function(i) lines(test2[,i], col=scales::alpha(palz[colrank][i], .5) ))

#do it with known parameters?
test3C <- replicate(40, sim_Oring(data.frame(V1=1, V2=.75, V3=.5), 1)$COstar)
plot(CO, type='o', pch=16, ylim=c(0,2))
sapply(1:40, function(i) lines(test3C[,i], col=scales::alpha('blue', .5) ))
sapply(1:40, function(i) SSE(CO, test3C[,i])) %>% dens

test3S <- replicate(40, sim_Oring(data.frame(V1=1, V2=.75, V3=.5), 1)$SOstar)
plot(SO, type='o', pch=16, ylim=c(0,2))
sapply(1:40, function(i) lines(test3S[,i], col=scales::alpha('blue', .5) ))
sapply(1:40, function(i) SSE(SO, test3S[,i]))




#figure out appropriate initial epsilons for CO and SO. 
sapply(1:Nsims, function(i) SSE(CO, test[[i]]$COstar)) %>% dens #1?
sapply(1:Nsims, function(i) SSE(SO, test[[i]]$SOstar)) %>% dens #4?

dens(gen5$V4)

SSE(CO, test[[3]]$COstar)

