#write simple simulation model to be evaluated by ABC fitting
#single species, fixed transmission coefficients, NN? (gaussian kernel might be ok), 1 replicate

rm(list=ls())
source('IBM/scripts/IBM_functions.R')


#functions
K <- function(D, sigma) exp(-sigma*D^2) #distance kernel
plot(seq(0,5,length.out = 100), K(seq(0,5,length.out = 100), 1), type = 'l')
headmat <- function(mat, n=10) mat[1:n,1:n]


#starting parameters
r <- 2 #interplanting distance in cm
pinoc <- .1 #percent inoculated
s <- 1 #number of spp
spp <- "sp1" #species name
tfinal <- 20 #time steps

#simulate communities in a 9.5in square tray in hexagonal grid
grid <- sample_community(1, pinoc, r) 
agents <- as.data.frame(grid) %>% mutate(ID = row_number())
N <- nrow(agents)

#distance matrix
d <- as.matrix(dist(cbind(agents$x, agents$y), method = 'euclidean', diag = TRUE, upper = TRUE))



#run simulation
f_sim1 <- function(beta, alpha, sigma){
  
  #keep track of all individuals' states
  states_matrix <- matrix(NA, nrow=N, ncol = tfinal)
  states_matrix[,1] <- as.character(agents$state0)
  
  #vectors denoting specific states
  challenged <- as.numeric(agents$state0=='C')
  susceptible <- 1-challenged
  infected <- rep(0, N)
  
  #define distance kernel
  K1 <- K(d, sigma)
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


#generate fake data----
beta <- .75 
alpha <- 1
sigma <- .5
D <- f_sim1(beta, alpha, sigma) #"True" data
Dstar <- f_sim1(1, .5, .6) #simulated data
plot_spread_map(grid, Dstar, animate = F)

#compare "data" with simulation
tmp <- data.frame(I = colSums(D=="I"),
           Istar = colSums(Dstar=="I"),
           t = 1:tfinal)
plot(tmp$t, tmp$I, type = 'o', col = 'blue')
lines(tmp$t, tmp$Istar, type = 'o', col = 'red')
legend('topleft', legend = c('D', 'D*'), fill = c('blue', 'red'))



#candidate summary statistics----

#difference functions
SSE <- function(S, Sstar){
  sqrt(sum(S - Sstar)^2)
}
chisq <- function(S, Sstar){ #pg 11 McKinley et al. 2009
  x <- ((S - Sstar)^2)/S
  x[is.na(x)] <- 0
  sum(x)
}


# 1. number of infecteds
nI <- function(x){
  sum(x=='I')
}

I <- apply(D, 2, nI)
Istar <- apply(Dstar, 2, nI)

SSE(I, Istar)
chisq(I, Istar)

sqrt((apply(D, 2, nI) - apply(Dstar, 2, nI))^2) %>% plot


# 2. nearest infected neighbor from Minter 2019
#not going to be good if most infections are from nearest neighbor!

uniq.dist <- sort(unique(round(c(d), 8)))
uniq.dist

nobs <- function(D){
  n.obs <- matrix(NA, nrow = length(uniq.dist), ncol = ncol(D))
  for(i in 1:ncol(D)){
    inf <- which(D[,i]=="I")
    if(length(inf)>0){
      NinfN <- apply(d[inf, inf], 1, function(x) round(min(x[x!=0]), 8))
    }else NinfN <- 0
    n.obs[,i] <- sapply(1:length(uniq.dist), function(a) sum(NinfN==uniq.dist[a]))
  }
  return(n.obs)
}

sqrt(colSums((nobs(D) - nobs(Dstar))^2)) %>% plot

NN <- nobs(D)
NNstar <- nobs(Dstar)







# 3. o-ring statistic
require(spatstat)

#create point pattern object for dataset
npp <- ppp(x=agents$x, y = agents$y, c(min(agents$x), max(agents$x)),
           c(min(agents$y), max(agents$y)))

#calculate O-ring for each time at specified distance
f_Oring <- function(ppp, Data, t, R, Spar=.7){
  marks(npp) <- as.factor(Data[,t]) #update time
  K12 <- Kcross(npp, 'I', 'I') #calculate o-ring stat 
  g12 <- pcf(K12, method="a", spar=Spar)
  lambda2 <- summary(npp)$marks['I',"intensity"]
  Oring <- eval.fv(lambda2*g12)
  res <- Oring$pcf[which.min(abs(Oring$r - 2))] 
  return(res)
}

#arbitrarily choose r=2
O <- sapply(2:20, function(t) f_Oring(npp, D, t, 2))
Ostar <- sapply(2:20, function(t) f_Oring(npp, Dstar, t, 2))

#doesn't really seem to differ from sum(I_t)
plot(O, type = 'o', col='blue', main = 'Oring')
points(Ostar, type = 'o', col='red')




# 4. ripley's k

#cumulative number of events across a range of distances. 
fK <- function(Marks, R, C=c("border", "isotropic", "Ripley", "translate")){
  marks(npp) <- as.factor(Marks)
  tmp <- Kest(subset(npp, marks == 'I'), r = c(0,R), correction=C) 
  tmp[[3]][2]
}

#arbitrarily choose r=2.5
rK <- sapply(2:20, function(t) fK(D[,t], R = 2.5, C = 'isotropic'))
rKstar <- sapply(2:20, function(t) fK(Dstar[,t], R = 2.5, C = 'isotropic'))
plot(rK, type = 'o', col='blue', main = 'Oring', ylim=c(min(c(rKstar, rK)), max(c(rKstar, rK))))
points(rKstar, type = 'o', col='red')



#compare all 4 summary statistics visually

par(mfrow=c(2,2))
plot(I, type = 'o', col='blue', main = '# infected' )
points(Istar, type = 'o', col='red')

plot(NN, type = 'o', col='blue', main = 'nearest infected neighbor', sub = 'doesnt work, matrix')
points(NNstar, type = 'o', col='red')

plot(O, type = 'o', col='blue', main = 'Oring', ylim=c(min(c(Ostar, O)), max(c(Ostar, O))))
points(Ostar, type = 'o', col='red')

plot(rK, type = 'o', col='blue', main = 'K function', ylim=c(min(c(rKstar, rK)), max(c(rKstar, rK))))
points(rKstar, type = 'o', col='red')



#plot difference function over time.
dI <- data.frame(d = sqrt((I - Istar)^2), t=1:length(I), metric = 'I')
dNN <-data.frame(d = sqrt(colSums((NN - NNstar)^2)), t=1:ncol(NN), metric = 'NN') 
dO <-data.frame(d = sqrt((O - Ostar)^2), t=2:(length(O)+1), metric = 'O-ring') 
drK <-data.frame(d = sqrt((rK - rKstar)^2), t=2:(length(rK)+1), metric = 'ripleys K') 
df <- bind_rows(dI, dNN, dO, drK)
head(df)
ggplot(df, aes(t, d, group = metric)) +
  geom_point() +
  geom_line(aes(color = metric)) +
  facet_wrap(~metric, scales = 'free_y')

