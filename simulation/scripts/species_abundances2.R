#what does it mean that communities conform to a lognormal distribution? what is log-normal? I think the histogram of the species pool is lognormal, and the ranked abundance curve is based off of that. I've been able to replicate what the canned function does. 
#draw random points from a lognormal distribution (the species pool) to form your local one. average across 1000 simulations.
library(mobsim)

par(mfrow=c(1,3))
par(mfrow=c(1,1))
#what you have been doing:
x <- dlnorm(1:6, meanlog = 1,sdlog = 1, log = F)
plot(x, type='o', col='blue')

#use a canned function to generate lognormal distribution. then figure out what the black box does. s_pool=#species, n_sim=#individuals
abund1 <- sim_sad(s_pool = 6, n_sim = 1000, sad_type = "lnorm",
                  sad_coef = list("meanlog" = 1, "sdlog" = .5))
abund1
sum(abund1)
abund2 <- as.vector(abund1)
plot(abund1) #distribution of the global species pool
plot(abund2, type='o') #ranked abundance curve.
lines(abund2, col='blue')
plot(abund2/sum(abund2), type='o')

N <- 6
sims <- 1000
m <- matrix(NA, nrow = sims, ncol = N)
for(i in 1:sims){
  abund1 <- sim_sad(s_pool = N, n_sim = 1000, sad_type = "lnorm",
                      sad_coef = list("meanlog" = 1, "sdlog" = .5))
  abund2 <- as.vector(abund1)
  m[i,] <- abund2/sum(abund2)
}
y <- apply(m, 2, mean) 
#plot(y, type='o')
lines(y, type='l')
#Try replicating it with basic R functions. it's identical. The right thing to do
N <- 6
sims <- 1000
m2 <- matrix(NA, nrow = sims, ncol = N)
for(i in 1:sims){
  y2 <- rlnorm(1:N, meanlog = 1, sdlog = .3) %>% sort(decreasing = T)
  m2[i,] <- y2
}
y2 <- apply(m2, 2, mean)/sum(apply(m2, 2, mean))
lines(y2, type='o', col="green")


#ok. so the numbers aren't off all that much, so I don't think it will affect your results. Still, you do need to fix it.
dat <- read.csv("simulation/outputs/hostabund.csv")
z <- dat$R4 %>% table*3
sum(z)

#compare numbers across methods
round(x*300)#current
round(y*300)#canned function
round(y2*300)#replicated with base R functions
plot(round(x*300), ylim=c(0, 100), type='o')
lines(round(y*300), type='o', col='blue')
lines(round(y2*300), type='o', col='green')
lines(z, type='o', col='red')

sum(x)#doesn't sum to 1.
sum(y)
sum(y2)


