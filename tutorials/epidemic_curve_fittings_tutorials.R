# https://kingaa.github.io/clim-dis/parest/parest.html --------------------

#fitting sir measles data. from Aaron King's workshop in 2017. Apparently assumed to have frequency-dependent transmission. 

#check out data
niamey <- read.csv("http://kingaa.github.io/clim-dis/parest/niamey.csv")
ggplot(niamey,mapping=aes(x=biweek,y=measles,color=community))+
  geom_line()+
  geom_point() 

#check out one community on the log scale
ggplot(niamey %>% filter(community == 'A'),mapping=aes(x=biweek,y=measles))+
  geom_line()+
  geom_point() + 
  scale_y_log10()

#the slope for the beginning parts of the outbreak (t<8) are proportional to R0-1 and the recovery rate. Get a cheap approximation of R0 with a linear regression.
fit1 <- lm(log(measles)~biweek,data=subset(niamey,biweek<=8&community=="A"))
summary(fit1) 

slope <- coef(fit1)[2]
slope 

#(R0-1)(gamma+mu). mu (death) can be ignored cuz mu<<gamma (recovery). Previous knowledge that infectious period is 2 weeks (1 unit of biweekly time.)
#R0 = slope/gamma + 1 = .43/1 + 1 = 1.4


R0finalsize <- function(f){
  R0 = -(log(1-f)) / f
  return(R0)
}
R0finalsize(.99)


#http://web.stanford.edu/class/earthsys214/notes/fit.html-------------

#uses ode solver and looks at fitting data to logistic curves

#first fit to a logistic curve with SSE----
#dN/dt = rN(1-N/K) #continuous time
#N(t) = (N0*exp(rt)) / (1+N0*(exp(rt) - 1)/K) #analytical solution

#custom function with SSE as loss function
fit.logistic <- function(par,y){
  r <- par[1]; k <- par[2]
  t <- y[,2]
  n <- y[,1]
  n0 <- y[1]
  tmp <- n[1] *exp(r*t)/(1 + n[1] * (exp(r*t)-1)/k)
  sumsq <- sum((n - tmp)^2)
}

#data (US population census)
usa <- c(3.929214,   5.308483,   7.239881,   9.638453,  12.866020,  
17.866020, 23.191876,  31.443321,  38.558371,  50.189209,
62.979766,  76.212168, 92.228496, 106.021537, 123.202624,
132.164569, 151.325798, 179.323175, 203.302031, 226.542199,
248.709873, 281.421906, 308.745538)
year <- seq(1790,2010,10) # decennial census
r.guess <- (log(usa[15])-log(usa[1]))/140
k.guess <- usa[15] #1930 US population
par <- c(r.guess,k.guess) #guesses for parameters. necessary to start `optim()`

census1930 <- cbind(usa[1:15], seq(0,140,by=10))
usa1930.fit <- optim(par,fit.logistic,y=census1930)
usa1930.fit$par

#plot predictions against observed
logistic.int <- expression(n0 * exp(r*t)/(1 + n0 * (exp(r*t) - 1)/k))
r <- usa1930.fit$par[1]
k <- usa1930.fit$par[2]
n0 <- usa[1]
t <- seq(0,220,by=10)
plot(seq(1790,2010,by=10), usa, type="n", xlab="Year", 
     ylab="Total Population Size (Millions)")
lines(seq(1790,2010,by=10), eval(logistic.int), lwd=2, col=grey(0.85))
points(seq(1790,2010,by=10),usa, pch=16, col="black")
abline(v=1930, lty=3, col="red") #point at which it deviates.




#Maximum likelihood-----
plague <- read.table("https://web.stanford.edu/class/earthsys214/data/us.plague.txt", header=TRUE)
names(plague)
plague <- plague %>% mutate(time = year - year[1] + 1)
head(plague)
table(plague$cases)
range(plague$cases)
var(plague$cases)/mean(plague$cases) #not 1, so poisson prob wouldn't be good.
plot(plague$time, plague$cases, type="l", xlab="Year", ylab="Plague Cases")

#estimating the "distribution of plague cases" (??) with a poisson model. 
f <- function(x,lambda) -sum(log(dpois(x,lambda))) 
aaa <- optimize(f, c(0,20), x=plague$cases)
aaa #minimum = estimated mean (lambda), objective = -loglikelihood

ll <- rep(0,20)
for(i in 1:20) ll[i] <- -sum(log(dpois(plague$cases,i)))
plot(1:20, ll, type="l", xlab="year", ylab="negative log-likelihood")
abline(h=aaa$objective, lty=3)


#fitting simple epidemic with ML-----
#numerically solve SIR model with `deSolve`

library(deSolve)

#define sir system of equations
sir <- function(t, x, parms){
  S <- x[1]
  I <- x[2]
  R <- x[3]
  with(as.list(parms),{
    dS = -beta*S*I
    dI = beta*S*I - nu*I
    dR = nu*I
    res <- c(dS, dI, dR)
    list(res)
  })
}

#shows how to simulate data----
#solve the equations with `ode` (wrapper of lsoda function)
#define parameters to simulate data
N <- 10000
parms <- c(N=N, beta=.0001, nu = 1/7)
times <- seq(0,30, .1)
x0 <- c(N, 1, 0) #size of initial S, I, R populations
stateMatrix <- ode(y=x0, times, sir, parms)
stateMatrix <- as.data.frame(stateMatrix) 
names(stateMatrix) <- c('time', 'S', 'I', 'R')

stateMatrix %>% pivot_longer(cols = S:R) %>% 
  ggplot(., aes(time, value, color = name)) +
  geom_line()


#now fit data to SIR model----
bombay <- c(0, 4, 10, 15, 18, 21, 31, 51, 53, 97, 125, 183, 292, 390, 448,
            641, 771, 701, 696, 867, 925, 801, 580, 409, 351, 210, 113, 65, 
            52, 51, 39, 33)
cumbombay <- cumsum(bombay)
weeks <- 0:31
plot(weeks, cumbombay, pch=16, xlab="Weeks", ylab="Cumulative Deaths")


require(bbmle)

# likelihood function
sirLL <- function(lbeta, lnu, logN, logI0) {
  parms <- c(beta=plogis(lbeta), nu=plogis(lnu)) #assume parameters between 0-1. plogis=inv_logit transformation
  x0 <- c(S=exp(logN), I=exp(logI0), R=0)
  out <- ode(y=x0, weeks, sir, parms)
  SD <- sqrt(sum( (cumbombay-out[,4])^2)/length(weeks) )
  -sum(dnorm(cumbombay, mean=out[,4], sd=SD, log=TRUE))
}
# minimize negative-log-likelihood
fit <- mle2(sirLL, 
            start=list(lbeta=qlogis(1e-5), #start=guessed values to get the optim started
                       lnu=qlogis(.2), 
                       logN=log(1e6), logI0=log(1) ),  
            method="Nelder-Mead",
            control=list(maxit=1E5,trace=2), #trace=2 see progress of output
            trace=FALSE)

summary(fit)

#estimated parameters/data on outcome scale (backtransformed): beta, nu, N, I0
theta <- as.numeric(c(plogis(coef(fit)[1:2]),
                      exp(coef(fit)[3:4])) )
theta


#plot predictions against observed
parms <- c(beta=theta[1], nu = theta[2])
times <- seq(0,30,0.1)
x0 <- c(theta[3],theta[4],0)
stateMatrix1 <- ode(y=x0, times, sir, parms)
colnames(stateMatrix1) <- c("time","S","I","R")
plot(stateMatrix1[,"time"], stateMatrix1[,"R"], type="l", lwd=2, 
     xaxs="i", xlab="Time", ylab="Cumulative Deaths")
points(weeks, cumbombay, pch=16, col="red")

#plot all the states
stateMatrix1 <- as.data.frame(stateMatrix1)
stateMatrix1 %>% pivot_longer(S:R) %>% 
  ggplot(., aes(time, value, color = name)) +
  geom_line()

parms[1]/parms[2] #R0, right?


#fit the data again, but keeping N and I0 fixed, as if you knew how large the pops were.
fit2 <- mle2(sirLL, 
             start=as.list(coef(fit)),
             fixed=list(logN=coef(fit)[3], 
                        logI0=coef(fit)[4]), 
             method="Nelder-Mead",
             control=list(maxit=1E5,trace=2),
             trace=TRUE)
summary(fit2)
inv_logit(coef(fit2)[1:2]) 



#try to fit mesocosm simulation data??
S_hi <- c(214, 214, 214, 214, 214, 211, 210, 210, 206, 196, 184, 166, 146, 124, 108, 93, 82, 69, 63, 54, 47, 44, 44, 41, 41)
S_med <- c(214,214,214,214,214,213,211,210,208,205,197,194,192,187,184,181,177,175,172,169,167,165,162,160,160)
S_lo <- rep(214, 25)
time <- 1:25
tmp_hi <- tibble(time, S = S_hi) 
tmp_med <- tibble(time, S = S_med) 
tmp_lo <- tibble(time, S = S_lo) 
plot(tmp_hi$time, (tmp_hi$S),  type = 'o', col = 'red3')
lines(tmp_hi$time, tmp_med$S, type = 'o', col = 'yellow4')
lines(tmp_hi$time, tmp_lo$S, type = 'o', col = 'blue4')

tmp2 <- tmp_hi

# likelihood function2. data=S, not R
sirLL2 <- function(lbeta, lnu, logN, logI0) {
  parms <- c(beta=exp(lbeta), nu=exp(lnu)) #assume parameters >0
  x0 <- c(S=exp(logN), I=exp(logI0), R=0)
  out <- ode(y=x0, tmp2$time, sir, parms)
  SD <- sqrt(sum( (tmp2$S-out[,2])^2)/length(tmp2$time) ) #loss function compares predictions against S data (out[,2])
  -sum(dnorm(tmp2$S, mean=out[,2], sd=SD, log=TRUE))
}

fit3 <- mle2(sirLL2, 
             start=list(lbeta = log(.005), 
                        lnu = log(.005), 
                        logN = log(214), logI0 = log(2)),
             fixed=list(logN=log(214)), 
             method="Nelder-Mead",
             control=list(maxit=1E5,trace=2),
             trace=TRUE)
summary(fit3)
theta <- as.numeric(exp(coef(fit3)))
theta
#plot predictions
parms <- c(beta=theta[1], nu = theta[2])
times <- seq(0,30,0.1)
x0 <- c(theta[3],theta[4],0)
stateMatrix1 <- ode(y=x0, times, sir, parms)
colnames(stateMatrix1) <- c("time","S","I","R")

#plot all the states
stateMatrix1 <- as.data.frame(stateMatrix1)
stateMatrix1 %>% pivot_longer(S:R) %>% 
  ggplot(., aes(time, value)) +
  geom_line(aes(color = name)) +
  geom_point(data = rename(tmp2, value = S), aes(time, value))


theta[1]/theta[2] #med=0.005399672, high=0.009804623, low=1.091307e-22


