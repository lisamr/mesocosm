#functions from Kleczkowski, Gilligan and Bailey (1997) used to describe the transitions between states, S and Y (Inf+latent)
#I think the r(x,t) values are equal to the rates of new infections per infected plant per time interval.
library(dplyr)
library(ggplot2)
library(plot3D)
#define parameter values from table 1 in the paper
pars <- data.frame(par=c("ap", "sp", "dp", "lp", "np", "tp", "as", "ss", "ds", "ls", "ts", "g"), mu=c(1.09, .277, .55, 0, .49, 4.898, 14.9, .00099, .36, .0524, 3.73, 3.89), se=c(.2, .027, .2, NA, .15, .082, 4, .00022, .19, .0058, .15 ,5))

#I'm going to reconstruct the Y_s(x,t) curves
Yp <- function(x, t){
  ap=pars$mu[pars$par=="ap"]
  np=pars$mu[pars$par=="ap"]
  sp=pars$mu[pars$par=="sp"]
  dp=pars$mu[pars$par=="dp"]
  lp=pars$mu[pars$par=="lp"]
  tp=pars$mu[pars$par=="tp"]
  H <- ifelse((t-tp-lp*x)>0, 1, 0)
  1-exp(-ap*(1+np*x) * exp(-sp*x) * (exp(-dp*lp*x) - exp(-dp*(t-tp))) * H)
}
Ys <- function(x, t){
  as=pars$mu[pars$par=="as"]
  ss=pars$mu[pars$par=="ss"]
  ds=pars$mu[pars$par=="ds"]
  ls=pars$mu[pars$par=="ls"]
  ts=pars$mu[pars$par=="ts"]
  #t=1:10
  #x=1
  ifelse((t-ts-ls*x)<0, 1, 0)
  H <- ifelse((t-ts-ls*x)>0, 1, 0)
  1-exp(-as*exp(-ss*(x^2)) * (exp(-ds*ls*x) - exp(-ds*(t-ts))) * H)
}

MYp <- sapply(1:20, function(x) Yp(x, seq(1,10, .1))) #time=col, dist=row
persp3D(z = (MYp), phi=0, theta = 290, xlab="time", ylab="dist", main="primary rates of infection")
MYs <- sapply(1:50, function(x) Ys(x, seq(1,10, .1))) #time=col, dist=row
persp3D(z = (MYs), phi=0, theta = 290, xlab="time", ylab="dist", main="secondary rates of infection")

#calculate derivate from these curves.
ys <- MYs[,20]
t=seq(1,10, .1)
plot(t, ys)
ys.prime <- diff(ys)/diff(t)
plot(ys.prime)
sum(ys.prime)*.1 

yp <- MYp[,2] #first distance
t=seq(1,10, .1)
plot(t, yp) #cumulative infections over time
yp.prime <- c(NA, c(diff(yp)/diff(t)))
plot(yp.prime) #rate of infections over time 
sum(yp.prime, na.rm = T)*.1 #sum of rates*dt=total infections

########################################
#how does dI/dt relate to beta??
########################################
yp#dI/dt (from primary transmission) at distance 2mm
I <- yp*10
S <- 10-I
plot(1:length(yp), S, ylim=c(0,10), main="S and I after prim. trans")
points(1:length(yp), I)

beta.p <- function(t){
  yp.prime[t]/(S[t]*I[t])
}
plot(1:100, beta.p(1:100), type='o', main="beta.prim")

ys#dI/dt (from primary transmission) at distance 20mm
I <- ys*10
S <- 10-I
plot(1:length(yp), S, ylim=c(0,10), main="S and I after sec. trans")
points(1:length(yp), I)

beta.s <- function(t){
  ys.prime[t]/(S[t]*I[t])
}
plot(1:100, beta.s(1:100), main="beta.prim", type='o')

################################################################
################################################################
#define parameter values from table 1 in the paper
pars <- data.frame(par=c("ap", "sp", "dp", "lp", "np", "tp", "as", "ss", "ds", "ls", "ts", "g"), mu=c(1.09, .277, .55, 0, .49, 4.898, 14.9, .00099, .36, .0524, 3.73, 3.89), se=c(.2, .027, .2, NA, .15, .082, 4, .00022, .19, .0058, .15 ,5))


#make a function to get beta streamlined
BetaP <- function(x, t){
  #I'm going to reconstruct the Y_s(x,t) curves
  Yp <- function(x, t){
    ap=pars$mu[pars$par=="ap"]
    np=pars$mu[pars$par=="ap"]
    sp=pars$mu[pars$par=="sp"]
    dp=pars$mu[pars$par=="dp"]
    lp=pars$mu[pars$par=="lp"]
    tp=pars$mu[pars$par=="tp"]
    H <- ifelse((t-tp-lp*x)>0, 1, 0)
    1-exp(-ap*(1+np*x) * exp(-sp*x) * (exp(-dp*lp*x) - exp(-dp*(t-tp))) * H)
  }
  Y <- Yp(x, t)
  dy <- c(NA, diff(Y)/diff(t))
  I <- Y
  S <- 1-I
  beta <- function(t){
    B <- dy[t]/(S[t]*I[t])
    B[is.na(B)] <- 0
    B
  }
  beta(t)
}

t=30

Bp <- sapply(1:20, function(x) BetaP(x, seq(1,t, by=1))) %>% t() #time=col, dist=row
persp3D(z =Bp, phi=25, facets=F, bty = "b2", theta = 235, ylab="time", xlab="dist", main="beta primary infection")

BetaS <- function(x, t){
  #I'm going to reconstruct the Y_s(x,t) curves
  Ys <- function(x, t){
    as=pars$mu[pars$par=="as"]
    #as=4
    ss=pars$mu[pars$par=="ss"]
    #ss=.000001
    ds=pars$mu[pars$par=="ds"]
    ls=pars$mu[pars$par=="ls"]
    ts=pars$mu[pars$par=="ts"]
    #t=1:10
    #x=1
    ifelse((t-ts-ls*x)<0, 1, 0)
    H <- ifelse((t-ts-ls*x)>0, 1, 0)
    1-exp(-as*exp(-ss*(x^2)) * (exp(-ds*ls*x) - exp(-ds*(t-ts))) * H)
  }
  Y <- Ys(x, t)
  dy <- c(NA, diff(Y*10)/diff(t))
  I <- Y*10
  S <- 10-I
  beta <- function(t){
    B <- dy[t]/(S[t]*I[t])
    B[is.na(B)] <- 0
    B
  }
  beta(t)
  #dy
}
Bs <- sapply(1:50, function(x) BetaS(x, seq(1,t, by=1))) %>% t() #time=col, dist=row

persp3D(z =Bs, phi=0,facets = F, theta = 170, ylab="time", xlab="dist", main="beta secondary infection", bty = "b2")
persp3D(z =Bs, phi=0,facets = F, theta = 300, ylab="time", xlab="dist", main="beta secondary infection",  bty = "b2")
persp3D(z =Bs, phi=0,facets = F, theta = 30, ylab="time", xlab="dist", main="beta secondary infection",  bty = "b2")

########################################################
########################################################
#Make dataframe with time, dist, species, transmissiontype, beta
#do this again for relevant distances: interplanting distances(cm) 1.75 1.90 2.28 3.00
t=33 #(going to chop off the first 3 days because too much lag)
Bp <- sapply(1:20, function(x) BetaP(x, seq(1,t, by=1))) %>% t() #time=col, dist=row
Bp <- Bp[,-(1:3)]
Bs <- sapply(seq(17.5, 30.0, by=.1), function(x) BetaS(x, seq(1,t, by=1))) %>% t() #time=col, dist=row
Bs <- Bs[,-(1:3)]

#varying magnitude of infection rates by species.
#assumes all parameters decay proportionally, which obviously is a big assumption.
mag=c(1, .6, .4, .2, 0, 0)
#primary rates
Bps <- list(NULL) #time=15 columns, distance=20 rows
for(i in 1:length(mag)){
  Bps[[i]] <- Bp*mag[i]
}
#secondary rates
Bss <- list(NULL) #time=15 columns, distance=50 rows
for(i in 1:length(mag)){
  Bss[[i]] <- Bs*mag[i]
}

#persp3D(z =Bss[[1]], phi=0,facets = F, theta = 170, ylab="time", xlab="dist", main="beta secondary infection", bty = "b2")
#persp3D(z =Bss[[1]], phi=0,facets = F, theta = 270, ylab="time", xlab="dist", main="beta secondary infection", bty = "b2")

#Making dataframe
d=seq(17.5, 30.0, by=.1)
t=30
BS <- function(i){
  data.frame(dist=rep(d, t), time=rep(1:t, each=length(d)), species=i, trans="sec", beta=as.vector(Bss[[i]])) 
}
B1 <- rbind(BS(1), BS(2), BS(3), BS(4), BS(5), BS(6))
d=20
BP <- function(i){
  data.frame(dist=rep(1:d, t), time=rep(1:t, each=d), species=i, trans="prim", beta=as.vector(Bps[[i]])) 
}
B2 <- rbind(BP(1), BP(2), BP(3), BP(4), BP(5), BP(6))

B3 <- rbind(B1, B2)
summary(B3)

head(B3)
ggplot(filter(B3, trans=="sec", dist==20), aes(time, beta, group=species)) +
  geom_point()+
  geom_line()
ggplot(filter(B3, trans=="sec", time==3), aes(dist, beta, group=species)) +
  geom_point()+
  geom_line()
ggplot(filter(B3, trans=="prim", dist==20), aes(time, beta, group=species)) +
  geom_point()+
  geom_line()
ggplot(filter(B3, trans=="prim", time==3), aes(dist, beta, group=species)) +
  geom_point()+
  geom_line()

write.csv(B3, "simulation/outputs/betas.csv", row.names = F)


