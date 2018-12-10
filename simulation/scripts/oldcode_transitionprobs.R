#GO DOWN TO THE BOTTOM FOR THE STREAMLINED CODE

#define parameter values from table 1 in the paper
pars <- data.frame(par=c("ap", "sp", "dp", "lp", "np", "tp", "as", "ss", "ds", "ls", "ts", "g"), mu=c(1.09, .277, .55, 0, .49, 4.898, 14.9, .00099, .36, .0524, 3.73, 3.89), se=c(.2, .027, .2, NA, .15, .082, 4, .00022, .19, .0058, .15 ,5))

########################################
#r(x,t)
########################################
#functions for primary transmission
#distance decay function
phi.p <- function(x){ 
  ap=pars$mu[pars$par=="ap"]
  dp=pars$mu[pars$par=="dp"]
  tp=pars$mu[pars$par=="tp"]
  np=pars$mu[pars$par=="np"]
  sp=pars$mu[pars$par=="sp"]
  Cp=ap*dp*exp(dp*tp)
  
  Cp*(1+np*x)*exp(-sp*x)
}
#time decay function
psi.p <- function(t){ 
  dp=pars$mu[pars$par=="dp"]
  exp(-dp*t)
}
#step function for delay in infection
H.p <- function(x, t){
  tp=pars$mu[pars$par=="tp"]
  lp=pars$mu[pars$par=="lp"]
  y <- t-tp-lp*x
  ifelse(y>0, 1, 0)
} 

#functions for primary transmission
phi.s <- function(x){
  as=pars$mu[pars$par=="as"]
  ds=pars$mu[pars$par=="ds"]
  ts=pars$mu[pars$par=="ts"]
  Cs=as*ds*exp(ds*ts)
  
  Cs*exp(-ss*(x^2))
}

psi.s <- function(t){
  ds=pars$mu[pars$par=="ds"]
  exp(-ds*t)
}

H.s <- function(x, t){
  ts=pars$mu[pars$par=="ts"]
  ls=pars$mu[pars$par=="ls"]
  y <- t-ts-ls*x
  ifelse(y>0, 1, 0)
} 

rp <- function(x, t){
  phi.p(x) * psi.p(t) * H.p(x, t)
}
rs <- function(x, t){
  phi.s(x) * psi.s(t) * H.s(x, t)
}
rp(10, 1:10)
rs(1, 1:10)

#values of primary transmission
dist=20
time=20
by.value=1
test.d <- phi.p(seq(1, dist, by =by.value), pars$mu[pars$par=="sp"])
test.t <- psi.p(seq(1, time, by =by.value), pars$mu[pars$par=="dp"])
test.H <- sapply(seq(1, time, by =by.value), function(t) H.p(t, x=seq(1, dist, by =by.value))) #time=columns, dist=rows
test.m <- t(sapply(seq(1, dist,by =by.value), function(x) test.d[x]*test.t))#time=columns, dist=rows
Mp <- test.H*test.m
plot(seq(1, dist, by =by.value), test.d, xlab="distance (mm)", ylab="phi.p")
plot(seq(1, time, by =by.value), test.t, xlab="time (d)", ylab="psi.p")
persp3D(z = Mp, phi=25, theta = 235, xlab="dist", ylab="time", main="primary rates of infection")


#values of secondary transmision
dist.start=17.5
dist=30
time=33
by.d=.1
by.t=1
test.d <- phi.s(seq(dist.start, dist, by =by.d))
test.t <- psi.s(seq(1, time, by =by.t))
test.H <- sapply(seq(1, time, by =by.t), function(t) H.s(t, x=seq(1, dist, by =by.d))) #time=columns, dist=rows
test.m <- t(sapply(seq(1, dist, by =by.d), function(x) test.d[x]*test.t))#time=columns, dist=rows
Ms <- test.H*test.m
plot(seq(dist.start, dist, by =by.d), test.d, xlab="distance (mm)", ylab="phi.p")
plot(seq(1, time, by =by.t), test.t, xlab="time (d)", ylab="psi.p")
persp3D(z = Ms, phi=0, theta = 40, facets = F, xlab="dist", ylab="time", main="secondary rates of infection")
persp3D(z = Ms, phi=0, theta = 135, facets = F, xlab="dist", ylab="time", main="secondary rates of infection")
persp3D(z = Ms, phi=0, theta = 245, facets = F, xlab="dist", ylab="time", main="secondary rates of infection")

#varying magnitude of infection rates by species.
#assumes all parameters decay proportionally, which obviously is a big assumption. 
mag=c(1, .5, .4, .2, 0, 0)

#primary transmission
t=30
dstart=1
dend=10
by.d=1
Mp <- sapply(seq(dstart, dend, by =by.d), function(d) rp(d, 1:t)) %>% t()

#secondary transmission
t=33
dstart=17.5
dend=30
by.d=.1
Ms <- sapply(seq(dstart, dend, by =by.d), function(d) rs(d, 1:t)) %>% t()
Ms <- Ms[,-(1:3)]#chop off the first 3 seconds

#primary rates
Mps <- list(NULL) #time=15 columns, distance=20 rows
for(i in 1:length(mag)){
  Mps[[i]] <- Mp*mag[i]
}

#secondary rates
Mss <- list(NULL) #time=15 columns, distance=50 rows
for(i in 1:length(mag)){
  Mss[[i]] <- Ms*mag[i]
}

persp3D(z = Mss[[1]], phi=25, theta = 120, xlab="dist", ylab="time", main="secondary rates of infection")

#turning it into a dataframe
#secondary
d=seq(17.5, 30.0, by=.1)
t=30
MS <- function(i){
  data.frame(dist=rep(d, t), time=rep(1:t, each=length(d)), species=i, trans="sec", beta=as.vector(Mss[[i]])) 
}
M1 <- rbind(MS(1), MS(2), MS(3), MS(4), MS(5), MS(6))

#primary
d=10
MP <- function(i){
  data.frame(dist=rep(1:d, t), time=rep(1:t, each=d), species=i, trans="prim", beta=as.vector(Mps[[i]])) 
}
M2 <- rbind(MP(1), MP(2), MP(3), MP(4), MP(5), MP(6))

M3 <- rbind(M1, M2)
summary(M3)

head(M3)
ggplot(filter(M3, trans=="sec", dist==20), aes(time, beta, group=species)) +
  geom_point()+
  geom_line()
ggplot(filter(M3, trans=="prim", dist==2), aes(time, beta, group=species)) +
  geom_point()+
  geom_line()

write.csv(M3, "simulation/outputs/rates.csv", row.names = F)


####
########################################
#streamlining code
########################################

#define parameter values from table 1 in the paper
pars <- data.frame(par=c("ap", "sp", "dp", "lp", "np", "tp", "as", "ss", "ds", "ls", "ts", "g"), mu=c(1.09, .277, .55, 0, .49, 4.898, 14.9, .00099, .36, .0524, 3.73, 3.89), se=c(.2, .027, .2, NA, .15, .082, 4, .00022, .19, .0058, .15 ,5))

#define functions
#r(x,t)

#functions for primary transmission
#distance decay function
phi.p <- function(x){ 
  ap=pars$mu[pars$par=="ap"]
  dp=pars$mu[pars$par=="dp"]
  tp=pars$mu[pars$par=="tp"]
  np=pars$mu[pars$par=="np"]
  sp=pars$mu[pars$par=="sp"]
  Cp=ap*dp*exp(dp*tp)
  
  Cp*(1+np*x)*exp(-sp*x)
}
#time decay function
psi.p <- function(t){ 
  dp=pars$mu[pars$par=="dp"]
  exp(-dp*t)
}
#step function for delay in infection
H.p <- function(x, t){
  tp=pars$mu[pars$par=="tp"]
  lp=pars$mu[pars$par=="lp"]
  y <- t-tp-lp*x
  ifelse(y>0, 1, 0)
} 

#functions for primary transmission
phi.s <- function(x){
  as=pars$mu[pars$par=="as"]
  ds=pars$mu[pars$par=="ds"]
  ts=pars$mu[pars$par=="ts"]
  ss=pars$mu[pars$par=="ss"]
  Cs=as*ds*exp(ds*ts)
  
  Cs*exp(-ss*(x^2))
}

psi.s <- function(t){
  ds=pars$mu[pars$par=="ds"]
  exp(-ds*t)
}

H.s <- function(x, t){
  ts=pars$mu[pars$par=="ts"]
  ls=pars$mu[pars$par=="ls"]
  y <- t-ts-ls*x
  ifelse(y>0, 1, 0)
} 

rp <- function(x, t){
  phi.p(x) * psi.p(t) * H.p(x, t)
}
rs <- function(x, t){
  phi.s(x) * psi.s(t) * H.s(x, t)
}

#varying magnitude of infection rates by species.
#assumes all parameters decay proportionally, which obviously is a big assumption. 
mag=c(1, .3, .2, .1, 0, 0)

#primary transmission
t=30
dstart=1
dend=10
by.d=1
Mp <- sapply(seq(dstart, dend, by =by.d), function(d) rp(d, 1:t)) %>% t()

#secondary transmission
t=33
dstart=17.5
dend=30
by.d=.1
Ms <- sapply(seq(dstart, dend, by =by.d), function(d) rs(d, 1:t)) %>% t()
Ms <- Ms[,-(1:3)]#chop off the first 3 seconds
#standardize to the maximum number
Ms <- Ms/max(Ms)

#primary rates
Mps <- list(NULL) #time=15 columns, distance=20 rows
for(i in 1:length(mag)){
  Mps[[i]] <- Mp*mag[i]
}

#secondary rates
Mss <- list(NULL) #time=15 columns, distance=50 rows
for(i in 1:length(mag)){
  Mss[[i]] <- Ms*mag[i]
}

#turning it into a dataframe
#secondary
d=seq(17.5, 30.0, by=.1)
t=30
MS <- function(i){
  data.frame(dist=rep(d, t), time=rep(1:t, each=length(d)), species=i, trans="sec", beta=as.vector(Mss[[i]])) 
}
M1 <- rbind(MS(1), MS(2), MS(3), MS(4), MS(5), MS(6))

#primary
d=10
MP <- function(i){
  data.frame(dist=rep(1:d, t), time=rep(1:t, each=d), species=i, trans="prim", beta=as.vector(Mps[[i]])) 
}
M2 <- rbind(MP(1), MP(2), MP(3), MP(4), MP(5), MP(6))

M3 <- rbind(M1, M2)
summary(M3)

head(M3)
ggplot(filter(M3, trans=="sec", dist==20), aes(time, beta, group=species)) +
  geom_point()+
  geom_line()
ggplot(filter(M3, trans=="sec", time==3), aes(dist, beta, group=species)) +
  geom_point()+
  geom_line()
ggplot(filter(M3, trans=="prim", dist==2), aes(time, beta, group=species)) +
  geom_point()+
  geom_line()

write.csv(M3, "simulation/outputs/rates.csv", row.names = F)

