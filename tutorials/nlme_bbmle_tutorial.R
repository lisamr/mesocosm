#load libraries
library(nlme)
#first try effect of treatment on logistic growth
Ks <- c(100,200,150) 
n0 <- c(5,5,6)
r <- c(0.15,0.2,0.15)
time <- 1:50
#this function returns population dynamics following
#a logistic curves
logF <- function(time,K,n0,r){
  d <- K * n0 * exp(r*time) / (K + n0 * (exp(r*time) - 1))
  return(d)
}
#simulate some data
dat <- data.frame(Treatment=character(),Time=numeric(),
                  Abundance=numeric())
for(i in 1:3){
  Ab <- logF(time = time,K=Ks[i],n0=n0[i],r=r[i])
  tmp <- data.frame(Treatment=paste0("T",i),Time=time,
                    Abundance=Ab+rnorm(time,0,5))
  #note that random deviates were added to the simulated 
  #population density values
  dat <-rbind(dat,tmp)
}
#the formula for the models
lF<-formula(Abundance~K*n0*exp(r*Time)/(K+n0*(exp(r*Time)-1)) | Treatment)

#fit the model
(m <- nlsList(lF,data=dat,start=list(K=150,n0=10,r=0.5)))


#continuous predictor-----

#use `bbmle::mle2``


#load libraries
library(bbmle)

#parameter for simulation
K <- 200
n0 <- 5
#the gradient in temperature
Temp <- ceiling(seq(0,20,length=10))


#simulate some data
mm <- sapply(Temp,function(x){ 
  rpois(50,logF(time = 1:50,K = K,n0=n0,r = 0.05 + 0.01 * x))})
#sample data from a poisson distribution with lambda parameter equal
#to the expected value, note that we do not provide one value for the
#parameter r but a vector of values representing the effect of temperature
#on r

#some reshaping
datT <- reshape2::melt(mm)
names(datT) <- c("Time","Temp","Abund")
datT$Temp <- Temp[datT$Temp]
head(datT)

#fit the model. Looks like just r is predicted. 
#Parameters estimated with a linear model are designated in the parameters argument. 
(mll <- mle2(Abund~dpois(K*n0*exp(r*Time)/(K+n0*(exp(r*Time)-1))),
             data=datT,parameters = list(r~Temp),
             start=list(K=100,n0=10,r=0.05)))
mll
summary(mll)

#check out predictors
datT$pred <- predict(mll)
ggplot(datT,aes(x=Time,y=Abund,color=factor(Temp)))+
  geom_point()+geom_line(aes(y=pred))


#repeat with a binomial likelihood
datT$N <- 300L
datT1 <- datT %>% filter(Temp==18)
head(datT1)
str(datT1)
sort(datT1$Abund/datT1$N)
ggplot(datT1, aes(Time, Abund)) +
  geom_point()

log_curve <- function(K, n0, r, Time){
  K*n0 / (n0 + (K - n0)*exp(-r*Time))
}
log_curve(1, 0, .2, 0:50)

mbinom <- mle2(Abund ~ dbinom(
  size = N, 
  prob = log_curve(K, n0, r, Time)
),
data=datT1,
start=list(K=.9,n0=.1,r=.2),
lower=c(K=1e-6, r=1e-6, n0=1e-6 ),
upper = c(K=1-1e-6), method = 'L-BFGS-B')
summary(mbinom)
datT1$pred2 <- predict(mbinom, newdata = datT1)

ggplot(datT1, aes(Time, Abund)) +
  geom_point() +
  geom_line(aes(y = pred2))


#non-monotonic relationship between K (pop max) and temperature----
#simulate some data. K and r now have functions attached to them. 
mm <- sapply(Temp,function(x){
  rpois(50,logF(time = 1:50,K = 100+20*x-x^2,n0=n0,r = 0.05 + 0.01 * x))})
#note that this time both K and r are given as vectors to represent their
#relation with temperature

#some reshaping
datT <- reshape::melt(mm)
names(datT) <- c("Time","Temp","Abund")
datT$Temp <- Temp[datT$Temp]
head(datT)

#the negative log-likelihood function
LL <- function(k0,k1,k2,n0,r0,r1){
  K <- k0 + k1*Temp + k2*Temp**2
  r <- r0 + r1*Temp
  lbd <- K*n0*exp(r*Time)/(K+n0*(exp(r*Time)-1))
  -sum(dpois(Abund,lbd,log=TRUE)) #calculate -LL. mle2 finds pars that minimize this. 
}

#notice response variable isn't predicted in the formula because it's included as an argument in LL
#I think I like this method better because I can clearly see how the model is built. 
(mll2 <- mle2(LL,data=datT,start=list(k0=100,k1=1,k2=-0.1,n0=10,r0=0.1,r1=0.1)))
mll2
summary(mll2)


cc <- coef(mll2)
#sorry for the awful coding, to get the model predicted values we need
#to code the relation by hand which is rather clumsy ...
datT$pred <- reshape::melt(sapply(Temp,function(x) {logF(time=1:50,K=cc[1]+cc[2]*x+cc[3]*x**2,n0=cc[4],r=cc[5]+cc[6]*x)}))[,3]
ggplot(datT,aes(x=Time,y=Abund,color=factor(Temp)))+
  geom_point()+geom_line(aes(y=pred))
