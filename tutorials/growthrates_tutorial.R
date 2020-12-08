#https://cran.r-project.org/web/packages/growthrates/vignettes/Introduction.html

#install.packages('growthrates')
library(growthrates)
library(tidyverse)
theme_set(theme_classic())

data(bactgrowth)
head(bactgrowth)
str(bactgrowth)

ggplot(bactgrowth, aes(time, value, group = interaction(strain, replicate))) +
  geom_line(aes(color = strain)) +
  facet_wrap( ~conc) +
  scale_color_viridis_d(option = 'D')

#fit single subsets to various models----

splitted.data <- multisplit(bactgrowth, c("strain", "conc", "replicate")) #splits data by replicate/treatment. 
dat <- splitted.data[[1]]
dat


#linear relationship----
fit <- fit_easylinear(dat$time, dat$value)
summary(fit)
coef(fit)#exp growth pars. "mumax" must be intrinsic growth rate.
#y0     y0_lm     mumax       lag 
#0.0180000 0.0123482 0.2048985 1.8392607 
plot(fit, log='y') #fits a linear function to the first part of the curve.


#logistic curve----
#dN/dt = mumax*N(1-N/K)
p     <- c(y0 = 0.01, mumax = 0.2, K = 0.1) #mean starting pars
lower <- c(y0 = 1e-6, mumax = 0,   K = 0) #lower bounds of pars
upper <- c(y0 = 0.05, mumax = 5,   K = 0.5)#upper bounds of pars

fit1 <- fit_growthmodel(FUN = grow_logistic, p = p, time = dat$time, y = dat$value,lower = lower, upper = upper)
summary(fit1)
#      Estimate Std. Error t value Pr(>|t|)    
#y0    0.017483   0.001581   11.06 9.98e-12 ***
#mumax 0.200069   0.013979   14.31 2.10e-14 ***
#K     0.099626   0.001850   53.87  < 2e-16 ***

#two-step logistic curve-----
#dy_i/dt = -kw*y_i
#dy_a/dt = kw*y_i + mumax*y_a(1-(y_a+y_i)/K)
#I think this better represents the dynamics in the trays, cuz step 1 represents primary transmission and step two is secondary transmission. 
p     <- c(yi = 0.02, ya = 0.001, kw = 0.1, mumax = 0.2, K = 0.1)
lower <- c(yi = 1e-6, ya = 1e-6, kw = 0,    mumax = 0,   K = 0)
upper <- c(yi = 0.05, ya = 0.05, kw = 10,   mumax = 5,   K = 0.5)

fit2 <- fit_growthmodel(FUN = grow_twostep, p = p, time = dat$time, y = dat$value,
                        lower = lower, upper = upper)
summary(fit2)
coef(fit2)
#yi          ya          kw       mumax           K 
#0.014595184 0.003855074 5.711012173 0.194235013 0.100192751 


#in case you want to not estimate all the parameters (like the starting pop sizes)
#specify which pars you want to estimate with `which`
fit3 <- fit_growthmodel(FUN = grow_twostep, p = p, time = dat$time, y = dat$value,
                        lower = lower, upper = upper, which = c("kw", "mumax", "K"))
summary(fit3)
coef(fit3) #yi and ya were determined by the values supplied in `p` (fixed at whatever you say they will be)

#compare fits
par(mfrow = c(1,3))
plot(fit1)
plot(fit2)
plot(fit3)



#nonparametric smoothing splines-----

#way to estimate max growth.
dat <- splitted.data[[2]]
time <- dat$time
y    <- dat$value

## automatic smoothing with cv
res <- fit_spline(time, y)

par(mfrow = c(1, 2))
plot(res, log = "y")
plot(res)
coef(res)


#fittingn multiple datasets-----

#use `all_growthmodels` to fit data from multiple trials and select a function to fit the model, ex/ `grow_logistic(time, parms)`. Include covariates in the formula after the base function.
#lets use the parametric models. start with the baranyi growth model?


#first attempt: fit all parameters in one model:
## initial parameters and box constraints
p   <- c(y0 = 0.03, mumax = .1, K = 0.1, h0 = 1)
lower   <- c(y0 = 0.001, mumax = 1e-2, K = 0.005, h0 = 0)
upper   <- c(y0 = 0.1,   mumax = 1,    K = 0.5,   h0 = 10)

## fit growth models to all data using log transformed residuals
many_baranyi1 <- all_growthmodels(
  value ~ grow_baranyi(time, parms) | strain + conc + replicate,
  data = bactgrowth,
  p = p, lower = lower, upper = upper,
  transform = "log", ncores = 2)
summary(many_baranyi1) #spits out independent? coefs for each trial. 
coef(many_baranyi1)

#second attempt: fix h0.
## use coefficients of first fit as new initial parameters
pp   <- coef(many_baranyi1)
## but set h0 to a (arbitrary) fixed value
pp[, "h0"] <- 0.65 #mean(pp[, "h0"] ) = .955
## re-fit models
many_baranyi2 <- all_growthmodels(
  value ~ grow_baranyi(time, parms) | strain + conc + replicate,
  data = bactgrowth,
  p = pp, lower = lower, upper = upper,
  which = c("y0", "mumax", "K"), transform = "log", ncores = 2)

par(mfrow = c(12, 6))
par(mar = c(2.5, 4, 2, 1))
plot(many_baranyi2)



many_baranyi2_res <- results(many_baranyi2)
xyplot(mumax ~ log(conc+1)|strain, data = many_baranyi2_res, layout = c(3, 1 ))





#try fitting the data again but with the logistic funciton
## fit growth models to all data using log transformed residuals
p     <- c(y0 = 0.01, mumax = 0.2, K = 0.1) #mean starting pars
lower <- c(y0 = 1e-6, mumax = 0,   K = 0) #lower bounds of pars
upper <- c(y0 = 0.05, mumax = 5,   K = 0.5)#upper bounds of pars
many_logistic <- all_growthmodels(
  value ~ grow_logistic(time, parms) | strain + conc,
  data = bactgrowth,
  p = p, lower = lower, upper = upper, ncores = 2)
summary(many_logistic) #spits out independent? coefs for each trial. 
coef(many_logistic)
plot(many_logistic)

predict(many_logistic)
