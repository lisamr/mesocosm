#comparing pI with polynomial linear model
library(dplyr)
library(lme4)
library(merTools)
library(ggplot2)

res <- read.csv("simulation/outputs/resplot3B3.csv")
#reduce the number of reps to just 10
set.seed(1)
s <- sample(unique(res$rep), 20)

res2 <- res %>%
  filter(species == "tot", rep %in% s) 
head(res2)

#####
#Models with NO random effects. doesn't treat communities as nested sets, which I don't think I want.

m1 <- lm(data = res2, pI ~ (rich + I(rich^2) + I(rich^3))*dens*rand)
m2 <- lm(data = res2, pI ~ (rich + I(rich^2))*dens*rand)
m3 <- lm(data = res2, pI ~ (rich)*dens*rand)

#try splines
library(splines)
fit<-lm(pI ~ bs(rich, df=3)*dens*rand,data = res2 )

anova(m2, m1) #m1 more likely
AIC(m1,m1.2, m2, m3, fit) #agrees with AIC
summary(m1)
# save predictions of the model in the new data frame 
# together with variable you want to plot against

#make new data frame to predict response from the model across full range of richness values
N=30
rich = seq(min(res2$rich), max(res2$rich), length.out=N)
d=unique(res2$dens)
r=unique(res2$rand)
rep=s
tmp <- expand.grid(rep=rep, rich=rich, dens=d, rand=r)
tmp$rep2 <- interaction(tmp$dens, tmp$rand, tmp$rep)
head(tmp)

#predict values based on model
pred <- predict(m1,tmp, se.fit = T, interval = "prediction")
tmp2 <- cbind(tmp, pred$fit)

#plot PI of model
ggplot(res2) +
  geom_point(aes((rich), pI, alpha = .2)) +
  geom_line(aes((rich), pI, alpha = .2, group=rep)) +
  geom_line(data=tmp2, aes(rich, fit), color = "blue") +
  geom_ribbon(data=tmp2, aes(rich, ymin = lwr, ymax = upr), fill = "blue", alpha = 0.2) +
  facet_wrap(~dens + rand)

#errors are not normally distributed
res3 <- cbind(res2, fit=fitted.values(m1), resid = m1$res)
hist(res3$resid)
ggplot(res3, aes(resid)) +
  geom_histogram()+
  facet_wrap(~dens + rand)

#display the fit
plot(density(res3$fit), ylim=c(0,4), col="red")
lines(density(res3$pI))
plot(res3$pI ~ res3$fit)+abline(0,1)

#just throw on a smoother
ggplot(res3, aes(jitter(rich), pI, alpha = .2)) +
  geom_point() +
  #geom_boxplot(aes(rich, pI, group = rich)) + 
  #geom_smooth() +
  geom_smooth(method="lm", formula = y ~ poly(x, 3))+
  #geom_smooth(method="loess",fullrange=F)+
  facet_wrap(~dens + rand)

####################################
####################################
#run the model with rep as a random effect?
head(res2)
res2 <- res2 %>% mutate(rep2=interaction(dens, rand, rep))
m4 <- lmer(data = res2, pI ~ poly(rich, 3)*dens*rand + (1|rep2))
m5 <- lmer(data = res2, pI ~ poly(rich, 3)*dens*rand + (rich|rep2))
m6 <- lmer(data = res2, pI ~ poly(rich, 3)*dens*rand + (dens|rep2))#errors
m7 <- lmer(data = res2, pI ~ poly(rich, 3)*dens*rand + (dens+rich|rep2))
m8 <- lmer(data = res2, pI ~ poly(rich, 3)*dens*rand + (dens*rich|rep2))#errors
m9 <- lmer(data = res2, pI ~ poly(rich, 3)*dens*rand + (rand+rich|rep2))
m10 <- lmer(data = res2, pI ~ poly(rich, 3)*dens*rand + (rand+dens|rep2))
m5.1 <- lmer(data = res2, pI ~ poly(rich, 3)*dens*rand + (poly(rich, 2)|rep2))
m5.2 <- lmer(data = res2, pI ~ poly(rich, 3)*dens*rand + (poly(rich, 3)|rep2)) #errors

AIC(m4, m5, m7, m9, m10)
anova(m7, m9)

#fit. purple or red is closest, but doesn't look that different.
plot(res2$pI~fitted.values(m9))+abline(0,1)
plot(density(res2$pI), ylim=c(0,4))
lines(density(fitted.values(m9)), col="grey") #random slope of dens + rich
lines(density(fitted.values(m7)), col="purple") #random slope of dens + rich
lines(density(fitted.values(m5.1)), col="red") #random slope with respect to rich 
lines(density(fitted.values(m5)), col="green") #random slope with respect to rich
lines(density(fitted.values(m1)), col="blue") #no random effects

#predict using the predict function, just to compare to the bootstrapped method
pred <- predict(m9, tmp)
tmp2 <- cbind(tmp, pred)
head(tmp2)

#bootstrap method. gonna get mu and 95% prediction interval
predFun <- function(fit) {
  predict(fit,tmp)
}
PI <- bootMer(m9,nsim=100,FUN=predFun,seed=101) #warning: failed to converge
dim(tmp2) #1200 rows
PI$t %>% str #100 rows (sims) and 1200 columns (data points)
#get metrics: mu, 95% CI of mu, 95% PI of values across the 100 simulations
mu <- apply(PI$t, 2, mean)
#mu.lwr <- apply(PI$t, 2, function(x) mean(x) - sd(x))
#mu.upr <- apply(PI$t, 2, function(x) mean(x) + sd(x))
mu.lwr <- apply(PI$t, 2, function(x) mean(x) - 1.96*sd(x)/sqrt(length(x)))
mu.upr <- apply(PI$t, 2, function(x) mean(x) + 1.96*sd(x)/sqrt(length(x)))
lwr <- apply(PI$t, 2, function(x) quantile(x, .025))
upr <- apply(PI$t, 2, function(x) quantile(x, .975))
tmp2 <- cbind(tmp, mu,mu.lwr, mu.upr, lwr, upr)
head(tmp2)
plot(m9) #something wierd going on at the ends.

ggplot(tmp2, aes(rich, pred, group=rep))+
  geom_point(data=res2,color="blue", aes( (rich), pI, alpha = .2))+
  geom_line(data=res2,color="blue", aes( (rich), pI, alpha = .2))+
  geom_line(color="blue", aes(rich, mu, alpha = .2))+
  geom_ribbon(aes(rich, ymin = mu.lwr, ymax = mu.upr), alpha = 0.9) +
  geom_ribbon(aes(rich, ymin = lwr, ymax = upr), alpha = 0.1) +
  facet_wrap(~dens+rand)
#these prediction intervals are off. the intercepts of mean lines should line up with the data. 

#I just predicted the fit for each replicate. Do I want the group mean, aka. the mean of the means?
tmp3 <- tmp2 %>% group_by(rich, dens, rand) %>% summarise(mu=mean(mu), mu.lwr=mean(mu.lwr), mu.upr=mean(mu.upr), lwr=mean(lwr), upr=mean(upr))
ggplot(tmp3, aes())+
  geom_point(data=res2,color="blue", aes(rich, pI, alpha = .2))+
  geom_line(data=res2,color="blue", aes( rich, pI, group=rep, alpha = .2))+
  geom_line(color="blue", aes(rich, mu, alpha = .2))+
  geom_ribbon(aes(rich, ymin = mu.lwr, ymax = mu.upr), alpha = 0.9) +
  geom_ribbon(aes(rich, ymin = lwr, ymax = upr), alpha = 0.5) +
  facet_wrap(~dens+rand)

#what happens when I predict using 1 random rep?
pred <- predict(m9, filter(tmp, rep==3))
tmp3 <- cbind(filter(tmp, rep==3), pred)
predFun <- function(fit) {
  predict(fit,filter(tmp, rep==3))
}
PI2 <- bootMer(m9,nsim=100,FUN=predFun,seed=101) #warning: failed to converge
mu <- apply(PI2$t, 2, mean)
#mu.lwr <- apply(PI2$t, 2, function(x) mean(x) - sd(x))
#mu.upr <- apply(PI2$t, 2, function(x) mean(x) + sd(x))
mu.lwr <- apply(PI2$t, 2, function(x) mean(x) - 1.96*sd(x)/sqrt(length(x)))
mu.upr <- apply(PI2$t, 2, function(x) mean(x) + 1.96*sd(x)/sqrt(length(x)))
lwr <- apply(PI2$t, 2, function(x) quantile(x, .025))
upr <- apply(PI2$t, 2, function(x) quantile(x, .975))
tmp3 <- cbind(filter(tmp, rep==3), mu,mu.lwr, mu.upr, lwr, upr)
head(tmp3)

ggplot(tmp3, aes(rich, mu))+
  geom_point(data=res2, aes(rich, pI, alpha = .2))+
  geom_line(data=res2, aes(rich, pI, alpha = .2, group=rep))+
  geom_line(color="blue", aes(rich, mu, alpha = .2))+
  geom_ribbon(aes(rich, ymin = mu.lwr, ymax = mu.upr),fill="blue", alpha = 0.3) +
  geom_ribbon(aes(rich, ymin = lwr, ymax = upr),fill="blue", alpha = 0.2) +
  facet_wrap(~dens+rand)

#model predictions obviously off.