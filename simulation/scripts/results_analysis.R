#comparing pI with polynomial linear model

res <- read.csv("simulation/outputs/resplot3B3.csv")
res2 <- res %>%
  filter(species == "tot") 

m1 <- lm(data = res2, pI ~ (rich + I(rich^2) + I(rich^3))*dens*rand)
m2 <- lm(data = res2, pI ~ (rich + I(rich^2))*dens*rand)
m3 <- lm(data = res2, pI ~ (rich)*dens*rand)

#try splines
library(splines)
fit<-lm(pI ~ bs(rich, df=3)*dens*rand,data = res2 )

anova(m2, m1) #m1 more likely
AIC(m1, m2, m3, fit) #agrees with AIC
summary(m1)
# save predictions of the model in the new data frame 
# together with variable you want to plot against
res3 <- data.frame(res2, predict(m1,res2, se.fit = T, interval = "confidence"), res=residuals(m1))

ggplot(res3) +
  geom_point(aes(jitter(rich), pI, alpha = .2)) +
  #geom_boxplot(aes(rich, pI, group = rich)) + 
  geom_line(aes(rich, fit.fit), color = "blue") + 
  geom_line(aes(rich, fit.lwr)) + 
  geom_line(aes(rich, fit.upr)) + 
  facet_wrap(~dens + rand)
ggplot(res3) +
  geom_point(aes(jitter(rich), res, alpha=.2)) +
  geom_abline(slope = 0, intercept =  0)+
  facet_wrap(~dens + rand)

#errors are not normally distributed
hist(res3$res)
ggplot(res3, aes(res)) +
  geom_histogram()+
  facet_wrap(~dens + rand)

#just throw on a smoother
ggplot(res3, aes(jitter(rich), pI, alpha = .2)) +
  geom_point() +
  #geom_boxplot(aes(rich, pI, group = rich)) + 
  geom_smooth() +
  facet_wrap(~dens + rand)



