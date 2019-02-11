#more time steps. how does that affect your results? does it allow you to fit it to a defined non-linear function? does that function have a parameter that measures "slope"?
library(emmeans)

design2 <- read.csv("simulation/outputs/design2.csv")
design2$r[design2$r==2.39] <- 2.4 #for some reason R sucks and hates 2.39. changing it slightly.

design.128 <- design2 %>% filter(d=="add"| d=="sub" & n==128)
design.210 <- design2 %>% filter(d=="add"| d=="sub" & n==200)
design.338 <- design2 %>% filter(d=="add"| d=="sub" & n==338)
design.392 <- design2 %>% filter(d=="add"| d=="sub" & n==392)

#when does disease stop spreading in the additive-deterministic treatment? Run for 40 time steps and get results
test <- results3(f.sim = s2,t = 40, B = B, nreps = 10, density = "add", Random = F, design = design.392)
testr <- test$responses 
ggplot(filter(testr, species=="tot"), aes(time, n.I, color=rich))+
  geom_point()+
  geom_line(aes(group=interaction(rich, rep)))
ggplot(filter(testr, species=="tot"), aes(time, pI))+
  geom_point()+
  geom_line(aes(group=interaction(rich, rep)))+
  geom_smooth(aes(group=rich))

moretime <- fb3(Beta = B, t = 30, n = 10, design = design.392)
moretimer <- moretime$responses 
moretimer.sum <- moretime$response.summary
#write.csv(moretimer, "simulation/outputs/moretime_responses.csv", row.names = F)
#write.csv(moretimer.sum, "simulation/outputs/moretime_response_summary.csv", row.names = F)

ggplot(filter(moretimer, species=="tot"), aes(time, pI))+
  #geom_point()+
  #geom_line(aes(group=interaction(rich, rep)))+
  geom_smooth(aes(group=rich,  color=as.factor(rich)))+
  facet_wrap(~dens+rand)
ggplot(filter(moretimer, species=="tot", time==30), aes(rich, (pI)))+
  geom_point()+
  geom_line(aes(group=interaction( rep)))+
  #geom_smooth(method='lm')+
  facet_wrap(~rand+dens)+
  background_grid(major = "xy", minor = "none")+
  labs(x="richness", y="proportion infected")+
  theme(strip.background = element_rect(fill="white",color="black"))

#another way to look at the interactions of density and assembly order
ggplot(filter(moretimer, species=="tot", time==30), aes(rich, (pI), color=interaction(rand, dens), group=interaction(rand, dens)))+
  geom_point()+
  geom_smooth(method='lm')
ggplot(filter(moretimer, species=="tot", time==30), aes(rich, (n.I), color=interaction(rand, dens), group=interaction(rand, dens)))+
  geom_point()+
  #geom_line(aes(group=interaction(rand, dens, rep)))+
  geom_smooth(method='lm')

#try fitting the data to a linear function and then try non-linear. look over old notes from ABG250.
d <- filter(moretimer, species=="tot", time==30)
m1 <- lm(pI ~ rich, data=d)
m2 <- lm(pI ~ -1+ (dens*rand)/rich, data=d)
m3 <- lm(pI ~ -1+ (dens+rand)*rich, data=d)
m4<- lm(pI ~ rich + dens + rand, data=d)
m5<- lm(pI ~ rich + dens*rand, data=d)
summary(m2)
#ANCOVA. compares the slopes.
#disease incidence
library(nlme)
d <- d %>% mutate(trt = interaction(dens, rand),
                  rep2 = interaction(trt, rep)) 
m1 <- lm(pI ~ -1 + trt/rich, data=d) #richness is nested within treatments
summary(m1)
lst <- lstrends(m1, "trt", var="rich")
pairs(lst)
plot(lst)

#total infections
m3 <- gls(n.I ~ -1+trt/rich, data=d, weights = varPower())
lst2 <- lstrends(m3, "trt", var="rich")#same as m1. SE the same, which it shouldnt be.
pairs(lst2)
plot(lst2)

#standard error doesn't match SE when I subset the data.
lm(n.I ~ rich, data=d %>% filter(dens=="sub" & rand==F)) %>% summary
lm(n.I ~ rich, data=d %>% filter(dens=="sub" & rand==T)) %>% summary

#problems with the models: 
#the standard error is not supposed to be the same for each slope coefficient. This is because variance around the fit is not homogeneous? (also not normal)
plot(resid(m1)~fitted(m1))
plot(resid(m3)~fitted(m3))

#respecify the model so account for the variation?
head(d)
d <- d %>% mutate( rep3 = case_when(rand==F~"det", T~as.character(rep2)) ) #make all the blocks the same for the deterministic trt
ggplot(d, aes(n.I))+
  geom_histogram()+
  facet_wrap(~dens+rand)
m4 <- glmer.nb(n.I ~ -1 + trt*rich + (1|rep2), data=d)
m5 <- MASS::glm.nb(n.I ~ -1 + (trt+rep2)*rich, data=d)
m5.1 <- MASS::glm.nb(n.I ~  (trt)*rich, data=d)
#m6 <- MASS::glm.nb(n.I ~ -1 + trt*rich + rich*rep3, data=filter(d, rand==T))

#model assessment. heteroscedastic!
model=m5.1
summary(model)
hist(resid(model))
qqnorm(resid(model))+qqline(resid(model))
plot(resid(model) ~ fitted(model))
plot(d$n.I~ fitted(model))+abline(0,1)#it's curvey. a non-linear curve probably better.
lstrends(model, "trt", var="rich") %>% plot

#plot predicted
pred <- predict(model, type="response")
d$pred <- pred
ggplot(d, aes(rich, n.I, group=rep2))+
  geom_point()+
  geom_line()+
  geom_point(aes(rich, pred), color='blue')+
  geom_line(aes(rich, pred), color='blue')+
  facet_wrap(~dens+rand)

#proportion infected
ggplot(d, aes(pI))+
  geom_histogram()+
  facet_wrap(~dens+rand)
#make proportion data into raw counts and cbind
head(d)
y <- cbind(d$n.I, d$n.S)
m6 <- glm(y ~ trt*rich, data=d, family = binomial)
m6.1 <- glm(y ~ trt*rich + rep2, data=d, family = binomial)
m7 <- glm(y ~ trt*rich + rich*rep2, data=d, family = binomial)
m8 <- glmer(y ~ trt*rich + (rich|rep2), data=d, family = binomial)
summary(m6.1)
summary(m7)
summary(m8)
#model eval 
model=m6.1
hist(resid(model))
qqnorm(resid(model))+qqline(resid(model))
plot(resid(model) ~ fitted(model))
plot(d$pI~fitted(model))+abline(0,1)#it's curvey. a non-linear curve probably better.
lstrends(model, "trt", var="rich") %>% plot

#plot predicted
pred <- predict(model, type="response")
d$pred <- pred
ggplot(d, aes(rich, pI, group=rep2))+
  geom_point()+
  geom_line()+
  geom_point(aes(rich, pred), color='blue')+
  geom_line(aes(rich, pred), color='blue')+
  facet_wrap(~dens+rand)
