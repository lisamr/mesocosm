#adding 3 more treatments: replicating substitutive design at each of the additive densities

design2 <- read.csv("simulation/outputs/design2.csv")
design2$r[design2$r==2.39] <- 2.4 #for some reason R sucks and hates 2.39. changing it slightly.

design.128 <- design2 %>% filter(d=="add"| d=="sub" & n==128)
design.210 <- design2 %>% filter(d=="add"| d=="sub" & n==200)
design.338 <- design2 %>% filter(d=="add"| d=="sub" & n==338)
design.392 <- design2 %>% filter(d=="add"| d=="sub" & n==392)

#load data. include density and random treatment
Random=F
density="sub"
data <- prepdata2(density, Random, design.392)

S2 <- s2(1, data, B, design.392, t = 20) #can ignore warnings. Has to do with binding results together
animate(S2$HPraster, S2$simulate)

#bind all the treatments together. resnormal is what has been orignally run and the res with numbers are those with additional substitutive treatments at different densities.
resnormal <- fb(B, 10)[[1]] %>% mutate(densN=NA)
res128 <- fb2(B,30, 10, design.128) %>% mutate(densN="sub128")
res210 <- fb2(B,30, 10, design.210) %>% mutate(densN="sub210")
res338 <- fb2(B,30, 10, design.338) %>% mutate(densN="sub338")
res392 <- fb2(B,30, 10, design.392) %>% mutate(densN="sub392")

l <- list(resnormal, res128, res210, res338, res392)
resbind <- bind_rows(l)
head(resbind)

results <- resbind %>% 
  mutate("comp"=Rename(c(1, .3, .2, .1, 0, 0), c(1:6), species)) %>%
  mutate("comp_Abund"=comp*n) %>% 
  mutate(rand=Rename(c("deterministic", "stochastic"), c("FALSE", "TRUE"), rand),
         dens = Rename(c("additive", "substitutive"), c("add", "sub"), dens)) %>% 
  filter(time==max(time)) %>%
  group_by(rep, dens, rand,rich, densN) %>% 
  summarise(N = n[species=="tot"],
            pI=pI[species=="tot"],
            n.I=n.I[species=="tot"],
            comp_Abund=sum(comp_Abund, na.rm = T),
            comp_Abund.N=sum(comp_Abund, na.rm = T)/N) 

ggplot(filter(results, is.na(densN)), aes(rich, pI))+
  geom_point(cex=2)+
  geom_line(aes(group=rep))+
  #geom_smooth(method = "lm", formula = y ~poly(x,3), color="black", lwd=.5)+
  facet_wrap(~rand+dens)+
  labs( y="proportion infected", x="richness") +
  scale_color_manual(values = viridis(2, begin = .25, end= .65, alpha=.8)) +
  theme(strip.background = element_rect(color="black", fill="white"))

#maybe # infected is better than percent infected.
ggplot(filter(results, is.na(densN)), aes(rich, n.I))+
  geom_point(cex=2)+
  geom_line(aes(group=rep))+
  #geom_point(data=results, aes(rich, pI), alpha=.5, color="black")+
  facet_wrap(~rand+dens)+
  labs( y="# infected", x="richness") +
  theme(strip.background = element_rect(color="black", fill="white"))
ggplot(filter(results), aes(comp_Abund, n.I, color=N, shape=rand))+
  geom_point(cex=2)

#relative community comp vs proportion infected, grouped by density.
head(results)
ggplot(filter(results, !N==300), aes(comp_Abund.N, pI, group=N, color=as.factor(N)))+
  geom_point()+
  geom_smooth(se = F, lwd=.5)

#PI generally increases with density, to a point. It increases, then plateues or goes back down at higher densities. how do you show that...
#bin comp_Abund.N into groups
results$comp_Abund.Ncut <-  cut_number(results$comp_Abund.N, 4)

ggplot(filter(results, !N==300), aes(N, pI, group=comp_Abund.Ncut, color=(comp_Abund.Ncut)))+
  geom_point()+
  geom_smooth(se=F, method='lm', formula = y~poly(x,2))
ggplot(filter(results, !N==300), aes(N, pI, color=comp_Abund.N))+
  geom_point()

#try fitting that plot to a parameter surface to see how density and community competency interact to affect incidence. THIS IS IT!!!!
parplot <- function(data, model){
  require(plot3D)
  N.pred <- seq(min(data$N), max(data$N), length.out = 30 )
  comp.pred <- seq(min(data$comp_Abund.N), max(data$comp_Abund.N), length.out = 30 )
  xy <- expand.grid(N = N.pred, comp_Abund.N = comp.pred)
  
  pI.pred <- matrix (nrow = 30, ncol = 30, 
                     data = predict(model, newdata = data.frame(xy), interval = "prediction"))
  
  # predicted z-values, fitted points for droplines to surface
  fitpoints <- predict(model) 
  
  #highlight the points that are additive and deterministic
  data$shape <- ifelse(data$dens=="add" & data$rand==F, 21,20)
  data$cex <- ifelse(data$dens=="add" & data$rand==F, 1.1,1)
  
  scatter3D(z = data$pI, x = data$N, y = data$comp_Abund.N,
            theta = 20, phi = 20, ticktype = "detailed",
            pch=data$shape, cex=data$cex, col = viridis(100, direction=1),
            xlab = "N", ylab = "comp_Abund.N", zlab = "pI", clab = "pI", 
            surf = list(x = N.pred, y = comp.pred, z = pI.pred, 
                        facets = NA, alpha=.7, fit = fitpoints),
            colkey = list(length = 0.8, width = 0.4),            
            main = "Proportion infected")
}


# linear fit
fit <- lm(pI ~ N*comp_Abund.N, data = results)
fit1 <- lm(pI ~ N+comp_Abund.N, data = results)
fit2 <- lm(pI ~ poly(comp_Abund.N, 3) + N, data = results)
fit3 <- lm(pI ~ poly(comp_Abund.N, 3) + poly(N,2), data = results)
fit4 <- lm(pI ~ poly(comp_Abund.N, 3) * poly(N,2), data = results)
fit5 <- lm(pI ~ poly(comp_Abund.N, 3) * poly(N,3), data = results)
AIC(fit, fit1, fit2, fit3, fit4, fit5)

#fit 5 is the best. residuals are normal and variance is equal.
plot(fit5)
hist(resid(fit5))
plot(resid(fit5))
plot(results$pI ~ fitted.values(fit5)) + abline(0,1) #overestimates pI at higher values by a little.
rmse <- function(p, o){
  sqrt(mean((p-o)^2))
}
rmse(fitted.values(fit5), results$pI)
r2 <- lm(results$pI ~ fitted.values(fit5))
summary(r2)

parplot(results, fit5)

#highlight the points in the additive and substitutive designs shown in the first figure.

  require(plot3D)
  N.pred <- seq(min(results$N), max(results$N), length.out = 30 )
  comp.pred <- seq(min(results$comp_Abund.N), max(results$comp_Abund.N), length.out = 30 )
  xy <- expand.grid(N = N.pred, comp_Abund.N = comp.pred)
  
  pI.pred <- matrix (nrow = 30, ncol = 30, 
                     data = predict(fit5, newdata = data.frame(xy), interval = "prediction"))
  
  # predicted z-values, fitted points for droplines to surface
  fitpoints <- predict(fit5) 
  v=viridis(100, direction=1)
  #results2 <- arrange(results, N)
  results$shape <- ifelse(results$dens=="additive" & results$rand=="deterministic", 21,20)
  results$cex <- ifelse(results$dens=="additive" & results$rand=="deterministic", 1.1,1)

  scatter3D(z = results$pI, x = results$N, y=results$comp_Abund.N,
            theta = 20, phi = 20, ticktype = "detailed",
            pch=results$shape, cex=results$cex, col=v,
            xlab = "# individuals", ylab = "Community competency", zlab = "Incidence", clab = "Incidence", 
            surf = list(x = N.pred, y = comp.pred, z = pI.pred, 
                        facets = NA, alpha=.7, fit = fitpoints),
            colkey = list(length = 0.8, width = 0.4))


######################################################################
######################################################################

#I need to be able to make this design empirically tractable. In order to get the parameter surface above, I did 6 times as many treatments. Obviously I cant do that in real life. 
#try filtering out the results and see if you can still come up with the same conclusions.
#keep 1, 2, 4 species at the same additive densities (128, 210, 338 individuals). those densities seem good enough to approximate the parameter surface.

results.f <- results %>% 
  mutate(densN = ifelse(is.na(densN), "orig", densN)) %>% 
  filter(densN!="sub392", rich!=6) 

#compare additive vs sub and assembly orders  
ggplot(filter(results.f, densN=="orig"), aes(rich, pI, color=dens, shape=rand))+
  geom_point()+
  geom_line(aes(group=rep))+
  #geom_smooth(method = "lm", formula = y ~poly(x,3), color="black", lwd=.5)+
  facet_wrap(~rand+dens)+
  labs(color="density", shape="stochastic?", y="proportion infected", x="richness")

#visualize parametric surface
ggplot(filter(results.f, densN!="orig"), aes(N, pI, color=comp_Abund.N))+
  geom_point()

results.f2 <- filter(results.f, densN!="orig")
# linear fit
fit <- lm(pI ~ N*comp_Abund.N, data = results.f2)
fit1 <- lm(pI ~ N+comp_Abund.N, data = results.f2)
fit2 <- lm(pI ~ poly(comp_Abund.N, 3) + N, data = results.f2)
fit3 <- lm(pI ~ poly(comp_Abund.N, 3) + poly(N,2), data = results.f2)
fit4 <- lm(pI ~ poly(comp_Abund.N, 3) * poly(N,2), data = results.f2)
fit5 <- lm(pI ~ poly(comp_Abund.N, 2) * poly(N,2), data = results.f2)
AIC(fit, fit1, fit2, fit3, fit4, fit5)

anova(fit5, fit4) 
#fit4 is the best. can't do 3rd degree polynomial for density because max is D-levels - 1
parplot(results.f2, fit4)

###
#HEY, CAN YOU DO THIS WITH YOUR ORIGINAL DESIGN?????
#not without adding a few more communities. there aren't enough communities with high density and community competency, producing spurious results.
results.orig <- results %>% filter(is.na(densN))

ggplot(results.orig, aes(N, pI, color=comp_Abund.N))+
  geom_point() #looks like the data is missing a few key communities for the surface to be accurately fitted, but maybe you can just augment the design a little.

# linear fit
fit <- lm(pI ~ N*comp_Abund.N, data = results.orig)
fit1 <- lm(pI ~ N+comp_Abund.N, data = results.orig)
fit2 <- lm(pI ~ poly(comp_Abund.N, 3) + N, data = results.orig)
fit3 <- lm(pI ~ poly(comp_Abund.N, 3) + poly(N,2), data = results.orig)
fit4 <- lm(pI ~ poly(comp_Abund.N, 2) * poly(N,2), data = results.orig)
fit5 <- lm(pI ~ poly(comp_Abund.N, 3) * poly(N,2), data = results.orig)
fit6 <- lm(pI ~ poly(comp_Abund.N, 2) * poly(N,3), data = results.orig)
fit7 <- lm(pI ~ poly(comp_Abund.N, 3) * poly(N,3), data = results.orig)
AIC(fit, fit1, fit2, fit3, fit4, fit5, fit6, fit7)

anova(fit7, fit5) 
anova(fit3, fit5) 
#fit5 is the best. lowest AIC. fit7 is what I used when I added those extra treatments, which is not significantly different from fit5.

parplot(results.orig, fit5)
################################################################
################################################################
#How about just make the substitutive design at 392 individuals?
dat <- read.csv("simulation/outputs/moretime_responses.csv")#substitutve has 392 individuals and is run for 30 time steps. See 'moretime_fitnls.R' for code

#treatment results
results.orig2 <- dat %>% 
  mutate("comp"=Rename(c(1, .3, .2, .1, 0, 0), c(1:6), species)) %>%
  mutate("comp_Abund"=comp*n) %>% 
  filter(time==max(time)) %>%
  group_by(rep, dens, rand,rich) %>% 
  summarise(N = n[species=="tot"],
            pI=pI[species=="tot"],
            n.I=n.I[species=="tot"],
            comp_Abund=sum(comp_Abund, na.rm = T),
            comp_Abund.N=sum(comp_Abund, na.rm = T)/N) 

ggplot(results.orig2, aes(rich, pI, color=dens, shape=rand))+
  geom_point()+
  geom_line(aes(group=rep))+
  #geom_smooth(method = "lm", formula = y ~poly(x,3), color="black", lwd=.5)+
  facet_wrap(~rand+dens)+
  labs(color="density", shape="stochastic?", y="proportion infected", x="richness")

#PI ~ density and competency
ggplot(results.orig2, aes(N, pI, color=comp_Abund.N))+
  geom_point()

#parameter surface
# linear fit
fit <- lm(pI ~ N*comp_Abund.N, data = results.orig2)
fit1 <- lm(pI ~ N+comp_Abund.N, data = results.orig2)
fit2 <- lm(pI ~ poly(comp_Abund.N, 3) + N, data = results.orig2)
fit3 <- lm(pI ~ poly(comp_Abund.N, 3) + poly(N,2), data = results.orig2)
fit4 <- lm(pI ~ poly(comp_Abund.N, 2) * poly(N,2), data = results.orig2)
fit5 <- lm(pI ~ poly(comp_Abund.N, 3) * poly(N,2), data = results.orig2)
fit6 <- lm(pI ~ poly(comp_Abund.N, 2) * poly(N,3), data = results.orig2)
fit7 <- lm(pI ~ poly(comp_Abund.N, 3) * poly(N,3), data = results.orig2)
AIC(fit, fit1, fit2, fit3, fit4, fit5, fit6, fit7)

anova(fit7, fit5) 
anova(fit3, fit5) 
#fit5 is the best. lowest AIC. fit7 is what I used when I added those extra treatments, which is not significantly different from fit5.

# predict on x-y grid, for surface
parplot(results.orig2, fit5)

#check assumptions of normal and homo.var resids
hist(resid(fit5))
plot(resid(fit5))
plot(fitted(fit5) ~ results.orig2$pI) + abline(0,1)

#how does diversity affect %infected of species1?
res.sum <- resnormal2[[2]]
head(res.sum)
test <- res.sum %>% filter(rand=="FALSE") %>%  group_by(rep, dens, rich) %>%  summarise(nI.1 = n.I[species==1], pI.1 = pI[species==1], N = n[species=="tot"])
ggplot(test, aes(rich, pI.1, color=dens, group=rep))+
  geom_point() +
  geom_line() +
  facet_wrap(~dens)+
  ylim(0,1)
ggplot(test, aes(rich, nI.1, group=rep))+
  geom_point() +
  geom_line()+
  facet_wrap(~dens)
head(res.sum)

test2 <- res.sum %>%  group_by(rep, dens, rand, rich) %>%  
  summarise(n.1 = ifelse(any(species==1), n[species==1], NA),
            nI.1 = ifelse(any(species==1), n.I[species==1], NA),
            nE.1 = ifelse(any(species==1), n.E[species==1], NA),
            pI.1 = ifelse(any(species==1), pI[species==1], NA),
            N = n[species=="tot"]) %>% 
  right_join(results.orig2)
head(test2)
ggplot(filter(test2, rand==F) , aes(rich, pI.1, group=rep, color=nE.1))+
  geom_point() +
  geom_line()+
  facet_wrap(~dens)
ggplot(filter(test2) , aes(rich, pI.1, group=rep))+
  geom_point() +
  geom_line()+
  facet_wrap(~dens+rand)
ggplot(filter(test2) , aes(comp_Abund.N, pI.1, group=rep, color=N))+
  geom_point() 
ggplot(filter(test2) , aes(N, comp_Abund, group=rep, color=pI.1))+
  geom_point() 
ggplot(filter(test2), aes(n.1, n.I, color=N))+
  geom_point(cex=2)

####
dat=results
#maybe # infected is better than percent infected.
ggplot(filter(dat, is.na(densN)), aes(rich, n.I,  shape=rand))+
  geom_point(cex=2)+
  geom_line(aes(group=rep))+
  facet_wrap(~dens+rand)+
  theme(strip.background = element_rect(color="black", fill="white"))
ggplot(filter(dat, is.na(densN)), aes(rich, comp_Abund,  shape=rand))+
  geom_point(cex=2)+
  geom_line(aes(group=rep))+
  facet_wrap(~dens+rand)+
  theme(strip.background = element_rect(color="black", fill="white"))
ggplot(filter(dat, is.na(densN)), aes(comp_Abund, n.I))+
  geom_point(cex=2)+
  geom_smooth(method = 'lm', formula = y ~ poly(x, 3))+
  labs(x="community competence", y="# infected", color="")

ggplot(filter(dat), aes(N, n.I, color=comp_Abund))+
  geom_point(cex=2)

m0 <- lm(n.I ~1, data = filter(dat, is.na(densN)))
m1 <- lm(n.I ~comp_Abund, data = filter(dat, is.na(densN)))
m2 <- lm(n.I ~poly(comp_Abund,2), data = filter(dat, is.na(densN)))
m3 <- lm(n.I ~poly(comp_Abund,3), data = filter(dat, is.na(densN)))
AIC(m0, m1, m2, m3)
anova(m0, m3)
plot(resid(m3)~fitted(m3))
hist(resid(m3))
summary(m3)
#3d plot


#parameter surface
# linear fit. go with fit7
fit <- lm(n.I ~ N*comp_Abund, data = dat)
fit0 <- lm(n.I ~ N, data = dat)
fit00 <- lm(n.I ~ comp_Abund, data = dat)
fit1 <- lm(n.I ~ N+comp_Abund, data = dat)
AIC(fit, fit0, fit00, fit1)
anova(fit, fit1)

fit2 <- lm(n.I ~ poly(comp_Abund, 3) + N, data = dat)
fit3 <- lm(n.I ~ poly(comp_Abund, 3) + poly(N,2), data = dat)
fit4 <- lm(n.I ~ poly(comp_Abund, 2) * poly(N,2), data = dat)
fit5 <- lm(n.I ~ poly(comp_Abund, 3) * poly(N,2), data = dat)
fit6 <- lm(n.I ~ poly(comp_Abund, 2) * poly(N,3), data = dat)
fit7 <- lm(n.I ~ poly(comp_Abund, 3) * poly(N,3), data = dat)
fit8 <- lm(n.I ~ poly(comp_Abund, 3), data = dat)
AIC(fit, fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8)
summary(fit5)
anova(fit5, fit8) #the addition of density helps

#new data to predict
N.pred <- seq(min(dat$N), max(dat$N), length.out = 30 )
comp.pred <- seq(min(dat$comp_Abund), max(dat$comp_Abund), length.out = 30 )
xy <- expand.grid(N = N.pred, comp_Abund = comp.pred)

nI.pred <- matrix (nrow = 30, ncol = 30, 
                   data = predict(fit7, newdata = data.frame(xy), interval = "prediction"))

# predicted z-values, fitted points for droplines to surface
fitpoints <- predict(fit7) 
v=viridis(100, direction=1)
#dat2 <- arrange(dat, N)
dat$shape <- ifelse(dat$dens=="add" & dat$rand==F, 21,20)
dat$cex <- ifelse(dat$dens=="add" & dat$rand==F, 1.1,1)
#plot
scatter3D(z = dat$n.I, x = dat$N, y=dat$comp_Abund,
          theta = 20, phi = 20, ticktype = "detailed",
          pch=dat$shape, cex=dat$cex, col=v,
          xlab = "# individuals", ylab = "Community competency", zlab = "# infected", clab = "# infected", 
          surf = list(x = N.pred, y = comp.pred, z = nI.pred, 
                      facets = NA, alpha=.7, fit = fitpoints),
          colkey = list(length = 0.8, width = 0.4))
scatter3D(z = dat$n.I, x = dat$N, y=dat$comp_Abund,
          theta = 20, phi = 20, ticktype = "detailed",
          pch=dat$shape, cex=dat$cex, col=v,
          xlab = "# individuals", ylab = "Community competency", zlab = "# infected", clab = "# infected", 
          colkey = list(length = 0.8, width = 0.4))



