#figures for proposal doc and presentation. important that all the figures are using the same simulated data.
#things that I want:
#1. design plot(get that from speciesdistribution.R)
#2. results (indidence and total infected vs. richness for 4 trts, community competence vs total infected, evidence of encounter reduction)
#3. boxplots of slopes for the treatments
#4. scatter3d plot of indidence ~ rel. community comp and density
 
rm(list=ls())
#load functions from 'all_sim_functions2.Rmd'

#load more functions
rmse <- function(p, o){
  sqrt(mean((p-o)^2))
}
parplot <- function(data, model, darkdens="additive"){
  require(plot3D)
  N.pred <- seq(min(data$n), max(data$n), length.out = 30 )
  comp.pred <- seq(min(data$comp_Abund.N), max(data$comp_Abund.N), length.out = 30 )
  xy <- expand.grid(n = N.pred, comp_Abund.N = comp.pred)
  
  pI.pred <- matrix (nrow = 30, ncol = 30, 
                     data = predict(model, newdata = data.frame(xy), interval = "prediction"))
  
  # predicted z-values, fitted points for droplines to surface
  fitpoints <- predict(model) 
  
  #highlight the points that are additive and deterministic
  data$shape <- ifelse(data$dens==darkdens & data$rand=="deterministic", 21,20)
  data$cex <- ifelse(data$dens==darkdens & data$rand=="deterministic", 1.1,1)
  
  scatter3D(z = data$pI, x = data$n, y = data$comp_Abund.N,
            theta = 20, phi = 20, ticktype = "detailed",
            pch=data$shape, cex=data$cex, col = viridis(100, direction=1),
            xlab = "N", ylab = "comp_Abund.N", zlab = "pI", clab = "pI", 
            surf = list(x = N.pred, y = comp.pred, z = pI.pred, 
                        facets = NA, alpha=.7, fit = fitpoints),
            colkey = list(length = 0.8, width = 0.4),            
            main = "Proportion infected")
}
library(dplyr)

#animate a community
design2 <- read.csv("simulation/outputs/design2.csv")
design2$r[design2$r==2.39] <- 2.4 #for some reason R sucks and hates 2.39. changing it slightly.
design.392 <- design2 %>% filter(d=="add"| d=="sub" & n==392)
#choose treatment
Random=F
density="sub"
data <- prepdata2(density, Random,L=50, design.392)
S2 <- s2(4, data, B, design.392, t = 25) #can ignore warnings. Has to do with binding results together
#animate(S2$HPraster, S2$simulate, saveplot = T, Name = "subR6.")
animate(S2$HPraster, S2$simulate)

#analyze data: 30 time steps, 10 reps per treatment. data acquired with code below.
#generate data.
#design2 <- read.csv("simulation/outputs/design2.csv")
#design2$r[design2$r==2.39] <- 2.4 #for some reason R sucks and hates 2.39. changing it slightly.
#design.392 <- design2 %>% filter(d=="add"| d=="sub" & n==392)
#moretime <- fb3(Beta = B, t = 30, n = 10, design = design.392)
#moretimer <- moretime$responses 
#moretimer.sum <- moretime$response.summary
#load data. run for 30 time steps.
dat <- read.csv("simulation/outputs/moretime_responses.csv")
datsum <- read.csv("simulation/outputs/moretime_response_summary.csv")
dat <- dat %>% 
  group_by(rep, dens, rand, rich) %>% 
  mutate(rel.n=n/n[species == "tot"]) %>% 
  mutate("comp"=Rename(c(1, .3, .2, .1, 0, 0), c(1:6), species)) %>%
  mutate("comp_Abund"=comp*n, "comp_relAbund"=comp*rel.n) 



#data at time final
dat30 <- dat %>% 
  filter(time==30, species=="tot") %>%  #look at total infections at time final
  mutate(trt = interaction(dens, rand),
         rep2 = interaction(trt, rep), #rep ID in each trt group
         rep3 = case_when(rand==F~"det", T~as.character(rep2)),
         comp=NULL, comp_Abund=NULL, comp_relAbund=NULL)

#get summary values for community competency and merge
sumstats <- dat %>% 
  filter(time==max(time)) %>% 
  group_by(rep, dens, rand, rich) %>% 
  summarise(comp_Abund=sum(comp_Abund, na.rm = T),
            comp_Abund.N=sum(comp_Abund, na.rm = T)/n[species=="tot"]) 
dat30 <- left_join(dat30, sumstats, by = c("rep", "dens", "rand", "rich"))
dat30 <- ungroup(dat30)

##################################################
#visualize results
##################################################
#1. disease vs richness
plotA=ggplot(dat30, aes(rich, n.I, group=interaction(rand, dens)))+
  geom_point()+
  geom_line(aes(group=rep))+
  facet_wrap(~rand+dens)+
  background_grid(major = "xy", minor = "none") +
  labs(x="richness", y="total infected")
plotB=ggplot(dat30, aes(rich, pI, group=interaction(rand, dens)))+
  geom_point()+
  geom_line(aes(group=rep))+
  facet_wrap(~rand+dens)+
  background_grid(major = "xy", minor = "none") +
  labs(x="richness", y="proportion infected")
plot_grid(plotA, plotB, labels = c("A", "B"))

#2. community competency vs disease
ggplot(dat30, aes(comp_Abund, n.I))+
  geom_point()+
  background_grid(major = "xy", minor = "none") +
  labs(x="community competency (p)", y="total infected")+
  geom_smooth(method='lm', formula = y ~ poly(x, 3))
#statistically test the relationship bewtwen com comp and total infections
m0 <- lm(n.I ~ 1, data=dat30)
m1 <- lm(n.I ~ comp_Abund, data=dat30)
m2 <- lm(n.I ~ poly(comp_Abund, 2), data=dat30)
m3 <- lm(n.I ~ poly(comp_Abund, 3), data=dat30)
AIC(m0, m1, m2, m3)
anova(m3, m2); anova(m3, m0)
#plot(m3) #model looks good
rmse(fitted(m3), dat30$n.I) #15.85785
lm(dat30$n.I~fitted(m3)) %>% summary #Adjusted R-squared:  0.9801 

#3. percent infected of species 1 vs richness. then exposure vs. %infected
head(datsum)
dat2 <- datsum %>% filter(rand==F, species==1) %>% mutate(pE=n.E/n)
ggplot(dat2, aes(rich, pI, group=rep))+
  geom_point()+
  geom_line()+
  facet_grid(~dens)+
  ylim(0,1)
ggplot(dat2, aes(rich, pE, group=rep))+
  geom_point()+
  geom_line()+
  facet_grid(~dens)+
  ylim(0,1)
ggplot(dat2, aes(pE, pI, color=rich))+
  geom_point()+
  facet_grid(~dens)+
  ylim(0,1)

#4. relative community competency vs richness. are they collinear? kinda.
ggplot(dat30, aes(rich, comp_Abund.N))+
  geom_point(aes(color=pI))+
  scale_colour_viridis_c()+
  theme_bw()
ggplot(dat30, aes(comp_Abund.N, pI))+
  geom_point(aes(color=rich))+
  scale_colour_viridis_c()+
  theme_bw()
ggplot(dat30, aes(rich, pI))+
  geom_point(aes(color=comp_Abund.N))+
  scale_colour_viridis_c()+
  theme_bw()
##################################################
#boxplots of slopes for the treatments
##################################################
#get slopes and se for each group
#proportion infected
require(nlme)
dat30 <- dat30 %>% mutate(rand=ifelse(rand==F, "deterministic", "stochastic"), dens=ifelse(dens=="sub", "substitutive", "additive"))
m <- lmList(pI ~ rich | rep2, data=dat30)
cm <- coefficients(m)
slopes.pI <- data.frame(rep2 = row.names(cm), slopes.pI = cm[,2])
dat30 <- left_join(dat30, slopes.pI)
head(dat30)

datslopes <- dat30 %>% filter(!duplicated(rep2))

ggplot(datslopes, aes(dens, slopes.pI, group=interaction(rand, dens), fill=rand))+
  geom_boxplot()+
  labs(x="density", y="slopes: incidence") +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf),fill = "red", alpha = 0.01)+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0),fill = "blue", alpha = 0.01)+
  geom_abline(slope = 0, intercept = 0)+
  geom_boxplot()+
  geom_point(position=position_dodge(.75), alpha=.5)+
  scale_fill_grey(start=.9, end=0.6) +
  theme_bw()

#total infected
m <- lmList(n.I ~ rich | rep2, data=dat30)
cm <- coefficients(m)
slopes.nI <- data.frame(rep2 = row.names(cm), slopes.nI = cm[,2])
dat30 <- left_join(dat30, slopes.nI)
head(dat30)
datslopes <- dat30 %>% filter(!duplicated(rep2))

ggplot(datslopes, aes(dens, slopes.nI, group=interaction(rand, dens), fill=rand))+
  geom_boxplot()+
  labs(x="density", y="slopes: total infections") +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf),fill = "red", alpha = 0.01)+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0),fill = "blue", alpha = 0.01)+
  geom_abline(slope = 0, intercept = 0)+
  geom_boxplot()+
  geom_point(position=position_dodge(.75), alpha=.5)+
  scale_fill_grey(start=.9, end=0.6) +
  theme_bw()

#are the slopes different from zero?
#linear model assumes normal distribution and equal variance of the groups: SE the same. bad.
#one-sample wilcoxon test (nonparametric t-test)
#incidence
wilcox.test(datslopes$slopes.pI[datslopes$trt=="sub.FALSE"], conf.int = T)#-0.127, p=0.002
wilcox.test(datslopes$slopes.pI[datslopes$trt=="sub.TRUE"], conf.int = T)#-0.052, 0.2324
wilcox.test(datslopes$slopes.pI[datslopes$trt=="add.FALSE"], conf.int = T)#-0.083, .002
wilcox.test(datslopes$slopes.pI[datslopes$trt=="add.TRUE"], conf.int = T)#0.023, .084

wilcox.test(datslopes$slopes.nI[datslopes$trt=="sub.FALSE"], conf.int = T)#.002
wilcox.test(datslopes$slopes.nI[datslopes$trt=="sub.TRUE"], conf.int = T)#0.2324
wilcox.test(datslopes$slopes.nI[datslopes$trt=="add.FALSE"], conf.int = T)#.002
wilcox.test(datslopes$slopes.nI[datslopes$trt=="add.TRUE"], conf.int = T)#.002

#bootstrap values to test distribution against zero
#p values = (number of times t_boot > t_obs)/N, N=number of simulated values
boot_test <- function(var, Median=T){
  var_centered <- var-mean(var) #data centered around zero
  if(Median==T){
    t_obs <- median(var) #observed median
    t_boot <- NA
    for(i in 1:5000){
      t_boot[i] <- median(sample(var_centered, size = length(var_centered), replace = T))
  }
  }else{
    t_obs <- mean(var) #observed mean
    t_boot <- NA
    for(i in 1:5000){
      t_boot[i] <- mean(sample(var_centered, size = length(var_centered), replace = T))
    }
  }

  #turn t_obs into neg
  t_obs2 <- ifelse(t_obs>0, t_obs*-1, t_obs)
  #pvalue
  p <- ( sum(t_boot<(t_obs2)) + sum(t_boot>(-t_obs2)) )/5000
  #plot hist
  range <- c(t_obs, t_boot)
  hist(t_boot, xlim=c(min(range), max(range)))
  abline(v=t_obs, col='red')
  return(list(obs.estimate = t_obs, p = p))
}

d <- c("additive", "substitutive")
r <- c("deterministic", "stochastic")
out <- data.frame(expand.grid(d = d, r = r), est=NA, pvalue=NA)
for(i in 1:nrow(out)){
  tmp <- datslopes %>% filter(dens==out[i,1], rand==out[i,2]) %>% pull(slopes.pI)
  tmp2 <- boot_test(tmp, Median = T)
  out$est[i] <- tmp2[[1]]
  out$pvalue[i] <- tmp2[[2]]
}
out

##################################################
#scatter3d plot
##################################################
#try fitting that plot to a parameter surface to see how density and community competency interact to affect incidence. THIS IS IT!!!!
ggplot(dat30, aes(n, pI, color=comp_Abund.N))+
  geom_point()

#parameter surface
# linear fit
fit <- lm(pI ~ n*comp_Abund.N, data = dat30)
fit0 <- lm(pI ~ n*poly(comp_Abund.N,3), data = dat30)
fit0.1 <- lm(pI ~ n+poly(comp_Abund.N,3), data = dat30)
fit1 <- lm(pI ~ n+comp_Abund.N, data = dat30)
fit2 <- lm(pI ~ poly(comp_Abund.N, 3) + n, data = dat30)
fit3 <- lm(pI ~ poly(comp_Abund.N, 3) + poly(n,2), data = dat30)
fit4 <- lm(pI ~ poly(comp_Abund.N, 2) * poly(n,2), data = dat30)
fit5 <- lm(pI ~ poly(comp_Abund.N, 3) * poly(n,2), data = dat30)
fit6 <- lm(pI ~ poly(comp_Abund.N, 2) * poly(n,3), data = dat30)
fit7 <- lm(pI ~ poly(comp_Abund.N, 3) * poly(n,3), data = dat30)
AIC(fit, fit0, fit0.1, fit1, fit2, fit3, fit4, fit5, fit6, fit7)

anova(fit7, fit5) 
anova(fit3, fit5) 
anova(fit0, fit5) #not sigificantly different technically.
anova(fit0, fit7)
#fit5 is the best. lowest AIC. fit7 is what I used when I added those extra treatments, which is not significantly different from fit5.

# predict on x-y grid, for surface
parplot(dat30, fit0, "additive")

#check assumptions of normal and homo.var resids
hist(resid(fit0))
plot(resid(fit0))
plot(fitted(fit0) ~ dat30$pI) + abline(0,1)

#calcuate rmse and r2
rmse(fitted.values(fit0), dat30$pI)
r2 <- lm(dat30$pI ~ fitted.values(fit0))
summary(r2)
