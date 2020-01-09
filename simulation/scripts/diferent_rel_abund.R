#need to see how relative community competency, richness, and percent infected vary. I don't have plots that are both high richness and high rel. C.C., so I'll need to add some plots witht that combo in order to show that CC controls infections, not something inherent about richness.  
 library(forcats)
 library(brms)
 library(purrr)
 library(tidyverse)
 library(modelr)
 library(tidybayes)
 
 rm(list=ls())
 
 #set plot theme
 #theme_set(ggthemes::theme_few())  
 theme_set(theme_bw())
 
 #load functions
 source('simulation/scripts/simulation_code_2.0.R')
 inv_logit <- function(x) exp(x)/(1+exp(x))
 logit <- function(x) log(x/(1-x))

 #run lots of replicates under the current design
 #design2 <- read.csv("simulation/outputs/design2.csv")
 #design2$r[design2$r==2.39] <- 2.4 #for some reason R sucks and hates 2.39. changing it slightly.
 #design.392 <- design2 %>% filter(d=="add"| d=="sub" & n==392)
 #moretime <- fb3(Beta = B, t = 30, n = 10, design = design.392)
 dat <- read.csv("simulation/outputs/moretime_responses.csv")
dat30 <- dat %>% 
   filter(time==30) %>%
   group_by(rep, dens, rand, rich) %>% 
   mutate(rel.n=n/n[species == "tot"]) %>% 
   mutate("comp"=Rename(c(1, .3, .2, .1, 0, 0), c(1:6), species)) %>%
   mutate("comp_Abund"=comp*n, "comp_relAbund"=comp*rel.n) %>% 
   summarise(
     pI= pI[which(species=='tot')] , 
     n=n[which(species=='tot')],
     nI=n.I[which(species=='tot')],
     comp_relAbund = sum(comp_relAbund, na.rm = T))

#relative community competency vs richness. are they collinear? kinda.
 plot1 <- ggplot(dat30, aes(rich, comp_relAbund))+
   geom_point(aes(color=pI))+
   scale_colour_viridis_c()+
   theme_bw()
 plot1
 ggplot(dat30, aes(comp_relAbund, pI))+
   geom_point(aes(color=rich))+
   scale_colour_viridis_c()+
   theme_bw()

 #look at the correlation of com comp and richness a little bit more traditionally. 
 pairs(dat30)
 cor(dat30[,4:6])
 
 #want more trays with high richness and high comp_relAbund
 #alter the relative abundances of the design. the function `prepdata2` is the one that assigns relative abundances.
 #can alter the sd of the species distributions to change how steep the relative abundacnes are. 0=flat, 2=mostly species1

#current design of communities
 des <- tibble(
    d=rep(c("add","sub"), each=4), 
    r=c(1.75, 1.92, 2.39, 3, 1.75, 1.75, 1.75, 1.75), 
    w=floor(25/r), 
    l=floor(50/r), 
    n=w*l, 
    sp=c( 6, 4, 2, 1, 6, 4, 2, 1),
    L=50) %>% 
    group_by(d) %>% 
    arrange(sp, .by_group=T)
  des$r[des$r==2.39] <- 2.4 
prep_D <- prepdata2(dens = 'sub', rand = F, L=unique(des$L), des, SD = 1)
#run simulation. works :)
S2 <- s2(Rich=4, prep_D, B, des, t=20)
animate(S2$HPraster, S2$simulate)

#create additional data with different rel abundances
more_d1 <- fb3(Beta = B, t = 30, L = 50, n=2, design = des, sd=2)
more_d2 <- fb3(Beta = B, t = 30, L = 50, n=2, design = des, sd=1)
more_d1_30 <- more_d1$responses %>% 
  filter(time==30) %>%
  group_by(rep, dens, rand, rich) %>% 
  mutate(rel.n=n/n[species == "tot"]) %>% 
  mutate("comp"=Rename(c(1, .3, .2, .1, 0, 0), c(1:6), species)) %>%
  mutate("comp_Abund"=comp*n, "comp_relAbund"=comp*rel.n) %>% 
  summarise(
    pI= pI[which(species=='tot')] , 
    n=n[which(species=='tot')],
    nI=n.I[which(species=='tot')],
    comp_relAbund = sum(comp_relAbund, na.rm = T))
more_d2_30 <- more_d2$responses %>% 
  filter(time==30) %>%
  group_by(rep, dens, rand, rich) %>% 
  mutate(rel.n=n/n[species == "tot"]) %>% 
  mutate("comp"=Rename(c(1, .3, .2, .1, 0, 0), c(1:6), species)) %>%
  mutate("comp_Abund"=comp*n, "comp_relAbund"=comp*rel.n) %>% 
  summarise(
    pI= pI[which(species=='tot')] , 
    n=n[which(species=='tot')],
    nI=n.I[which(species=='tot')],
    comp_relAbund = sum(comp_relAbund, na.rm = T))

#bind more_bind30 and original dat30
datbind <- bind_rows(more_d1_30, more_d2_30, dat30)

#relative community competency vs richness. are they collinear? nah.
plot1 <- plot1 + theme(legend.position = "none") +
  labs(title='Original Design', y='Relative community competency', x='Richness')
plot2 <- ggplot(datbind, aes(rich, comp_relAbund))+
  geom_point(aes(color=pI))+
  scale_colour_viridis_c()+
  theme_bw() +
  labs(title = 'Additional communities', y='Relative community competency', x='Richness', color='Percent infected')

scatterplot <- cowplot::plot_grid(plot1, plot2, rel_widths = c(.8, 1), labels = c('A', 'B'))
ggsave('simulation/images/different_rel_abund/scatterplot_rich_comcomp.jpg', scatterplot, 'jpeg')
#maybe do a linear regression to see how commuty comp and richness relate to pi?
#2. community competency vs disease
ggplot(datbind, aes(comp_relAbund, pI, color=rich))+
  geom_point()+
  background_grid(major = "xy", minor = "none") +
  labs(x="community competency (p)", y="percent infected")+
  geom_smooth(method='lm', formula = y ~ poly(x, 3))

#statistically test the relationship bewtwen com comp and total infections

#center predictors
center <- function(x){
  x-mean(x)
}
Scale <- function(x){
  (x-mean(x))/sd(x)
}
dens <- function(x, ...) plot(density(x), ...)#density plot function

datbind <- datbind %>% mutate(
  comp_relAbund_s = Scale(comp_relAbund),
  rich_s = Scale(rich))
dat30 <- dat30 %>% mutate(
  comp_relAbund_s = Scale(comp_relAbund),
  rich_s = Scale(rich))

#set formula
f1 <- bf(nI | trials(n) ~ comp_relAbund + rich, family = binomial)
f1s <- bf(nI | trials(n) ~ comp_relAbund_s + rich_s, family = binomial)
f2 <- bf(nI | trials(n) ~ comp_relAbund + rich + dens*rand, family = binomial)

#get prior
get_prior(f1, datbind)
get_prior(f1s, datbind)
get_prior(f2, datbind)
#figure out prior. want something flat.
rnorm(1000, 0, 1.5) %>% inv_logit %>% dens

#set prior
priors1 <- c(set_prior('normal(0,1.5)', class="Intercept"),
  set_prior('normal(0,.5)', class="b"))
priors2 <- c(set_prior('normal(0,1.5)', class="Intercept"))

#binomial model
fit1 <- brm(formula = f1, data = datbind, prior = priors2, family = binomial,chains = 4, cores = 4)
fit1s <- brm(formula = f1s, data = datbind, prior = priors2, family = binomial,chains = 4, cores = 4)
fit2 <- brm(formula = f2, data = datbind, prior = priors2, family = binomial,chains = 4, cores = 4)
fit1_dat30 <- brm(formula = f1, data = dat30, prior = priors2, family = binomial,chains = 4, cores = 4)
fit1s_dat30 <- brm(formula = f1s, data = dat30, prior = priors2, family = binomial,chains = 4, cores = 4)

#save models 
#saveRDS(fit1, 'simulation/outputs/model_dif_rel_abund_fit1.RDS')
#saveRDS(fit1s, 'simulation/outputs/model_dif_rel_abund_fit1s.RDS')
#saveRDS(fit2, 'simulation/outputs/model_dif_rel_abund_fit2.RDS')
#saveRDS(fit1_dat30, 'simulation/outputs/model_dif_rel_abund_fit1_dat30.RDS')
#saveRDS(fit1s_dat30, 'simulation/outputs/model_dif_rel_abund_fit1s_dat30.RDS')

#read models
fit1 <- readRDS('simulation/outputs/model_dif_rel_abund_fit1.RDS')
fit2 <- readRDS('simulation/outputs/model_dif_rel_abund_fit2.RDS')
fit1s <- readRDS('simulation/outputs/model_dif_rel_abund_fit1s.RDS')
fit1_dat30 <- readRDS('simulation/outputs/model_dif_rel_abund_fit1_dat30.RDS')
fit1s_dat30 <- readRDS('simulation/outputs/model_dif_rel_abund_fit1s_dat30.RDS')


model=fit1_dat30
plot(model)
pp_check(model)#not predicting correctly. scaled model way more off. not sure why.
pp_check(model, type = 'intervals')#looks like the model is too confident 

#coef plot
get_variables(fit1s)
coef30 <- fit1s_dat30 %>% 
  gather_draws(b_Intercept, b_comp_relAbund_s, b_rich_s) %>% 
  mutate(Data="status_quo") 
coefbind <- fit1s %>% 
  gather_draws(b_Intercept, b_comp_relAbund_s, b_rich_s) %>% 
  mutate(Data="more_trays") 
coefs <- bind_rows(coef30, coefbind)
#create summary table
coefs_sum <- coefs %>% group_by(.variable, Data) %>% 
  mean_hdci()
coefs_sum$.variable <- factor(coefs_sum$.variable, levels = c("b_Intercept", "b_comp_relAbund_s",  "b_rich_s"))
#plot
coefplot <- ggplot(coefs_sum, aes(.variable, .value, color=Data)) +
  geom_pointinterval(position = position_dodge(width = .2), size=1) +
  geom_hline(yintercept = 0, lty=2) +
  labs(x='variable', y='parameter estimate (predicting Pr(inf))') +
  theme(text=element_text(size=16)) +
  scale_x_discrete(labels=c("b_Intercept"='Intercept', "b_comp_relAbund_s"="com. competency",  "b_rich_s"='richness'))
ggsave('simulation/images/different_rel_abund/rich_vs_comcomp_coefplot.jpg', coefplot, 'jpeg', width = 7, height = 4)


#predictions
model=fit1_dat30
newd_CC <- expand.grid(comp_relAbund=seq(0,1,length.out = 10000), n=100, rich=1)
fitted_CC <- add_fitted_draws(newd_CC, model)
pred_CC <- add_predicted_draws(newd_CC, model)

#manually calculate HDPI of the predictions and fit, then plot. will help diagnose and be faster.
fitbounds_CC <- median_hdci(fitted_CC)
predbounds_CC <- median_hdci(pred_CC)

#plot
predplot <- ggplot(fitbounds_CC, aes(comp_relAbund, .value/100)) +
  geom_line() +
  geom_ribbon(data=predbounds_CC, aes(y=.prediction/100, ymin = .lower/100, ymax = .upper/100), alpha=.2) +
  geom_ribbon(aes(ymin = .lower/100, ymax = .upper/100), alpha=.5, fill='blue4') +
  geom_point(data=datbind, aes(comp_relAbund, pI), alpha=.2) +
  labs(x='relative community competency', y='predicted proporition infected')
ggsave('simulation/images/different_rel_abund/predplot.jpg', predplot, 'jpeg', width = 7, height = 5)

#predictions for richness
newd_R <- expand.grid(comp_relAbund=seq(0,1,by=.1), n=100, rich=1:6)
fitted_R <- add_fitted_draws(newd_R, model)
pred_R <- add_predicted_draws(newd_R, model)

#manually calculate HDPI of the predictions and fit, then plot. will help diagnose and be faster.
fitbounds_R <- fitted_R %>% filter(comp_relAbund==.5) %>% 
  median_hdci() 
predbounds_R <- pred_R %>% 
  ungroup() %>% group_by(rich) %>% 
  median_hdci() 

#plot
ggplot(fitbounds_R, aes(rich, .value)) +
  geom_line() +
  geom_ribbon(data=predbounds_R, aes(y=.prediction, ymin = .prediction.lower, ymax = .prediction.upper), alpha=.2) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha=.5, fill='blue4') +
  geom_point(data=datbind, aes(rich, pI*100), alpha=.2)

 