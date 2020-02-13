#IBM of my experimental design. "IBM_functions.R" is the source code for all the functions. Here, I'll put those to use on the same compostions that will be in the greenhouse.

#use community compositions from experimental design (take output from 'species_distributions.R')

rm(list=ls())

#load source code----
source('IBM/scripts/IBM_functions.R')

#parameters----
#tray dimensions
width <- 9.5*2.54
length <- width 

#DESIGN (keep the names the same. simplifies downstream functions)
pinoc <- .1 #percent inoculated at start of experiemnt
s <- 6 #max number of species
spp <- c(paste0('sp_', rep(1:s))) #names of species
tfinal <- 25 #how many time steps
comp <- c(1, .5, .3, .2, 0, 0) #vector of relative "competencies"

#DISEASE 
#transmission curves
#values approximated from Otten et al. (2003)
beta_curve <- function(x) .2*exp(-3*(log(x/11))^2)
alpha_curve <- function(x) .4*(1-.3)^x
#distance decay in probability of secondary transmission from Kleczkowski et al. (1997)
dist_decay <- function(x, a=14.9, d=.36, tau=3.73, sigma=.00099){
  C <- a*d*exp(d*tau)
  y <- C*exp(-sigma*x^2)
  #standardize to 0-1, relative to a distance of 0
  y0 <- C*exp(-sigma*0^2)
  y/y0
}

### transmission from C -> S ###
#inoculum decay rate invariant of species.
delta <- 1/5 #1/average number of days inoc stays around

### transmission from S -> I ###
#create a matrix of the beta_ij values. Assume that pairwise values will control the amplitude of beta(t). these will come from empirical data, but will just make them up for now.

#create matrix of amplitudes of the beta_ij 
beta_ij_t <- make_beta_ij_t(comp)

### transmission from C -> I ###
#rate of infection from inoculum to plant. varies by species and time
alpha_i_t <- make_alpha_i_t(comp)

#create community----
#read in list of spatial polygons dataframes for each tray
spdf_list <- readRDS('GH_output/species_distributions/spdf_list.RDS')

#plot one of them
plot_maps(spdf_list[[11]]) #can take ~15sec. if error, try again.

#run epidemic----

#for every individual, if state==C, then C->C, C->S, or C->I; if state==S, then S->S, or S->I; if state==I, then always stays I. Each transition has a probability and the fate of the transition determined by a value drawn from a uniform distribution between 0 and 1. Will need to loop through every individual every time step. At each time step, record the state of every individual. Output should be a matrix of states, with rows equalling # individuals and cols equalling # time steps. 

ptm <- proc.time()# Start the clock!
IBM_list_NN <- lapply(spdf_list, function(x) IBM(x, "NN") )
howlongIBM_NN <- proc.time() - ptm# 67 seconds

ptm <- proc.time()# Start the clock!
IBM_list_Kernel <- lapply(spdf_list, function(x) IBM(x, "Kernel", spatialdecay = .001) )
howlongIBM_Kernel <- proc.time() - ptm# 38 seconds

#save IBM output
saveRDS(IBM_list_NN, 'IBM/outputs/IBM_list_NN.RDS')
saveRDS(IBM_list_Kernel, 'IBM/outputs/IBM_list_Kernel.RDS')

#read IBM output
IBM_list_NN <- readRDS('IBM/outputs/IBM_list_NN.RDS')
IBM_list_Kernel <- readRDS('IBM/outputs/IBM_list_Kernel.RDS')

#visualize single tray----

#### check out a single tray ####
trayIDs <- suppressWarnings(lapply(1:length(spdf_list), function(i) (spdf_list[[i]]@data)[1,4:9]) %>% bind_rows())  
head(trayIDs)

#first plot summary of S and I
plotS_I(IBM_list_NN[[1]])
plotS_I(IBM_list_Kernel[[2]])
plotS_I(IBM_list_Kernel[[1]])

#now plot spatial map of the spread (animation about a minute)
plot_spread_map(spdf_list[[1]], IBM_list_NN[[1]], animate = F)
plot_spread_map(spdf_list[[15]], IBM_list_Kernel[[15]], animate = T)
plot_spread_map(spdf_list[[17]], IBM_list_Kernel[[17]], animate = T)

#anim_save('IBM/plots/spread_map_detsub6.gif') #saves last animation

#how do treatments affect D-D relationship----

#need to associate the state changes to treatments, which are contained in the spatial polygons dataframes. Do it for each community. Then you can see how the treatments affect disease across all trays.

trt_states_list <- lapply(1:length(spdf_list), function(i) bind_treatment_to_states(spdf_list[[i]], IBM_list_Kernel[[i]]))

#turn it into one dataframe
trt_states_df <- bind_rows(trt_states_list)
head(trt_states_df)

#plot changes in I (#inf and %inf) over time
filter(trt_states_df, SD==.5) %>% 
  ggplot(., aes(time, I, group=trayID)) +
  geom_line(alpha=.5, color='dodgerblue4') +
  facet_grid(cols = vars(richness),
             rows = vars(rand, dens))
ggsave('IBM/plots/I_over_time_Kernel.pdf')

filter(trt_states_df, SD==.5) %>% 
  ggplot(., aes(time, percI, group=trayID)) +
  geom_line(alpha=.5, color='dodgerblue4') +
  facet_grid(cols = vars(richness),
             rows = vars(rand, dens))
ggsave('IBM/plots/percI_over_time_Kernel.pdf')

#plot I (#inf and %inf) at final time step
p1 <- trt_states_df %>% filter(time==tfinal, SD==.5) %>% 
  ggplot(., aes(richness, I, group = rep)) +
  geom_line(alpha=.5, color='dodgerblue4') +
  facet_grid(cols = vars(dens),
             rows = vars(rand))

p2 <- trt_states_df %>% filter(time==tfinal, SD==.5) %>% 
  ggplot(., aes(richness, percI, group = rep)) +
  geom_line(alpha=.5, color='dodgerblue4') +
  facet_grid(cols = vars(dens),
             rows = vars(rand))

cowplot::plot_grid(p1, p2, labels = c("A", "B"))
ggsave('IBM/plots/disease_tfinal_Kernel.pdf', width = 10, height = 4.5)

#Community competency----
#community competency can be calculated before the simulation runs
#turn spdf_list into dataframes, bind the rows, then summarize by tray and include treatment, what it's community competency is. can then compare to disaese too.
binded_spdf <- suppressWarnings(bind_rows(lapply(1:length(spdf_list), function(x)spdf_list[[x]]@data))) 
CCsummary <- binded_spdf %>% 
  group_by(trayID, rand, dens, richness) %>% 
  summarise(nind = length(spID), CC = sum(comp), relCC =CC/nind)

#merge Com Comp summary with infection data
CCsummary2 <- trt_states_df %>% 
  filter(time==tfinal) %>% 
  select(trayID, I, percI, SD) %>% 
  left_join(CCsummary, ., by = "trayID") 
head(CCsummary2)

#plot CC vs infections
ggplot(CCsummary2, aes(relCC, percI, color=nind)) +
  geom_point(aes(shape=as.factor(richness))) +
  scale_color_viridis_c() +
  labs(x='avg host competency', y='% infected')
ggplot(CCsummary2, aes(CC, I, color=richness)) +
  geom_point() +
  labs(x='community competency', y='number infected')
ggplot(CCsummary2, aes(CC, percI, color=richness)) +
  geom_point() +
  labs(x='community competency', y='% infected') +
  scale_color_viridis_c()
ggsave('IBM/plots/CC_percI.pdf')

#plot CC vs richness
ggplot(CCsummary2, aes(richness, relCC, color=percI, shape=as.factor(SD))) +
  geom_point() 
ggplot(CCsummary2, aes(richness, CC, color=percI, shape=as.factor(SD))) +
  geom_point() 

#plot density vs richness for both treatments.
CCsummary2 %>% 
  filter(SD==.5) %>% 
  ggplot(., aes(richness, nind, color=dens)) + 
  geom_point() + 
  geom_line() + 
  scale_y_continuous(limits=c(0,400)) 

#statistics----
#load stuff
dens <- function(x, ...) plot(density(x), ...)
inv_logit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))
Scale <- function(x){
  (x-mean(x))/sd(x)
}
library(brms)
library(tidybayes)

#analyze relative effects of richness, host competency, and density on probability of infections
head(CCsummary2)

#scale the predictors too.
CCsummary3 <- CCsummary2
CCsummary3$relCC <- Scale(CCsummary3$relCC)
CCsummary3$richness <- Scale(CCsummary3$richness)
CCsummary3$nind_scaled <- Scale(CCsummary3$nind)
CCsummary3$CC <- Scale(CCsummary3$CC)

#set formula
f1 <- bf(I | trials(nind) ~ relCC + richness + nind, family = binomial)
f1s <- bf(I | trials(nind) ~ relCC + richness + nind_scaled, family = binomial)
f2s <- bf(I | trials(nind) ~ CC + richness, family = binomial)

#get prior
get_prior(f1, CCsummary2)
#figure out prior. want something flat.
rnorm(1000, 0, 1.5) %>% inv_logit %>% dens

#set prior
priors1 <- c(set_prior('normal(0,1.5)', class="Intercept"),
             set_prior('normal(0,.5)', class="b"))

#binomial model
fit1 <- brm(formula = f1, data = CCsummary2, prior = priors1, family = binomial,chains = 4, cores = 4)
fit1s <- brm(formula = f1s, data = CCsummary3, prior = priors1, family = binomial,chains = 4, cores = 4)
fit2s <- brm(formula = f2s, data = CCsummary3, prior = priors1, family = binomial,chains = 4, cores = 4) #what I used in my ppt

#check out model
model=fit2s
pp_check(model, nsamples = 100)#fit2s is best, but still not great
pp_check(model, type = 'intervals')#looks like the model is too confident 

#coef plot
get_variables(model)
coefs <- model %>% 
  gather_draws(b_Intercept, b_CC, b_richness)
#create summary table
coefs_sum <- coefs %>% mean_hdci()
#plot
ggplot(coefs_sum, aes(.variable, .value)) +
  geom_pointinterval(position = position_dodge(width = .2), size=1) +
  geom_hline(yintercept = 0, lty=2) +
  labs(x='variable', y='parameter value') +
  theme(text=element_text(size=16))

#predictions
newd_CC <- expand.grid(CC=seq(-2,2,length.out = 100), nind=100, richness=c(unique(CCsummary3$richness)))
fitted_CC <- add_fitted_draws(newd_CC, model)
pred_CC <- add_predicted_draws(newd_CC, model)

#manually calculate HDPI of the predictions and fit, then plot. will help diagnose and be faster.
fitbounds_CC <- median_hdci(fitted_CC)
predbounds_CC <- median_hdci(pred_CC)
#plot
ggplot(fitbounds_CC, aes(CC, .value/100, group=richness)) +
  geom_line(aes(color=richness))+
  geom_ribbon(data=predbounds_CC, aes(y=.prediction/100, ymin = .lower/100, ymax = .upper/100, fill=richness), alpha=.2) +
  #geom_ribbon(aes(ymin = .lower/100, ymax = .upper/100), alpha=.5, fill='black') +
  #geom_point(data=CCsummary2, aes(Scale(CC), percI, color=richness), alpha=.8) +
  labs(x='scaled community competency', y='pr(infected)')+
  scale_color_viridis_c()+
  scale_fill_viridis_c()
ggsave('IBM/plots/predictions_infections_brms.pdf')








