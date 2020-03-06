#explore how I'm going to set up the competency tests. Not sure if I should randomly or uniformly disperse inoculations. Also, not sure if I should calculate secondary infection probablity with aggregated data or disagregated. 
#simulate both distrubions of inoculation and see how analysis and distributions affects conclusions. 

#LOAD SOURCE CODE----
rm(list=ls())
source('IBM/scripts/IBM_functions.R')

#functions to be added to source code...

test_track_individuals <- function(spatialdataframe, IBM_output){
  spatialdataframe=grid_random
  IBM_output=IBM_random
  
  #get final state
  spatialdataframe$state_tf <- IBM_output[,ncol(IBM_output)]
  
  #summarize things about the tray: density, density of species, avg host competency
  spdf <- as.data.frame(spatialdataframe)
  trayinfo1 <- c('density'=nrow(spdf), 'avgCC'=mean(spdf$comp))
  trayinfo2 <- sapply(1:length(spp), function(i) sum(spdf$spID==spp[i]))
  names(trayinfo2) <- spp
  trayinfo <- c(trayinfo1, trayinfo2)
  trayinfo_mat <- matrix(rep(trayinfo, nrow(spdf)), ncol = length(trayinfo), byrow = T)
  colnames(trayinfo_mat) <- names(trayinfo)
  
  #bind trayinfo to individual info
  output <- cbind(as.data.frame(spatialdataframe), trayinfo_mat) 
  return(output)
}

#DESIGN PARAMETERS----
#distances wanted
Dist <- c(1.7)

#spp used
spp <- c('radish', 'arugula', 'pac_choy', 'romaine', 'basil', 'clover')
comp <- c(.6, .3, .1, .1, 0, 0) #vector of relative "competencies"
names(comp) <- spp

#create tray
#tray dimensions
width <- 9.5*2.54
length <- width 
tray <- make_tray(width, length) #contraining the hexes to 9.5 inches.

#MAKE TRAYS----
#interplanting distance is just 1.7cm

#create spatial grid
grid_random <- sample_community(1, perc_inoc = .1, planting_dist = Dist)
#add in column for competencies
grid_random$comp <- comp[grid_random$spID]
plot_maps(grid_random)



#reassign which individuals are inoculated
grid_uniform <- grid_random
grid_uniform$state0 <- c("S")
grid_uniform$state0[c(31, 36, 41, 101, 106, 111, 171, 176, 181)] <- "C"
plot_maps(grid_uniform)

#SIMULATION----

#model parameters
tfinal <- 25 #how many time steps
beta_curve <- function(x) .2*exp(-3*(log(x/11))^2)
alpha_curve <- function(x) .4*(1-.3)^x
delta <- 1/5 #1/average number of days inoc stays around
beta_ij_t <- make_beta_ij_t(comp) #matrix of amplitudes of the beta_ij 
alpha_i_t <- make_alpha_i_t(comp) #rate of infection from inoculum to plant

#random or uniform
IBM_random <- IBM(grid_random, "Kernel", .0015)
IBM_uniform <- IBM(grid_uniform, "Kernel", .0015)

#plot
plotS_I(IBM_random)
plot_spread_map(grid_random, IBM_random, animate = F)
plot_spread_map(grid_uniform, IBM_uniform, animate = F)

#thoughts: uniform limits how many plants can be inoculated. also, if the pathogen takes off, it very well may go beyond the two rows of susceptibles, making the assumption of "distict neighborhood" incorrect. Go with random, stuff more infecteds in there, and don't hav assumptions of neighborhood. 

#analysis----
library(brms)
library(tidybayes)

#do a logistic regression to predict the probability of a secondary infection. do on invididuals. bernoulli trial. 
ind_trials <- test_track_individuals(grid_random, IBM_random)
ind_trials$infected <- ifelse(ind_trials$state_tf=="I", 1, 0)

#ind_trials$trayID <- as.factor(ind_trials$trayID)

#set formula
#f1 <- bf(infected ~ state0 + spID + density + (1|trayID), family = bernoulli)
f1 <- bf(infected ~ state0 + density, family = bernoulli)
f0 <- bf(infected ~ state0, family = bernoulli)

#set priors (probably need to fiddle)
get_prior(f1, ind_trials)
get_prior(f0, ind_trials)
#priors1 <- c(set_prior('normal(0,1.5)', class="Intercept"),
#             set_prior('normal(0,.5)', class="b"),
#             set_prior('exponential(1)', class='sd'))
priors1 <- c(set_prior('normal(0,1.5)', class="Intercept"),
             set_prior('normal(0,.5)', class="b"))

#run the model
ptm <- proc.time()# Start the clock!
fit1 <- brm(formula = f1, data = ind_trials, prior = priors1, chains = 3, cores = 4)
howlongfit1 <-proc.time() - ptm # 89sec

fit0 <- brm(formula = f0, data = ind_trials, prior = priors1, chains = 3, cores = 4)

#assess fit
model=fit0
pp_check(model, nsamples = 100)#looks pretty spot on. takes a while.

#look at coefficents
#coef plot
get_variables(model)
coefs <- model %>% gather_draws(b_Intercept, b_state0S, b_density)
coefs <- model %>% gather_draws(b_Intercept, b_state0S)

#create summary table
coefs_sum <- coefs %>% mean_hdci()

#plot
ggplot(coefs_sum, aes(.variable, .value)) +
  geom_pointinterval(position = position_dodge(width = .2), size=1) +
  geom_hline(yintercept = 0, lty=2) +
  labs(x='variable', y='parameter value') +
  theme(text=element_text(size=16))

#predictions
newd_CC <- expand.grid(state0=c('C', 'S'))
fitted_CC <- add_fitted_draws(newd_CC, model)

#get a table of coef predictions
fitbounds_CC <- median_hdci(fitted_CC) 

#plot. shows probability of infection with primary and secondary transmission. 
fitted_CC %>% 
  ggplot(., aes(y = state0, x = .value)) +
  geom_halfeyeh(.width = .9, size=.1, fatten_point = .1, relative_scale = 2, point_interval = median_hdi) 


