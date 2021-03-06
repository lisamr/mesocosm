#explore how I'm going to set up the competency tests. Not sure if I should randomly or uniformly disperse inoculations. Also, not sure if I should calculate secondary infection probablity with aggregated data or disagregated. 
#simulate both distrubions of inoculation and see how analysis and distributions affects conclusions. 
#CONCLUSION: uniform limits how many plants can be inoculated. also, if the pathogen takes off, it very well may go beyond the two rows of susceptibles, making the assumption of "distict neighborhood" incorrect. Go with random, stuff more infecteds in there, and don't hav assumptions of neighborhood

#LOAD SOURCE CODE----
rm(list=ls())
source('IBM/scripts/IBM_functions.R')

#functions to be added to source code...

test_track_individuals <- function(spatialdataframe, IBM_output){
  
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


#IMPORT DESIGN AND SPATIAL DATAFRAME----

design <- read.csv('GH_output/real_experiment/host_competency_design_03092020.csv', row.names = 'X')
spdf_list <- readRDS('GH_output/real_experiment/host_competency_03092020.RDS')

#GENERATE DATAFRAME----

#ID numbers start in bottom left corner and move up L to R
#trays missing last row. fill in with NAs.


add_axes <- function(spatial_object){
  #convert spatial object into dataframe
  spdf <- spatial_object@data
  
  #create axis columns that appear on maps
  #get x and y coordinates
  xs <- unique(sort(round(spdf$x, 4)))
  ys <- spdf$y %>% round(4) %>% sort %>% unique
  
  #translate coordinates into letters or integers
  horiz.ax <- rep(1:ceiling(length(xs)/2), each=2)
  #have to correct for rows that have odd numbers.
  horiz.ax <- if(length(xs)/2 != ceiling(length(xs)/2)){
    horiz.ax[-length(horiz.ax)]
  }else{
    horiz.ax
  }
  horiz.ax <- horiz.ax + rep(c(0,.5), length.out=length(horiz.ax))
  vert.ax <- LETTERS[rev(1:length(ys))]
  names(horiz.ax) <- xs
  names(vert.ax) <- ys
  
  #add axis columns
  spdf <- spdf %>% 
    mutate(x = round(x, 4),
           y = round(y, 4),
           horiz = recode(x, !!!horiz.ax),
           vert = recode(y, !!!vert.ax))
  
  return(spdf)
}

#choose the trays you used. Didn't plant trays 10, 24, 25, 29, 30, 34, 35, 37, 38, 39, 40
design

whichones <- c(1:40)[-c(10, 24, 25, 29, 30, 34, 35, 37, 38, 39, 40)]

#turn grid_list into dataframe to record data
traydf <- suppressWarnings(bind_rows(lapply(whichones, function(i) add_axes(spdf_list[[i]])))) 

traydf <- traydf %>% mutate(day_infected = NA) %>% select(trayID, horiz, vert, everything(), -comp)

write.csv(traydf, 'GH_output/real_experiment/host_competency_03092020_datasheet.csv', row.names = F)

#get state matrix----

#after doing your experiment, enter data and save.
data <- read_csv('GH_data/real_experiment/host_competency_03092020_datasheet.csv')
print(data, width=Inf)

#determine final state of each individual. 
data <- data %>%
  mutate(state_final = case_when(
    is.na(state0) ~ NA_character_, #can't use NA. must designate which NA to use. (NA_real_, NA_character_, NA_integer_, NA_complex_)
    is.na(day_infected) ~ "S",
    day_infected>0 ~ "I"
  )) 

#put into the form of the IBM simulation. 
#generate a matrix of state changes. That can be plotted. Also bind state matrix to the original attribute data. 
tfinal <- max(data$day_infected, na.rm = T) #number of days of monitoring
spp <- unique(data$spID) #names of species

#fill in the state matrix
state_mat <- matrix("S", ncol=tfinal, nrow=nrow(data))
for(i in 2:(tfinal)){
  state_mat[,i] <- state_mat[,i-1]#current time same as last time
  I <- which(data$state_final=='I' & data$day_infected==i)#find infecteds
  state_mat[I,i] <- "I"#update state if infected
}
#replace first column with challenged status and rows as NA
state_mat[,1] <- data$state0
state_mat[is.na(data$state_final),] <- NA

#split state matrix by tray
rows <- data %>%
  mutate(row=row_number()) %>%  
  group_by(trayID) %>% 
  summarise(firstrow=first(row), lastrow=last(row))

state_mat_list <- list(NULL)
for(i in rows$trayID){
  whichtray <- rows %>% filter(trayID == i)
  state_mat_list[[i]] <- state_mat[c(whichtray$firstrow:whichtray$lastrow),]
}

#update the spatial dataframes with any changes in NAs and states (if inoculations moved, species identity changed, or no germination)
for(i in 1:length(state_mat_list)){
  spdf_list[[i]]$state0 <- state_mat_list[[i]][,1]
}

#plot a tray
plotS_I(state_mat_list[[1]])
plot_spread_map(spdf_list[[1]], state_mat_list[[1]], animate = F)


#ANALYSIS----
library(brms)
library(tidybayes)
dens <- function(x, ...) plot(density(x), ...)
inv_logit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))
Scale <- function(x){
  (x-mean(x))/(2*sd(x))
}

#do a logistic regression to predict the probability of a secondary infection. do on invididuals. bernoulli trial. 
ind_trials_list <- lapply(1:length(state_mat_list), function(x) test_track_individuals(spdf_list[[x]], state_mat_list[[x]]))
ind_trials <- bind_rows(ind_trials_list)
ind_trials$infected <- ifelse(ind_trials$state_tf=="I", 1, 0)


#PART 1: SINGLE SPECIES----

#look at single species in each model, just secondary transmission
ind_trials_ss <- filter(ind_trials, spID=="radish")
table(ind_trials_ss$state_tf)

#priors
rnorm(10000, -1, 1.5) %>% density %>% plot
rnorm(10000, -1, 1.5) %>% inv_logit %>%  density %>% plot
f0 <- bf(infected ~ 1 + state0 + (1|trayID), family = bernoulli)
get_prior(f0, ind_trials_ss)
priors0 <- c(set_prior('normal(-1,1.5)', class="Intercept"),
             set_prior('normal(0,.5)', class="b"),
             set_prior('exponential(1)', class='sd'))

#run model
fit0 <- brm(formula = f0, data = ind_trials_ss, prior = priors0, chains = 3, cores = 4)

#look at model 
fit0
plot(fit0)
pp_check(fit0)

#look at random effects
ranef <- fit0 %>% spread_draws(r_trayID[trayID,Intercept]) %>% median_hdci
ggplot(ranef, aes(y=trayID, x=r_trayID)) +
  geom_pointintervalh(fatten_point = .5, size_range = c(.1, .5))

#predictions
newd_CC <- expand.grid(state0=c('C', 'S'), trayID=1)
fitted_CC <- add_fitted_draws(newd_CC, fit0, re_formula = NA)
#get a table of coef predictions
fitbounds_CC <- median_hdci(fitted_CC) 
#plot. shows probability of infection with primary and secondary transmission. 
fitted_CC %>% 
  ggplot(., aes(y = state0, x = .value)) +
  geom_halfeyeh(.width = .9, size=.1, fatten_point = .1, relative_scale = .5, point_interval = median_hdi)



#PART 2: ALL SPECIES----

#look at all species in one model

#set formula
f1 <- bf(infected ~ spID + state0 + (1|trayID), family = bernoulli) #might want density in there for the real data. screwing up the model right now though. 

#set priors (probably need to fiddle)
get_prior(f1, ind_trials)
rnorm(10000, 0, .5) %>% inv_logit %>% density %>% plot(xlim=c(0,1))
priors1 <- c(set_prior('normal(-1,1.5)', class="Intercept"),
             set_prior('normal(0, .5)', class="b"),
             set_prior('exponential(.5)', class='sd'))

#run the model
ptm <- proc.time()# Start the clock!
fit1 <- brm(formula = f1, data = ind_trials, prior = priors1, chains = 3, cores = 4)
howlongfit1 <-proc.time() - ptm # 89sec
beepr::beep(8)

#saveRDS(fit1, 'IBM/outputs/fit1_host_competency.RDS')

#check out model
fit1
plot(fit1)

#look at coefficents
#coef plot
get_variables(fit1)
coefs <- fit1 %>% spread_draws(
  b_Intercept, b_spIDbasil, b_spIDclover, b_spIDpac_choy, b_spIDradish, b_spIDromaine, b_state0S) %>% 
  mutate(arugula = b_Intercept, 
         basil = b_Intercept + b_spIDbasil,
         clover = b_Intercept + b_spIDclover,
         pac_choy = b_Intercept + b_spIDpac_choy,
         radish = b_Intercept + b_spIDradish,
         romaine = b_Intercept + b_spIDromaine) %>% 
  gather_variables()

#create summary table
coefs_sum <- coefs %>% mean_hdci() %>% 
  filter(.variable %in% c('arugula', 'basil', 'clover', 'pac_choy', 'radish', 'romaine', 'b_state0S'))

#plot
ggplot(coefs_sum, aes(.variable, .value)) +
  geom_pointinterval(position = position_dodge(width = .2), size=1) +
  geom_hline(yintercept = 0, lty=2) +
  labs(x='variable', y='parameter value') +
  theme(text=element_text(size=16))

#check out random effects. looks like the random effects are picking up the variation from the species. You shouldn't see patterns. adjust the priors doesn't help. 
ranef <- fit1 %>% spread_draws(r_trayID[trayID,Intercept]) %>% median_hdci
#plot
ggplot(ranef, aes(y=trayID, x=r_trayID)) +
  geom_pointintervalh(fatten_point = .5, size_range = c(.1, .5))

#predictions
newd_CC <- expand.grid(state0=c('C', 'S'), spID = c('arugula', 'basil', 'clover', 'pac_choy', 'radish', 'romaine'), trayID=1)
fitted_CC <- add_fitted_draws(newd_CC, fit1, re_formula = NA)

#get a table of coef predictions
fitbounds_CC <- median_hdci(fitted_CC) 

#plot. shows probability of infection with primary and secondary transmission. 
fitted_CC %>% 
  ggplot(., aes(y = spID, x = .value)) +
  geom_halfeyeh(.width = .9, size=.1, fatten_point = .1, relative_scale = 2, point_interval = median_hdi) +
  facet_grid(~ state0)


#PART 3: ALL SPECIES, RANDOM TRAY INTERCEPT VARYING BY SPECIES----

#look at all species in one model

#set formula
f2 <- bf(infected ~ state0 + (1|trayID) + (1|spID), family = bernoulli) #might want density in there for the real data. screwing up the model right now though. 

#set priors (probably need to fiddle)
get_prior(f2, ind_trials)
rnorm(10000, 0, .5) %>% inv_logit %>% density %>% plot(xlim=c(0,1))
priors1 <- c(set_prior('normal(-1,1.5)', class="Intercept"),
             set_prior('normal(0, .5)', class="b"),
             set_prior('exponential(.5)', class='sd'))



#run the model
ptm <- proc.time()# Start the clock!
fit2 <- brm(formula = f2, data = ind_trials, prior = priors1, chains = 3, cores = 4)
howlongfit1 <-proc.time() - ptm # 89sec
beepr::beep(8)

#saveRDS(fit2, 'IBM/outputs/fit2_host_competency.RDS')

#check out model
fit2
plot(fit2)
pp_check(fit2, nsamples = 100)#looks pretty spot on. takes a while.

#check out random effects
#species
ranef_sp <- fit2 %>%
  spread_draws(r_spID[species,Intercept]) %>%
  median_hdci()
ggplot(ranef_sp, aes(y = fct_rev(species), x = r_spID)) +
  geom_pointintervalh(size_range = c(.1,.5))
#tray
ranef_tray <- fit2 %>%
  spread_draws(r_trayID[tray,Intercept]) %>%
  median_hdci()
ggplot(ranef_tray, aes(y = tray, x = r_trayID)) +
  geom_pointintervalh(size_range = c(.1,.5)) #looks better

#predictions
newd_CC <- expand.grid(state0=c('C', 'S'), spID = c('arugula', 'basil', 'clover', 'pac_choy', 'radish', 'romaine'), trayID=1)
fitted_fit2 <- add_fitted_draws(newd_CC, fit2, re_formula = ~ (1|spID), scale = 'response')

#get a table of coef predictions
fitbounds_fit2 <- median_hdci(fitted_fit2) 

#plot. shows probability of infection with primary and secondary transmission. visually seems to agree with model run with single species.
fitted_fit2 %>% 
  ggplot(., aes(y = spID, x = .value)) +
  geom_halfeyeh(.width = .9, size=.1, fatten_point = .1, relative_scale = 4, point_interval = median_hdi) +
  facet_grid(rows = vars(state0))
#ggsave('IBM/plots/host_competency_estimates_fit2.pdf')


