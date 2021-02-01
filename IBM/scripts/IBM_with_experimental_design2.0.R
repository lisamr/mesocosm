#v2.0 = analysis involves probability of infection at tray- and individidual level and rate of growth (at just tray level i think). Also will try to parse additive and non-additive effects.

#IBM of my experimental design. "IBM_functions.R" is the source code for all the functions. Here, I'll put those to use on the same compostions that will be in the greenhouse.

#use community compositions from experimental design (take output from 'species_distributions.R')

rm(list=ls())
library(rethinking)
library(bayesplot)

#load source code and functions----
source('IBM/scripts/IBM_functions.R')
theme_set(theme_bw() + theme(panel.grid = element_blank()))

unzscore <- function(x, x_z){
  mean(x) + x_z*sd(x)
}
zscore <- function(x){
  (x - mean(x)) / sd(x)
}

predictions <- function(newdata, predictions, observed){ #for the trt models
  preds <- newdata %>% cbind( data.frame(
    mean = apply(predictions, 2, mean),
    lower = apply(predictions, 2, HPDI, .95)[1,],
    upper = apply(predictions, 2, HPDI, .95)[2,]
  ) ) %>% 
    rename(richness_std = richness) %>% 
    mutate(richness = unzscore(observed$richness, richness_std))
  trt <- observed %>% distinct(trt, dens, rand)
  preds <- preds %>% left_join(trt)
  return(preds)
}

predictions2 <- function(newdata, predictions, observed, trt=NULL){ #for the disease risk models
  x <- colnames(newdata)
  y <- as.data.frame(sapply(1:length(x), function(i) unzscore(pull(observed[,x[i]]), newdata[,x[i]])) )
  names(y) <- x
  preds <- y %>% cbind( data.frame(
    mean = apply(predictions, 2, mean),
    lower = apply(predictions, 2, HPDI, .95)[1,],
    upper = apply(predictions, 2, HPDI, .95)[2,]
  ) ) 
  return(preds)
}

fdat_list_binomGLM <- function(observed, covs, newdat){
  dat_list <- list(
    N = nrow(observed),
    K = ncol(covs),
    X = covs,
    I = observed$I,
    n = observed$n,
    Nsims = nrow(newdat),
    Xsim = newdat
  )
  return(dat_list)
}

#parameters for simulating epidemics----
#tray dimensions
width <- 9.5*2.54
length <- width 

#DESIGN (keep the names the same. simplifies downstream functions)
pinoc <- .1 #percent inoculated at start of experiemnt
s <- 6 #max number of species
#spp <- c(paste0('sp_', rep(1:s))) #names of species
spp <- c('radish', 'arugula', 'basil', 'green_rom', 'red_rom', 'butter')
tfinal <- 25 #how many time steps
comp <- c(.6, .3, .2, .1, 0, 0) #vector of relative "competencies"

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


### SECONDARY TRANSMISSION ###
#create a matrix of the beta_ij values. Assume that pairwise values will control the amplitude of beta(t). these will come from empirical data, but will just make them up for now.

#create matrix of amplitudes of the beta_ij 
beta_ij_t <- make_beta_ij_t(comp)

### PRIMARY TRANSMISSION ###
#rate of infection from inoculum to plant. varies by species and time
alpha_i_t <- make_alpha_i_t(comp)

#create community----
#read in list of spatial polygons dataframes for each tray
#spdf_list <- readRDS('GH_output/species_distributions/spdf_list.RDS')
#spdf_list <- readRDS('GH_output/species_distributions/spdf_list_tighterdens.RDS')
spdf_list <- readRDS('GH_output/real_experiment/big_experiment_spdf_list_12272020.RDS')

#plot one of them
plot_maps(spdf_list[[41]]) #can take ~15sec. if error, try again.


#run epidemic----

#for every individual, if state==C, then C->C or C->I; if state==S, then S->S, or S->I; if state==I, then always stays I. the fate of the transition determined by a probabilistic binomial process. Will need to loop through every individual every time step. At each time step, record the state of every individual. Output should be a matrix of states, with rows equalling # individuals and cols equalling # time steps. 

#ptm <- proc.time()# Start the clock!
#IBM_list_NN <- lapply(spdf_list, function(x) IBM(x, "NN") )
#howlongIBM_NN <- proc.time() - ptm# 77 seconds

#decide which state to use
#for(i in 1:length(spdf_list)){
#  spdf_list[[i]]$state0 <- spdf_list[[i]]$state0v2
#  spdf_list[[i]]$state0v2 <- NULL
#}

#ptm <- proc.time()# Start the clock!
#IBM_list_Kernel <- lapply(spdf_list, function(x) IBM(x, "Kernel", spatialdecay = .001) )#spatialdecay=.001 in paper, but more likely .0025 in your system.
#howlongIBM_Kernel <- proc.time() - ptm# 55 seconds

#save IBM output
#saveRDS(IBM_list_NN, 'IBM/outputs/IBM_list_NN.RDS')
#saveRDS(IBM_list_Kernel, 'IBM/outputs/IBM_list_Kernel_ninoc.RDS') #10 ind inoculated proportial to composition
#saveRDS(IBM_list_Kernel, 'IBM/outputs/IBM_list_Kernel_ninocmostcomp.RDS')#10 ind inoculated of most competent species

#read IBM output
#IBM_list_NN <- readRDS('IBM/outputs/IBM_list_NN.RDS')
#IBM_list_Kernel <- readRDS('IBM/outputs/IBM_list_Kernel.RDS')
#IBM_list_Kernel <- readRDS('IBM/outputs/IBM_list_Kernel_ninoc.RDS')
IBM_list_Kernel <- readRDS('IBM/outputs/IBM_list_Kernel_ninocmostcomp.RDS')


#visualize single tray----

#### check out a single tray ####
trayIDs <- suppressWarnings(lapply(1:length(spdf_list), function(i) (spdf_list[[i]]@data)[1,4:9]) %>% bind_rows())  
head(trayIDs)

#first plot summary of S and I
#plotS_I(IBM_list_NN[[1]])
#plotS_I(IBM_list_Kernel[[1]])

#now plot spatial map of the spread (animation about a minute)
#plot 21, 64, 107, 150 showing det/sub series
#plot 1, 41, 84, 127 showing det/add series
#plot_spread_map(spdf_list[[21]], IBM_list_Kernel[[21]], animate = F)
#plot_spread_map(spdf_list[[64]], IBM_list_Kernel[[64]], animate = F)
#plot_spread_map(spdf_list[[107]], IBM_list_Kernel[[107]], animate = F)
#plot_spread_map(spdf_list[[150]], IBM_list_Kernel[[150]], animate = F)


#create individual- and tray-level dataframes-----
#IBM_list_Kernel <- IBM_list_NN
#attributes of inds + state at tfinal
tf <- 25
tmplist <- list(NULL)
for(i in 1:length(spdf_list)){
  state_tf <- IBM_list_Kernel[[i]][,tf]
  res <- spdf_list[[i]]@data
  res$state_tf <- state_tf
  tmplist[[i]] <- res
}
df_inds <- suppressWarnings(bind_rows(tmplist) ) %>% 
  mutate(I = as.numeric(state_tf == 'I'),
         richness_std = (richness - mean(richness)) / sd(richness),
         trt = as.integer(as.factor(interaction(rand, dens))))
head(df_inds)

#attributes of trays + infections over time
tmplist <- list(NULL)
for(i in 1:169){
  I <- colSums(IBM_list_Kernel[[i]] == 'I') 
  N <- colSums(!is.na(IBM_list_Kernel[[i]]))
  trayID <- i
  time <- 1:length(N)
  tmplist[[i]] <- cbind(trayID, time, N, I)
}
I_time <- as.data.frame(do.call(rbind, tmplist))
df_tray <- left_join(trayIDs, I_time) %>% 
  mutate(percI = I/N,
         richness_std = (richness - mean(richness)) / sd(richness),
         trt = as.integer(as.factor(interaction(rand, dens))))
head(df_tray)

#visualize spread over time
ggplot(df_tray, aes(time, I/N, group = trayID)) +
  geom_line(alpha = .2) +
  facet_grid(cols = vars(richness),
             rows = vars(rand, dens)) 



#AUC as disease spread------
#AUC ~ gamma(mu, phi)
#log(mu) = a0 + beta*richness
#phi = 

#function to calculate AUC. treats segments as trapezoids
AUC <- function(time, y){
  require(zoo)
  id <- order(time)
  AUC <- sum(diff(time[id])*rollmean(y[id],2))
  return(AUC)
}

df_AUC <- df_tray %>% 
  filter(SD == .5) %>% 
  group_by(trayID, rand, dens, richness, richness_std, trt) %>% 
  summarise(AUC = AUC(time, percI)) %>% 
  ungroup()



#playing with priors
N =10000
tmp <- rnorm(N, 0, 1)
exp(rnorm(N, .5, 1) ) %>% dens

trts <- df_AUC %>% distinct(trt, dens, rand)
#trt dens  rand 
#<int> <fct> <fct>
#1     1 add   det  
#2     2 add   stoch
#3     3 sub   det  
#4     4 sub   stoch
dat_listAUC <- list(
  AUC = df_AUC$AUC,
  richness = df_AUC$richness_std,
  trt = df_AUC$trt
)

mAUC <- ulam(
  alist(
    AUC ~ gamma(mu, phi),
    mu <- exp(a0[trt] + beta[trt]*richness),
    a0[trt] ~ normal(.5, 1),
    beta[trt] ~ normal(0, 1),
    phi ~ exponential(.5)
    ), data=dat_listAUC, chains=1
)
precis(mAUC, depth = 2, prob = .95)
pp <- sim(mAUC)
ppc_dens_overlay(dat_listAUC$AUC, pp[1:50,])

#look at predictions 
newdat <- expand_grid(trt = 1:4, richness = seq(-1.17, 1.429, length.out = 10))
preds <- predictions(newdat, link(mAUC, newdat), df_AUC)

ggplot(df_AUC, aes(richness, AUC)) +
  geom_jitter(height = 0, width = .1, alpha = .5) +
  geom_line(data = preds, aes(richness, mean), color = grey(.5)) +
  geom_ribbon(data = preds, aes(y = mean, ymin = lower, ymax = upper), alpha = .2) +
  facet_wrap(~trt)


#run this model for each treatment seperately----
#might want to do this if I'm borrowing trays across treatments
gamma_stan <- stan_model('post_quarantine_trials/Stan_models/gamma_GLM.stan')

#generate data lists
ND <- expand_grid(richness_std = unique(df_AUC$richness_std))
datlist_GLM <- function(df){
  Xmat <- cbind(df$richness_std)
  datlist <- list(
    N = nrow(df),
    K = ncol(Xmat),
    y = df$AUC,
    X = Xmat,
    Nsim = nrow(ND),
    Xsim = ND
  )
  return(datlist)
}
datlists_gamma <- lapply(split(df_AUC, df_AUC$trt), datlist_GLM)

# iterate over datalists
AUC_trt_mods <- list(NULL)
for(i in 1:length(datlists_gamma)){
  AUC_trt_mods[[i]] <- sampling(gamma_stan, data = datlists_gamma[[i]], chains = 1, iter = 2000)
}

#plot predictions
pred <- list(NULL)
for(i in 1:length(datlists_gamma)){
  pred[[i]] <- as.data.frame(precis(AUC_trt_mods[[i]], pars = "mu_sim", depth = 2, prob = .90)) %>% 
    select(1:4) %>% 
    cbind(ND) %>% mutate(trt = i)
  names(pred[[i]]) <- c('mean', 'sd', 'lower', 'upper', 'richness_std', 'trt')
}
pred <- bind_rows(pred)
ggplot(pred, aes(richness_std, mean)) +
  geom_point(data = df_AUC, aes(richness_std, AUC), alpha = .4) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2) +
  facet_wrap(~trt)


#AUC vs. Perc Infected-----
head(df_tray)


tmp <- df_tray %>% 
  filter(SD == .5) %>% 
  group_by(trayID) %>% 
  mutate(cumAUC = cumsum(percI)) %>% 
  ungroup() %>% 
  mutate_at(vars(cumAUC, percI), list(std = zscore))

ggplot(tmp, aes(percI, cumAUC )) +
  geom_point(aes(color = time)) +
  facet_wrap(~trt)+
  scale_color_viridis_c()


#animate differences in AUC and percI
tmp2 <- tmp %>% 
  pivot_longer(c(percI_std, cumAUC_std), values_to = "y", names_to = "response") %>% 
  select(y, response, time, richness, trt)
tmp2
p1 <- ggplot(tmp2, aes(richness, y, color = response)) +
  geom_point(alpha = .4, position = position_dodge(width = .4)) +
  facet_wrap(~trt) +
  transition_states(time, .5, .1) +
  ggtitle('time step {closest_state}')




#evaluating individual-level percent infected-------

#tray-level infection prevalence looks identical to AUC
df_tray %>% 
  filter(time == tf, SD == .5) %>% 
  ggplot(., aes(richness, percI)) +
  geom_jitter(height = 0, width = .1, alpha = .5) +
  facet_grid(rows = vars(rand), cols = vars(dens))

#model individual-level infection risk controlling for species identity
df_inds2 <- df_inds %>% 
  filter(SD == .5, state0 != "C") %>% 
  mutate(sp = as.integer(as.factor(spID)), 
         trayID = as.integer(as.factor(trayID)))
dat_listind <- list(
  I = df_inds2$I,
  richness = df_inds2$richness_std,
  trt = df_inds2$trt,
  spID = df_inds2$sp,
  trayID = df_inds2$trayID
)

#CAUTION: TAKES ABOUT 30 min TO RUN (1 chain of 1000 iter)
m_indrisk <- ulam(
  alist(
    I ~ bernoulli(p),
    p <- inv_logit(a0[trt] + beta[trt]*richness + as[spID] + aj[trayID]),
    a0[trt] ~  normal(0, 1),
    beta[trt] ~ normal(0, 1),
    as[spID] ~ dcauchy(0, sigmas),
    aj[trayID] ~ normal(0, sigmaj),
    sigmas ~ exponential(1),
    sigmaj ~  exponential(1)
  ), data=dat_listind, chains=1
)
precis(m_indrisk, depth = 2, pars = c('a0', 'beta', 'as', 'sigmas', 'sigmaj'))

#evaluate fit and predictions
post <- extract.samples(m_indrisk)
pp <- sim(m_indrisk)
ppc_dens_overlay(dat_listind$I, pp[1:50,])

newdat <- expand_grid(trt = 1:4, richness = seq(-1.39, 1.42, length.out = 10), spID = 1:6, trayID = 1)
df_inds3 <- df_inds2 %>% 
  group_by(trayID, sp, richness, trt, rand, dens) %>% 
  summarise(percI = sum(I)/n()) %>% 
  rename(spID = sp) %>% 
  ungroup()
preds1 <- predictions(newdat, link(m_indrisk, newdat), df_inds2)

#plot
ggplot(df_inds3, aes(richness, percI)) +
  geom_jitter(height = 0, width = .1, alpha = .6, size = 2.5, aes(color = as.factor(spID))) +
  geom_line(data = preds1, aes(richness, mean, color = as.factor(spID), group = spID )) +
  geom_ribbon(data = preds1, aes(y = mean, ymin = lower, ymax = upper, fill = as.factor(spID)), alpha = .2) +
  scale_color_brewer(palette = 'RdYlBu') +
  scale_fill_brewer(palette = 'RdYlBu') +
  facet_grid(rows = vars(rand), cols = vars(dens)) 






#more complexity....need to figure out the priors for this. 
m_indrisk_complex <- ulam(
  alist(
    I ~ bernoulli(p),
    p <- inv_logit(a0[trt,spID] + beta[trt,spID]*richness + aj[trayID]),    
    a0 ~  normal(0, 1),
    beta[trt] ~ normal(0, 1),
    as[spID] ~ dcauchy(0, sigmas),
    aj[trayID] ~ normal(0, sigmaj),
    sigmas ~ exponential(1),
    sigmaj ~  exponential(1)
  ), data=dat_listind, chains=1
)




#In the meantime, fit the models seperately for sp1 and sp2 specifically. (the "amplifier" and the "amplified"?) intercept and slope need to come from a MVN distribution. 
#I ~ binomial(n, p)
#p[i] = a0[trt[i]] + br[trt[i]]*richness[i];

#prep the data
df_sp1 <- df_inds2 %>% 
  #filter(spID == 'sp_1') %>% 
  filter(spID == 'radish') %>% 
  group_by(trayID, rand, dens, richness, trt) %>% 
  summarize(I = sum(I),
            n = n()) %>% 
  ungroup() %>% 
  mutate(richness_std = zscore(richness))
df_sp2 <- df_inds2 %>% 
  #filter(spID == 'sp_2') %>% 
  filter(spID == 'arugula') %>% 
  group_by(trayID, rand, dens, richness, trt) %>% 
  summarize(I = sum(I),
            n = n()) %>% 
  ungroup() %>% 
  mutate(richness_std = zscore(richness))
sim_sepspp <- expand_grid(trt = 1:4, richness = seq(-2, 2, length.out = 10))

dat_listsp1 <- list(
  N = nrow(df_sp1),
  `T` = max(df_sp1$trt),
  I = as.integer(df_sp1$I),
  n = as.integer(df_sp1$n),
  richness = df_sp1$richness_std,
  trt = df_sp1$trt,
  Nsims = nrow(sim_sepspp),
  trt_sim = sim_sepspp$trt,
  richness_sim = sim_sepspp$richness
)
dat_listsp2 <- list(
  N = nrow(df_sp2),
  `T` = max(df_sp2$trt),
  I = as.integer(df_sp2$I),
  n = as.integer(df_sp2$n),
  richness = df_sp2$richness_std,
  trt = df_sp2$trt,
  Nsims = nrow(sim_sepspp),
  trt_sim = sim_sepspp$trt,
  richness_sim = sim_sepspp$richness
)

#run model
mod_sepspp <- stan_model('post_quarantine_trials/Stan_models/trt_effects_separatespp.stan')
fit_sp1 <- sampling(mod_sepspp, data = dat_listsp1, chains = 1, iter = 2000, cores=4)
fit_sp2 <- sampling(mod_sepspp, data = dat_listsp2, chains = 1, iter = 2000, cores=4)

#check out
fit_sp1
precis(fit_sp1, pars = c("a0", "br", "sigma_pars", 'rho[2,1]'), depth = 3)
precis(fit_sp2, pars = c("a0", "br", "sigma_pars", 'rho[2,1]'), depth = 3)
#traceplot(fit_sp1, pars = c("a0", "br", "sigma_pars", 'rho[2,1]'))
post_sp1 <- extract.samples(fit_sp1)
post_sp2 <- extract.samples(fit_sp2)
ppc_dens_overlay(dat_listsp1$I, post_sp1$y_rep[1:50,]) #beta-binomial is MUCH better :)
ppc_dens_overlay(dat_listsp1$I, post_sp1$y_rep[1:50,])

#predictions
df_spp <- bind_rows(
  df_sp1 %>% mutate(sp = 'sp1'),
  df_sp2 %>% mutate(sp = 'sp2')
  ) %>% 
  mutate(percI = I/n)

preds_sp1 <- predictions(sim_sepspp, post_sp1$p_sim, df_sp1)
preds_sp2 <- predictions(sim_sepspp, post_sp2$p_sim, df_sp2)
pred_spp <- bind_rows(
  preds_sp1 %>% mutate(sp = 'sp1'),
  preds_sp2 %>% mutate(sp = 'sp2')
) 
pal <- RColorBrewer::brewer.pal(6, 'RdYlBu')
ggplot(df_spp, aes(richness, percI, group = sp)) +
  geom_jitter(height = 0, width = .1, alpha = .5, aes( color = sp)) +
  geom_line(data = pred_spp, aes(richness, mean, color = sp)) +
  geom_ribbon(data = pred_spp, aes(y = mean, ymin = lower, ymax = upper, fill = sp), alpha = .2) +
  facet_grid(rows = vars(rand), cols = vars(dens)) +
  coord_cartesian(xlim = c(1,6))+
  scale_color_manual(values = pal[1:2]) +
  scale_fill_manual(values = pal[1:2])






#Disease risk models-------

#what best explains infection risk? Can species composition alone explain it or is there something special (i.e. non-additive effects) about richness? 
df_inds_DR <- df_inds %>% 
  #filter( state0 != "C") %>% #include extra trays (remove SD filter)
  mutate(trt =  as.integer(as.factor(interaction(rand,dens))),
         sp = as.integer(as.factor(spID)), 
         trayID = as.integer(as.factor(trayID)))
#update the values of comp 
names(comp) <- spp 
df_inds_DR <- df_inds_DR %>% 
  mutate(comp = recode(spID, !!!comp))

df_trays <- df_inds_DR %>% 
  group_by(trayID, rand, dens, richness) %>% 
  summarise(CC = sum(comp),
            #nsp_1 = sum(spID == 'sp_1'),
            nsp_1 = sum(spID == 'radish'),
            I = as.integer(sum(I)),
            n = n()) %>% 
  ungroup %>% 
  mutate(percI = I/n) %>% 
  mutate_at(c("richness", "CC", "n", 'nsp_1'), list(std = zscore) )
head(df_trays)    

ggplot(df_trays, aes(richness, percI)) +
  geom_point(aes(color = CC))
ggplot(df_trays, aes(richness, CC)) +
  geom_point(aes(color = CC))

#prep data
covs1 <- df_trays %>% select(richness_std, CC_std) %>% as.matrix()
covs2 <- df_trays %>% select(richness_std, CC_std, n_std) %>% as.matrix()
covs3 <- df_trays %>% select(richness_std, nsp_1_std)%>% as.matrix()
newdat_trays1 <- expand_grid(richness = 0, CC = seq(-2.5, 3, length.out = 50)) %>% as.matrix()
newdat_trays2 <- expand_grid(richness = 0, CC = seq(-2.5, 3, length.out = 50), n = unique(df_trays$n_std)) %>% as.matrix()
newdat_trays3 <- expand_grid(richness = seq(-2.5, 2.5, length.out = 50),nsp_1 = unique(df_trays$nsp_1_std)) %>% as.matrix()

dat_listtrays1 <- fdat_list_binomGLM(df_trays, covs1, newdat_trays1)
dat_listtrays2 <- fdat_list_binomGLM(df_trays, covs2, newdat_trays2)
dat_listtrays3 <- fdat_list_binomGLM(df_trays, covs3, newdat_trays3)


#run model
mod_binom <- stan_model('post_quarantine_trials/Stan_models/betabinomial_GLM.stan')
fit_tray1 <- sampling(mod_binom, data = dat_listtrays1, chains = 1, iter = 1000)
fit_tray2 <- sampling(mod_binom, data = dat_listtrays2, chains = 1, iter = 1000)
fit_tray3 <- sampling(mod_binom, data = dat_listtrays3, chains = 1, iter = 1000)
precis(fit_tray1, pars = c("a0", "beta"), depth = 2, prob = .95); colnames(covs1)
precis(fit_tray2, pars = c("a0", "beta"), depth = 2, prob = .95); colnames(covs2)
precis(fit_tray3, pars = c("a0", "beta"), depth = 2, prob = .95); colnames(covs3)

post_tray1 <- extract.samples(fit_tray1)
post_tray2 <- extract.samples(fit_tray2)
post_tray3 <- extract.samples(fit_tray3)
ppc_dens_overlay(dat_listtrays1$I, post_tray1$y_rep[1:50,])#good
ppc_dens_overlay(dat_listtrays2$I, post_tray2$y_rep[1:50,]) #good
ppc_dens_overlay(dat_listtrays3$I, post_tray3$y_rep[1:50,])#good

lootray1 <- loo(fit_tray1)
lootray2 <- loo(fit_tray2)
lootray3 <- loo(fit_tray3)
loo::loo_compare(lootray1, lootray2, lootray3)


#check out predictions against observed
predstrays2 <- predictions2(newdat_trays2, post_tray2$psim, df_trays) %>% mutate(n = round(n, 0))
ggplot(df_trays, aes(CC, percI)) +
  geom_point(aes(color = n)) +
  geom_line(data = predstrays2, aes(CC, mean, group = n, color = n)) +
  geom_ribbon(data = predstrays2, aes(y=mean, ymin = lower, ymax = upper, group = n, fill = n), alpha = .2) +
  scale_color_viridis_c() +
  scale_fill_viridis_c() #+ facet_wrap(~n)




#Disease risk models, seperate species-------

#what best explains infection risk? Can species composition alone explain it or is there something special (i.e. non-additive effects) about richness? 
df_inds_DR <- df_inds %>% 
  mutate(trt =  as.integer(as.factor(interaction(rand,dens))),
         sp = as.integer(as.factor(spID)), 
         trayID = as.integer(as.factor(trayID))) 
#update the values of comp 
names(comp) <- spp 
df_inds_DR <- df_inds_DR %>% 
  mutate(comp = recode(spID, !!!comp)) %>% 
  group_by(trayID) %>% 
  mutate(CC = sum(comp))
print(df_inds_DR, width = Inf)

#calculate FOI[i] (decays with distance)
f <- function(rho, dist){
  exp(-1/rho * dist^2)
}
x <- 0:20
plot(x, f(20, x), type = 'l')#play with the range
f_FOI <- function(data, tray, rho = 20){
  tmp <- data %>% filter(trayID == tray)
  dmat <- as.matrix(dist(cbind(tmp$x, tmp$y), upper = T, diag = T))
  dmat2 <- exp(-1/rho*(dmat^2))
  diag(dmat2) <- 0
  FOI = colSums(dmat2*tmp$comp) 
  return(FOI)
}
FOI <- sapply(unique(df_inds_DR$trayID), function(x) f_FOI(df_inds_DR, x))
df_inds_DR$FOI <- unlist(FOI)
ggplot(df_inds_DR, aes(CC, FOI)) +geom_point(alpha = .2)


#calculate FOI--NN (decays with distance, but blocked by NN)
#fix this function. need to generate adjacency matrix. dont feel like turning this into spatial object so need to figure out next move. 
f_FOI_NN <- function(data, tray, rho = 20){
  tmp <- data %>% filter(trayID == tray)
  dmat <- as.matrix(dist(cbind(tmp$x, tmp$y), upper = T, diag = T))
  dmatround <- round(dmat, 2)
  NN_dist <- unique(sort(dmatround))[2]

  dmat2 <- exp(-1/rho*(dmatround^2))
  diag(dmat2) <- 0
  FOI = colSums(dmat2*tmp$comp) 
  return(FOI)
}
FOI <- sapply(unique(df_inds_DR$trayID), function(x) f_FOI(df_inds_DR, x))
df_inds_DR$FOI <- unlist(FOI)


df_trays <- df_inds_DR %>% 
  group_by(trayID, rand, dens, richness) %>% 
  summarise(CC = sum(comp),
            #nsp_1 = sum(spID == 'sp_1'),
            #I = as.integer(sum(I[spID == 'sp_1'])),
            nsp_1 = sum(spID == 'radish'),
            I = as.integer(sum(I[spID == 'radish'])),
            n = n()) %>% 
  ungroup %>% 
  mutate(percI = I/nsp_1) %>% 
  filter(!is.na(percI)) %>% 
  mutate_at(c("richness", "CC", "n", 'nsp_1'), list(std = zscore) )
head(df_trays)    

ggplot(df_trays, aes(richness, percI)) +
  geom_point(aes(color = CC))
ggplot(df_trays, aes(richness, CC)) +
  geom_point(aes(color = CC))
ggplot(df_trays, aes(CC, percI)) +
  geom_point(aes(color = CC))

#prep data
covs1 <- df_trays %>% select(richness_std, CC_std) %>% as.matrix()
covs2 <- df_trays %>% select(richness_std, CC_std, n_std) %>% as.matrix()
covs3 <- df_trays %>% select(richness_std, nsp_1_std)%>% as.matrix()
newdat_trays1 <- expand_grid(richness = 0, CC = seq(-2.5, 3, length.out = 50)) %>% as.matrix()
newdat_trays2 <- expand_grid(richness = 0, CC = seq(-2.5, 3, length.out = 50), n = unique(df_trays$n_std)) %>% as.matrix()
newdat_trays3 <- expand_grid(richness = 0, nsp_1 = unique(df_trays$nsp_1_std)) %>% as.matrix()

fdat_list_binomGLM2 <- function(observed, covs, newdat){
  dat_list <- list(
    N = nrow(observed),
    K = ncol(covs),
    X = covs,
    I = observed$I,
    n = observed$nsp_1,
    Nsims = nrow(newdat),
    Xsim = newdat
  )
  return(dat_list)
}
dat_listtrays1 <- fdat_list_binomGLM2(df_trays, covs1, newdat_trays1)
dat_listtrays2 <- fdat_list_binomGLM2(df_trays, covs2, newdat_trays2)
dat_listtrays3 <- fdat_list_binomGLM2(df_trays, covs3, newdat_trays3)


#run model
fit_tray1 <- sampling(mod_binom, data = dat_listtrays1, chains = 1, iter = 1000)
fit_tray2 <- sampling(mod_binom, data = dat_listtrays2, chains = 1, iter = 1000)
fit_tray3 <- sampling(mod_binom, data = dat_listtrays3, chains = 1, iter = 1000)
precis(fit_tray1, pars = c("a0", "beta"), depth = 2, prob = .95);colnames(covs1)
precis(fit_tray2, pars = c("a0", "beta"), depth = 2, prob = .95);colnames(covs2)
precis(fit_tray3, pars = c("a0", "beta"), depth = 2, prob = .95);colnames(covs3)

post_tray1 <- extract.samples(fit_tray1)
post_tray2 <- extract.samples(fit_tray2)
post_tray3 <- extract.samples(fit_tray3)
ppc_dens_overlay(dat_listtrays1$I, post_tray1$y_rep[1:50,])#good
ppc_dens_overlay(dat_listtrays2$I, post_tray2$y_rep[1:50,]) #good
ppc_dens_overlay(dat_listtrays3$I, post_tray3$y_rep[1:50,])#good

lootray1 <- loo(fit_tray1)
lootray2 <- loo(fit_tray2)
lootray3 <- loo(fit_tray3)
loo::loo_compare(lootray1, lootray2, lootray3)


#check out predictions against observed
predstrays2 <- predictions2(newdat_trays1, post_tray1$psim, df_trays) 
ggplot(df_trays, aes(CC, percI)) +
  geom_point() +
  geom_line(data = predstrays2, aes(CC, mean)) +
  geom_ribbon(data = predstrays2, aes(y=mean, ymin = lower, ymax = upper), alpha = .2) +
  scale_color_viridis_c() +
  scale_fill_viridis_c() 


predstrays2 <- predictions2(newdat_trays2, post_tray2$psim, df_trays) 
ggplot(df_trays, aes(CC, percI, group = n)) +
  geom_point(aes(color = n)) +
  geom_line(data = predstrays2, aes(CC, mean, color = n)) +
  geom_ribbon(data = predstrays2, aes(y=mean, ymin = lower, ymax = upper, fill = n), alpha = .2) +
  labs(y = 'percI radish') +
  scale_color_viridis_c() +
  scale_fill_viridis_c() #+ facet_wrap(~n)




#species-specific disease risk-------
df_N <- df_inds_DR %>% 
  filter(sp <= 2) %>% 
  select(trayID, sp, I) 

newdatbern_1 <- expand_grid(as.data.frame(newdat_trays1), sp = 1:max(df_N$sp)) 
newdatbern_2 <- expand_grid(as.data.frame(newdat_trays2), sp = 1:max(df_N$sp)) 
newdatbern_3 <- expand_grid(as.data.frame(newdat_trays3), sp = 1:max(df_N$sp)) 

fdat_list_bernGLMM <- function(covs, newdat){
  covs <- cbind(1, covs[df_N$trayID,])
  list(
    N = nrow(df_N),
    S = max(df_N$sp),
    J = max(df_N$trayID),
    K = ncol(covs),
    X = covs,
    I = df_N$I,
    spID = df_N$sp,
    trayID = df_N$trayID,
    Nsims = nrow(newdat),
    Xsims = cbind(1, as.matrix(newdat)[,1:(ncol(newdat)-1)]),
    spIDsims = newdat$sp
  )
}

dat_listbern1 <- fdat_list_bernGLMM(covs1, newdatbern_1)
dat_listbern2 <- fdat_list_bernGLMM(covs2, newdatbern_2)

mod_bern <- stan_model('post_quarantine_trials/Stan_models/bernoulli_GLMM.stan') #probably needs non-centering
fitbern1 <- sampling(mod_bern, data = dat_listbern1, iter = 1000, chains = 1)
fitbern2 <- sampling(mod_bern, data = dat_listbern2, iter = 1000, chains = 1)
precis(fitbern1, pars = c('B', 'Bbar', 'sdj', 'sigma_pars'), depth = 3)
precis(fitbern2, pars = c('B', 'Bbar', 'sdj', 'sigma_pars'), depth = 3)
beepr::beep()

post_bern1 <- extract.samples(fitbern1)
post_bern2 <- extract.samples(fitbern2)

#check out predictions against observed
tmp <- df_trays[df_N$trayID,] 
tmp2 <- newdatbern_2[,-ncol(newdatbern_2)]
tmp3 <- df_N %>% 
  group_by(trayID, sp) %>% 
  summarise(percI = sum(I) / n()) %>% 
  left_join(df_trays %>% select(-I, -percI)) %>% 
  rename(spID = sp) 
predsbern <- predictions2(tmp2, post_bern2$psim, tmp) %>% mutate(spID = rep(1:3, length.out = nrow(tmp2)))

ggplot(tmp3, aes(CC, percI)) +
  geom_point(aes(color = n)) +
  geom_line(data = predsbern, aes(CC, mean, group = n, color = n)) +
  geom_ribbon(data = predsbern, aes(y=mean, ymin = lower, ymax = upper, group = n, fill = n), alpha = .2) +
  scale_color_viridis_c() +
  scale_fill_viridis_c() + facet_grid(~spID) 









#logistic curves------


df_tray_filt <- filter(df_tray, rand == 'det', dens == 'sub') %>% 
  mutate(richness_std = (richness - mean(richness))/ sd(richness))
head(df_tray_filt)


#define baseline S-curve function
log_curve <- function(K, n0, r, time){
  K*n0 / (n0 + (K - n0)*exp(-r*time))
}
log_curve(1, 0.001, .3, 0:25)

library(rethinking)
library(bayesplot)

#playing with priors
N=10000
dens(rbeta(N, 2, 2))#K
dens(rlnorm(N, log(.1), 1)) #n0
dens(rlnorm(N, log(.25), .6))
tmp <- rnorm(N, 0, .5)
dens(inv_logit(tmp))
tmp <- rnorm(N, -3, .75)
dens(exp(tmp))
dens(rgamma2(N, .05, .1))

#create data list
dat_list <- list(
  I = as.integer(df_tray_filt$I),
  N = as.integer(df_tray_filt$N),
  time = df_tray_filt$time,
  richness = df_tray_filt$richness_std
)
str(dat_list)

#1. Binomial likelihood with logistic function as baseline
mlogistic <- ulam(
  alist(
    I ~ binomial(N, p),
    p <- K*n0 / (n0 + (K - n0)*exp(-r*time)),
    K <- inv_logit(k0 + bk*richness),
    r <- exp(r0 + br*richness),
    k0 ~ normal(0, 1),
    bk ~ normal(0, 1),
    r0 ~ normal(-1.5, .5),
    br ~ normal(0, 1),
    n0 ~ dgamma2(.05, .1)
  ), data=dat_list, chains=1 )
precis(mlogistic, digits = 2)

#1.5. rate is non-monotonic (better overall fit...)
mlogistic1.5 <- ulam(
  alist(
    I ~ binomial(N, p),
    p <- K*n0 / (n0 + (K - n0)*exp(-r*time)),
    K <- inv_logit(k0 + bk*richness),
    r <- exp(r0 + br*richness + br2*richness^2),
    k0 ~ normal(0, 1),
    bk ~ normal(0, 1),
    r0 ~ normal(-1.5, .5),
    br ~ normal(0, 1),
    br2 ~ normal(0, 1),
    n0 ~ dgamma2(.05, .1)
  ), data=dat_list, chains=1 )
precis(mlogistic1.5, digits = 2)


#2. tray-level varying intercept?
dat_list2 <- list(
  I = as.integer(df_tray_filt$I),
  N = as.integer(df_tray_filt$N),
  time = df_tray_filt$time,
  richness = df_tray_filt$richness_std,
  trayID = as.integer(as.factor(df_tray_filt$trayID))
)
str(dat_list2)

mlogistic2 <- ulam(
  alist(
    #model
    I ~ binomial(N, p),
    p <- K*n0 / (n0 + (K - n0)*exp(-r*time)),
    K <- inv_logit(k0 + bk*richness + ak[trayID]),
    r <- exp(r0 + br*richness + + br2*richness^2 + ar[trayID]),
    
    #priors
    #parameter k
    k0 ~ normal(0, 1),
    bk ~ normal(0, 1),
    ak[trayID] ~ normal(0, sigmak),
    sigmak ~ exponential(1),
    #parameter r
    r0 ~ normal(-1.5, .5),
    br ~ normal(0, 1),
    br2 ~ normal(0, 1),
    ar[trayID] ~ normal(0, sigmar),
    sigmar ~ exponential(1),
    #parameter n0
    n0 ~ dgamma2(.05, .1)
  ), 
  data=dat_list2, chains=1 )
precis(mlogistic2, depth = 1)





#Evaluate models--------

#goodness of fit
model = mlogistic1.5
pp <- sim(model) #post. predictions
ppc_dens_overlay(df_tray_filt$I, pp[1:50,])

#plot predictions against observed points
pm <- link(model)#posterior mean
df_tray_filt$pmean <- apply(pm$p, 2, mean)
df_tray_filt$plower <- apply(pm$p, 2, HPDI, .9)[1,]
df_tray_filt$pupper <- apply(pm$p, 2, HPDI, .9)[2,]

ggplot(df_tray_filt, aes(time, I/N, group = trayID)) +
  geom_line(aes(y = pmean), alpha = 1, color=' red3') +
  geom_ribbon(aes(ymin = plower, ymax = pupper), alpha = .3) +
  geom_point(alpha = .1) +
  facet_wrap(~richness)



#check out relationship to k and r and richness
newdat1 <- expand_grid(
  time = 0,
  richness = seq(-1.1711142, 1.43, length.out = 10),
  N = 238)
newdat2 <- expand_grid(
  time = 0,
  richness = seq(-1.1711142, 1.43, length.out = 10),
  N = 238, 
  trayID = 1:40)
newdat <- newdat1
pm2 <- link(model, newdat)
K_r <- data.frame(Kmean = apply(pm2$K, 2, mean), 
                  Klower = apply(pm2$K, 2, HPDI, .9)[1,],
                  Kupper = apply(pm2$K, 2, HPDI, .9)[2,],
                  rmean = apply(pm2$r, 2, mean), 
                  rlower = apply(pm2$r, 2, HPDI, .9)[1,],
                  rupper = apply(pm2$r, 2, HPDI, .9)[2,])
K_r <- cbind(newdat, K_r)
head(K_r)

ggplot(K_r, aes(richness, Kmean, group = trayID)) +
  geom_line() +
  geom_ribbon(aes(ymin = Klower, ymax = Kupper), alpha = .2)
ggplot(K_r, aes(richness, rmean, group = trayID)) +
  geom_line() +
  geom_ribbon(aes(ymin = rlower, ymax = rupper), alpha = .2)


