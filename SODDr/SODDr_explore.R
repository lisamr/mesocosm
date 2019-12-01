#Exploring the SODDr model to see if it can be adapted for mesocosm. at the moment, i've created 2 models--an IBM and ODE model. The IBM is very specific, but slow, hard to adapt, and likely has syntax errors.the ODE model is very short, thus easy to adapt and less error prone, but doesn't quite capture what's going on spatially. The SODDr model from Noam Ross/Rich Cobb might be what I want. 
#https://www.noamross.net/archives/2012-11-16-sod-dynamics-1/

#devtools::install_github("noamross/SODDr", force = T)#make sure to NOT update dependent packages
#library(SODDr)
setwd("/Users/lisarosenthal/Box/mesocosm expt/mesocosm.git")
rm(list=ls())

source('SODDr/SODDr_code_modified.R')
library(dplyr)
library(ggplot2)
library(forcats)
library(gganimate)

theme_set(theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#read in tree params from package
treeparms.df <- read.csv((system.file('paper_tree_parms_eq.csv', package="SODDr")), stringsAsFactors = F)
#3 species, only 1 (tanoak) has multiple age classes. each group has parms for SIR and the transitions, space (resource competition) requirements, kernels for dispersal, WAIF (like B_ij). not sure about the difference between kernel.par1 and kernel.par2. assume there's a kernel function with 2 arguments.

################################################
#RUN SIMULATION WITHOUT DISEASE
#several values are missing on purpose. We can calculate them in the absense of disease.
#calculations follow equation 8 in Cobb et al. 2012.
#calculate missing parameters
treeparms.df <- within(treeparms.df, {
  # Calculate the relative space requirements of tanoak size classes based on initial conditions
  space[1:4] <- 0.25 * (sum(S.init[1:4])/S.init[1:4])#so what are the units of space?
  
  # Set recruitment rates to steady-state levels.  For Redwood and Bay, this is simply the mortality rate divided by the density-dependence coefficient at simulation start.  For tanoak, which has multiple size classes, it's a but more involved
  S.recruit[5] <- S.mortality[5]/(1 - sum(S.init * space))
  I.recruit[5] <- I.mortality[5]/(1 - sum(S.init * space))
  S.recruit[6] <- S.mortality[6]/(1 - sum(S.init * space))
  I.recruit[6] <- I.mortality[6]/(1 - sum(S.init * space))
  
  A2 <- S.transition[1]/(S.transition[2] + S.mortality[2])
  A3 <- (S.transition[2]/(S.transition[3] + S.mortality[3])) * A2
  A4 <- (S.transition[3]/S.mortality[4]) * A3
  S.recruit[1] <- (S.transition[1] + S.mortality[1])/(1 - sum(S.init * space)) - (S.recruit[2] * A2 + S.recruit[3] * A3 + S.recruit[4] * A4)
  I.recruit[1] <- S.recruit[1]
  A2 <- NULL
  A3 <- NULL
  A4 <- NULL
})

#dispersal function is defined by `adjacent.dispersal`, which depends on 3 arguments, 1 for distance and 2 for dispersal to local and adjacent?
print(adjacent.dispersal)
adjacent.dispersal(seq(.1,5,by=.1), 1, 2)#depending on the distance, dispersal is either local (self? "within the cell"), to neighbor, or nothing.

#make a grid of trees
locations <- MakeLattice(nx = 20, ny = 20, dist = 1)
head(locations)
nrow(locations)#400 locations. each location can be thought of as a plot containing multiple individuals. This will need to change for your uses because each location will have only 1 individual.

#create matrix of initial population values of S and I. there should be n_locations (400) x n_species-classes x 2 (S, I).
initial.vec <- as.vector(rbind(treeparms.df$S.init, treeparms.df$I.init)) #alternates S and I
init <- matrix(data=initial.vec, nrow=nrow(locations), ncol=2*nrow(treeparms.df), byrow = T)
#fun for desired time (discrete time, supply time steps)
time.steps <- 1:100

#run the model
pop.df <- SODModel(treeparms.df, locations, time.steps, init)
str(pop.df)
#change the following to factors (I think the package should've done that.)
pop.df <- pop.df %>% 
  mutate_at(vars(Species, AgeClass, Disease), as.factor) %>% 
  rename(Population=value)
head(pop.df)#provides S or I. I would want R also.

#check out the results
#plot total population against time. distinguish between ageclass, species, disease.
#average across locations
pop.df.totals <- pop.df %>% group_by(Time, Species, AgeClass, Disease) %>% 
  summarise(TotPop=mean(Population))
ggplot(pop.df.totals, aes(Time, TotPop)) +
  geom_line(aes(lty=fct_rev(Disease), color=AgeClass, group=interaction(AgeClass, Disease))) +
  facet_grid(~Species) +
  scale_color_viridis_d(begin = .1, end=.9)
#plot shows dynamics are not at steady state after 100 time steps.


#make plot like the one in the publication
#time vs percent of small and big tanoaks in total population
paper.df <- pop.df.totals %>% 
  mutate(Size = ifelse(AgeClass %in% c(1, 2), "Small", "Large")) %>%
  group_by(Time, Species, Size) %>% 
  dplyr::summarise(sizePop = sum(TotPop)) %>% 
  ungroup() %>% 
  mutate(Pct = 100 * sizePop/sum(sizePop))
#plot it
ggplot(subset(paper.df, Species == 1), aes(x = Time, y = Pct, lty = Size)) + 
  geom_line(lwd = 1) + 
  scale_y_continuous(limits = c(0, 0.6),expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) + scale_linetype_manual(values = c(1,5)) 

################################################
#RUN SIMULATION WITH DISEASE
#Now, let’s change the initial conditions to include disease and see what happens. I change the populations of tanoak and bay in one pixel to infectious instead of susceptible and run the model:
  
#make one of the S individuals an I for tanoak and bay laurel
init[190, 2] <- init[190, 1]
init[190, 1] <- 0
init[190, 10] <- init[190, 9]
init[190, 9] <- 0
pop2.df <- SODModel(treeparms.df, locations, time.steps, init)
str(pop2.df)

#manipulate data
paper.df <- pop2.df %>% 
  mutate_at(vars(Species, AgeClass, Disease), as.factor) %>% 
  rename(Population=value) %>% 
  group_by(Time, Species, AgeClass, Disease) %>% 
  summarise(TotPop=mean(Population)) %>% 
  mutate(Size = ifelse(AgeClass %in% c(1, 2), "Small", "Large")) %>%
  group_by(Time, Species, Size) %>% 
  summarise(sizePop = sum(TotPop)) %>% 
  ungroup() %>% 
  mutate(Pct = 100 * sizePop/sum(sizePop))

#plot it
ggplot(subset(paper.df, Species == 1), aes(x = Time, y = Pct, lty = Size)) + 
  geom_line(lwd = 1) + 
  scale_y_continuous(limits = c(0, 0.6),expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_linetype_manual(values = c(1,5)) 

################################################
#change up the composition: mostly tanoak, some redwood, no bay laurel. Still has disease. 

#change composition
treeparms.df$S.init <- c(0.7 * treeparms.df$S.init[1:4]/sum(treeparms.df$S.init[1:4]), 0, 0.19)
initial.vec <- as.vector(rbind(treeparms.df$S.init, treeparms.df$I.init))
init <- matrix(data = initial.vec, nrow = nrow(locations), ncol = 2 * nrow(treeparms.df), byrow = TRUE)

#introduce an infected tanoak individual
init[190, 2] <- init[190, 1]
init[190, 1] <- 0


treeparms.df <- within(treeparms.df, {
  # Calculate the relative space requirements of tanoak size classes based on initial conditions
  space[1:4] <- 0.25 * (sum(S.init[1:4])/S.init[1:4])
  
  # Set recruitment rates to steady-state levels.  For Redwood and Bay, this is simply the mortality rate divided by the density-dependence coefficient at simulation start.  For tanoak, which has multiple size classes, it's a but more involved
  
  S.recruit[5] <- S.mortality[5]/(1 - sum(S.init * space))
  I.recruit[5] <- I.mortality[5]/(1 - sum(S.init * space))
  S.recruit[6] <- S.mortality[6]/(1 - sum(S.init * space))
  I.recruit[6] <- I.mortality[6]/(1 - sum(S.init * space))
  
  A2 <- S.transition[1]/(S.transition[2] + S.mortality[2])
  A3 <- (S.transition[2]/(S.transition[3] + S.mortality[3])) * A2
  A4 <- (S.transition[3]/S.mortality[4]) * A3
  S.recruit[1] <- (S.transition[1] + S.mortality[1])/(1 - sum(S.init * space)) - 
    (S.recruit[2] * A2 + S.recruit[3] * A3 + S.recruit[4] * A4)
  
  I.recruit[1] <- S.recruit[1]
  A2 <- NULL
  A3 <- NULL
  A4 <- NULL
})

pop3.df <- SODModel(treeparms.df, locations, time.steps, init)

#transform data
paper3.df <- pop3.df %>% 
  group_by(Time, Species, AgeClass, Disease) %>% 
  summarise(TotPop=mean(value)) %>% 
  mutate(Size = ifelse(AgeClass %in% c(1, 2), "Small", "Large")) %>%
  group_by(Time, Species, Size) %>% 
  dplyr::summarise(sizePop = sum(TotPop)) %>% 
  ungroup() %>% 
  mutate(Pct = 100 * sizePop/sum(sizePop))

#plot it
ggplot(subset(paper3.df, Species == 1), aes(x = Time, y = Pct, lty = Size)) + 
  geom_line(lwd = 1) + 
  scale_y_continuous(limits = c(0, 0.6),expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_linetype_manual(values = c(1,5)) 

########################################################
#mostly redwood
treeparms.df$S.init <- c(0.08 * treeparms.df$S.init[1:4]/sum(treeparms.df$S.init[1:4]), 0, 0.69)

initial.vec = as.vector(rbind(treeparms.df$S.init, treeparms.df$I.init))
init <- matrix(data = initial.vec, nrow = nrow(locations), ncol = 2 * nrow(treeparms.df),  byrow = TRUE)

#introduce infected tanoak and bay laurel
init[190, 2] <- init[190, 1]
init[190, 1] <- 0
init[190, 10] <- init[190, 9]
init[190, 9] <- 0

#update space requirements and recruitment rates
treeparms.df <- within(treeparms.df, {
  # Calculate the relative space requirements of tanoak size classes based
  # on initial conditions
  space[1:4] <- 0.25 * (sum(S.init[1:4])/S.init[1:4])
  
  # Set recruitment rates to steady-state levels.  For Redwood and Bay, this
  # is simply the mortality rate divided by the density-dependence
  # coefficient at simulation start.  For tanoak, which has multiple size
  # classes, it's a but more involved
  
  S.recruit[5] <- S.mortality[5]/(1 - sum(S.init * space))
  I.recruit[5] <- I.mortality[5]/(1 - sum(S.init * space))
  S.recruit[6] <- S.mortality[6]/(1 - sum(S.init * space))
  I.recruit[6] <- I.mortality[6]/(1 - sum(S.init * space))
  
  A2 <- S.transition[1]/(S.transition[2] + S.mortality[2])
  A3 <- (S.transition[2]/(S.transition[3] + S.mortality[3])) * A2
  A4 <- (S.transition[3]/S.mortality[4]) * A3
  S.recruit[1] <- (S.transition[1] + S.mortality[1])/(1 - sum(S.init * space)) - 
    (S.recruit[2] * A2 + S.recruit[3] * A3 + S.recruit[4] * A4)
  
  
  I.recruit[1] <- S.recruit[1]
  A2 <- NULL
  A3 <- NULL
  A4 <- NULL
})

#run model
pop4.df <- SODModel(treeparms.df, locations, time.steps, init)
#transform data
paper4.df <- pop4.df %>% 
  group_by(Time, Species, AgeClass, Disease) %>% 
  summarise(TotPop=mean(value)) %>% 
  mutate(Size = ifelse(AgeClass %in% c(1, 2), "Small", "Large")) %>%
  group_by(Time, Species, Size) %>% 
  dplyr::summarise(sizePop = sum(TotPop)) %>% 
  ungroup() %>% 
  mutate(Pct = 100 * sizePop/sum(sizePop))
#plot it
ggplot(subset(paper4.df, Species == 1), aes(x = Time, y = Pct, lty = Size)) + 
  geom_line(lwd = 1) + 
  scale_y_continuous(limits = c(0, 0.6),expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_linetype_manual(values = c(1,5)) 

########################################################
#create animation
#TotTan adds up pops of big tanoaks at each location and time
#FracDisease quantifies pops of infected big tanoaks divided by total population
tanoak_anim <- pop2.df %>% 
  filter(Species == 1) %>% 
  select(-Species) %>% 
  group_by(Time, Location) %>% 
  summarise(TotTan = sum(value[which(AgeClass %in% c(3,4))]), FracDisease = sum(value[which(Disease=="I" & AgeClass %in% c(3,4))])/TotTan) %>% 
  left_join(as.data.frame(locations), by=c("Location"='location'))

#animate!
anim <- ggplot(data = tanoak_anim, aes(x = x, y = y, fill = FracDisease, alpha = TotTan)) + 
  geom_tile() + 
  scale_alpha(name = "TanoakPopulation",limits = c(0, 0.4), breaks = seq(0, 0.4, 0.1)) + 
  scale_fill_continuous(name = "Fraction Diseased", limits = c(0, 1), breaks = seq(0, 1, 0.25), low = "#408000", high = "#FF8000") +
  theme(rect = element_blank(), line = element_blank()) +
  # Here comes the gganimate specific bits
  labs(title = 'Time step: {frame_time}') +
  transition_time(Time) +
  ease_aes('linear')
#anim_save("SODDr/animation.gif", anim) 

########################################################
#use the model to examine epidemics on the local scale. i.e. each cell is 1 individual as opposed to composite of multiple individuals.
#also, no tanoak. it has complicated ageclasses that don't apply for mesocosm.

#create initial population matrix with randomly placed individuals
init_ind <- matrix(0, nrow=nrow(locations), ncol=2*nrow(treeparms.df), byrow = T)
cell <- sample(c(9,11), 400, replace = T, prob = c(3,1))#more bays than redwood
for(i in 1:nrow(init_ind)){
  init_ind[i,cell[i]] <- .7
}

#introduce an infected bay laurel individual
init_ind[190,] <- 0
init_ind[190, 10] <- .7


treeparms.df <- within(treeparms.df, {
  # Calculate the relative space requirements of tanoak size classes based on initial conditions
  space[1:4] <- 0.25 * (sum(S.init[1:4])/S.init[1:4])
  
  # Set recruitment rates to steady-state levels.  For Redwood and Bay, this is simply the mortality rate divided by the density-dependence coefficient at simulation start.  For tanoak, which has multiple size classes, it's a but more involved
  
  S.recruit[5] <- S.mortality[5]/(1 - sum(S.init * space))
  I.recruit[5] <- I.mortality[5]/(1 - sum(S.init * space))
  S.recruit[6] <- S.mortality[6]/(1 - sum(S.init * space))
  I.recruit[6] <- I.mortality[6]/(1 - sum(S.init * space))
  
  A2 <- S.transition[1]/(S.transition[2] + S.mortality[2])
  A3 <- (S.transition[2]/(S.transition[3] + S.mortality[3])) * A2
  A4 <- (S.transition[3]/S.mortality[4]) * A3
  S.recruit[1] <- (S.transition[1] + S.mortality[1])/(1 - sum(S.init * space)) - 
    (S.recruit[2] * A2 + S.recruit[3] * A3 + S.recruit[4] * A4)
  
  I.recruit[1] <- S.recruit[1]
  A2 <- NULL
  A3 <- NULL
  A4 <- NULL
})

pop_ind.df <- SODModel(treeparms.df, locations, time.steps, init_ind)

#plot total population against time. distinguish between ageclass, species, disease.
#average across locations
pop_ind.totals <- pop_ind.df %>% group_by(Time, Species, AgeClass, Disease) %>% 
  summarise(TotPop=mean(value))
ggplot(pop_ind.totals, aes(Time, TotPop)) +
  geom_line(aes(lty=fct_rev(Disease), color=AgeClass, group=interaction(AgeClass, Disease))) +
  facet_grid(~Species) +
  scale_color_viridis_c(begin = .1, end=.9)

#create animation
#TotTan adds up pops of big tanoaks at each location and time
#FracDisease quantifies pops of infected big tanoaks divided by total population
pop_ind.df %>% 
  filter(!Species==1, value>0) %>% 
  select(-AgeClass) %>% 
  group_by(Time, Location, Disease) %>% 
  summarise(bay=sum(value[which(Species==2)]),
            redwood=sum(value[which(Species==3)])) %>% 
  head

#visualize spread with animation
animdf <- pop_ind.df %>% 
  filter(!Species==1) %>% 
  select(-AgeClass) %>% 
  group_by(Time, Location, Species) %>% 
  summarise(pop=sum(value), #S+I
            pinf=ifelse(pop==0, 0, value[which(Disease=="I")]/sum(value))) %>% 
  left_join(as.data.frame(locations), by=c("Location"='location')) %>% 
  filter(pop>0, Time<31)
#plot
anim <- ggplot(animdf, aes(x, y)) +
  geom_tile(aes(fill=as_factor(Species))) +
  geom_point(aes(size=pinf), alpha=.5) +
  scale_fill_manual(values = c('#ef8a62', '#67a9cf')) +
  theme(rect = element_blank(), line = element_blank()) +
  # Here comes the gganimate specific bits
  labs(title = 'Time step: {frame_time}') +
  transition_time(Time) +
  ease_aes('linear')
#play
animate(anim, fps=4)
  
