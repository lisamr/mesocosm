#Exploring the SODDr model to see if it can be adapted for mesocosm. at the moment, i've created 2 models--an IBM and ODE model. The IBM is very specific, but slow, hard to adapt, and likely has syntax errors.the ODE model is very short, thus easy to adapt and less error prone, but doesn't quite capture what's going on spatially. The SODDr model from Noam Ross/Rich Cobb might be what I want. 
#https://www.noamross.net/archives/2012-11-16-sod-dynamics-1/

devtools::install_github("noamross/SODDr", force = T)#make sure to NOT update dependent packages
library(SODDr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(forcats)

#read in tree params from package
treeparms.df <- read.csv((system.file('paper_tree_parms_eq.csv', package="SODDr")), stringsAsFactors = F)
#3 species, only 1 (tanoak) has multiple age classes. each group has parms for SIR and the transitions, space (resource competition) requirements, kernels for dispersal, WAIF (like B_ij). not sure about the difference between kernel.par1 and kernel.par2. assume there's a kernel function with 2 arguments.

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
init <- matrix(data=initial.vec, nrow=nrow(locations), ncol=2*nrow(treeparms.df))
#fun for desired time (discrete time, supply time steps)
time.steps <- 1:100

#run the model
pop.df <- SODModel(treeparms.df, locations, time.steps, init)
str(pop.df)
#change the following to factors (I think the package should've done that.)
pop.df <- pop.df %>% mutate_at(vars(Species, AgeClass, Disease), as.factor)
head(pop.df)#provides S or I. I would want R also.

#check out the results
#plot total population against time. distinguish between ageclass, species, disease.
#average across locations
pop.df.totals <- pop.df %>% dplyr::group_by(Time, Species, AgeClass, Disease) %>% 
  dplyr::summarise(TotPop=mean(Population))
ggplot(pop.df.totals, aes(Time, TotPop)) +
  geom_line(aes(lty=fct_rev(Disease), color=AgeClass, group=interaction(AgeClass, Disease))) +
  facet_grid(~Species) +
  scale_color_viridis_d(begin = .1, end=.9)
#plot shows dynamics are not at steady state after 100 time steps.

#check out total populations without distinguishing age
pop.df.totals2 <- pop.df %>% group_by(Time, Species, Disease) %>% 
  summarise(TotPop=mean(Population))
ggplot(pop.df.totals2, aes(Time, TotPop)) +
  geom_line(aes(color=Disease, group=Disease)) +
  facet_grid(~Species) 

#make plot like the one in the publication
#time vs percent of small and big tanoaks in total population

