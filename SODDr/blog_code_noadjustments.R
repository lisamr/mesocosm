#run code as is in the blog use the package. not getting the plots like in the blog. I think there weren't supposed to be infections the first time around. 
setwd("/Users/lisarosenthal/Box/mesocosm expt/mesocosm.git")

library(SODDr)
library(reshape)
library(plyr)
library(ggplot2)
rm(list=ls())

#read in tree params from package
treeparms.df <- read.csv((system.file('paper_tree_parms_eq.csv', package="SODDr")), stringsAsFactors = F)

#several values are missing on purpose. We can calculate them in the absense of disease.
#calculations follow equation 8 in Cobb et al. 2012.
#calculate missing parameters
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


pop.df.totals <- ddply(pop.df, c("Time", "Disease", "Species", "AgeClass"), summarise, TotPop = mean(value))
head(pop.df.totals)
ggplot(pop.df.totals, aes(x = Time, y = TotPop, fill=Disease)) +
  geom_area(position = 'stack', alpha = 0.6) + 
  facet_grid(~Species)
