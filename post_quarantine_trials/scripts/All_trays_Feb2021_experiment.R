#compiling all of the trays for the Feb 2021 experiment. 

#1. relabeling species from the Dec 2020 experiment and also changing which ones get inoculated (to 12 each tray). need to preserve locations of species in each tray.
#2. removing some trays that are duplicated.
#3. adding some extra single-species and control trays. 

rm(list = ls())
library(raster)
library(tidyverse)
source('IBM/scripts/IBM_functions.R')

#load existing spatial maps
spdf_list <- readRDS('GH_output/real_experiment/big_experiment_spdf_list_03302020.RDS') 

#define some stuff
ninoc <- 12
spp <- c('radish', 'arugula', 'basil','red_rom', 'green_rom',  'butter')
tmp <- spdf_list[[100]]
#relabel species and reassign challenged individuals
f_ninoc <- function(tmp){
  #first recode the species names 
  spkey <- levels(tmp$spID)
  names(spkey) <- spp
  tmp$spID <- fct_recode(tmp$spID, !!!spkey)
  
  #determine how much of each species gets inoculated
  tab <- table(tmp$spID)
  ninoc_sp <- floor(ninoc/nrow(tmp)*tab)
  x <- ninoc - sum(ninoc_sp) #remainder
  x <- sample(tmp$spID, x, replace = T) %>% table() 
  key <- data.frame(ninoc_sp + x) %>% rename('spID' = 'Var1')
  
  #choose which ones get inoculated
  tmp$state0 <- rep('S', nrow(tmp))
  for(i in 1:nrow(key)){
    ii <- key$spID[i]
    challenged <- sample(which(tmp$spID==ii), key$Freq[key$spID==ii])
    tmp$state0[challenged] <- 'C'
  }
  
  #inoculate 10 individuals of the most competent instead?
  most_comp <- which(tab>0)[1]
  key$Freq2 <- rep(0, nrow(key))
  key$Freq2[most_comp] <- ninoc
  tmp$state0v2 <- rep('S', nrow(tmp))
  for(i in 1:nrow(key)){
    ii <- key$spID[i]
    challenged <- sample(which(tmp$spID==ii), key$Freq2[key$spID==ii])
    tmp$state0v2[challenged] <- 'C'
  }
  return(tmp)
}
set.seed(2020)
spdf_list2 <- lapply(spdf_list, f_ninoc)



#remove some trays---------
#I don't actually plant all of these trays since some of them are replicates. As long as I don't contrast treatments, it's fine. It's mostly to save time and seeds.

#single-species radish trays from stochastic treatments and the sub-det-rich 4 treatment.
omit <- c(33, 11, 17, 107:116)
sapply(spdf_list2[omit], function(x) x$trayID[1])#make sure they're the right trays
spdf_list3 <- spdf_list2[-omit]
length(spdf_list3) #now only 156 trays, from 169


#import extra trays-----------
#I have some extra trays to plant--some control trays and extra single-species trays for estimating competence during the big experiment, rather than seperately.
source('post_quarantine_trials/scripts/extra_trays_Feb2021_experiment.R')
spdf_list4 <- c(spdf_list3, extra_trays)



#export------

#spatial dataframe
saveRDS(spdf_list4, 'post_quarantine_trials/output/big_experiment_spdf_list_02082021.RDS')


#data of the maps
tmp <- lapply(spdf_list4, function(x) x@data)

#binds columns, fills in any missing variables with NA
mapdf <- data.table::rbindlist(tmp, fill = TRUE) %>% 
  arrange(trayID)
head(mapdf); tail(mapdf)
write_csv(mapdf, 'post_quarantine_trials/output/big_experiment_design_02082021.csv')


#Maps for binder
spdf_list4 <- readRDS('post_quarantine_trials/output/big_experiment_spdf_list_02082021.RDS') 

maps_all <- lapply(1:length(spdf_list4), function(x) plot_maps(spdf_list4[[x]], point_shape = 16, point_cex = 2.5))
pdf('post_quarantine_trials/figures/maps_all_for_binder02082021.pdf')
maps_all
dev.off()
