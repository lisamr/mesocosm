#relabeling species and also changing which ones get inoculated. need to preserve locations of species in each tray.
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



#export------

#spatial dataframe
saveRDS(spdf_list2, 'GH_output/real_experiment/big_experiment_spdf_list_02082021.RDS')
spdf_list2[[2]]@data

#data of the maps
tmp <- lapply(spdf_list2, function(x) x@data)
mapdf <- do.call(rbind, tmp)
head(mapdf)
write_csv(mapdf, 'GH_output/real_experiment/big_experiment_design_02082021.csv')

plot_maps(spdf_list2[[1]], point_shape = 16, point_cex = 2.5)
#Maps for binder
maps_all <- lapply(1:length(spdf_list2), function(x) plot_maps(spdf_list2[[x]], point_shape = 16, point_cex = 2.5))
pdf('GH_plots/maps/maps_all_for_binder02082021.pdf')
maps_all
dev.off()
