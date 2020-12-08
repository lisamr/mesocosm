#blank trays.

#First inoculation since restarting in October 2020. just getting bearings and getting a coarse idea of plant competency. Planted 3 trays per species at 1.7 cm apart ~50% capacity, 10% inoculated.

#load stuff----
rm(list=ls())
source('IBM/scripts/IBM_functions.R')

#parameters----
#tray dimensions
width <- 9.5*2.54
length <- width 

#DESIGN (keep the names the same. simplifies downstream functions)
pinoc <- 0 #setting at 0, assigned manually.
spp <- c('') #leaving blank.

df <- expand_grid(spp, rep = 1:3)

#create trays----
spdf_list <- list(NULL)
for(i in 1:nrow(df)){
  spdf_list[[i]] <- sample_community(which_spp = i, perc_inoc = pinoc, planting_dist = 1.7)
  spdf_list[[i]]$trayID <- i
  spdf_list[[i]]$rep <- df$rep[i]
}

#Visualize trays----

plot_maps(spdf_list[[1]]) 
