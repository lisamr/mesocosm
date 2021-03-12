#extra trays to be planted for the experiment in Feb 2021


source('IBM/scripts/IBM_functions.R')
library(tidyverse)
theme_set(theme_classic())


#DESIGN PARAMETERS--------------------------------
#distances wanted
Dist <- c(1.7)

#spp used (high, med, low competency species)
spp <- c('radish', 'arugula','basil', 'red_rom', 'green_rom', 'butter')

#design
#6 trays for high/med competency species, 4 for low
design <- rbind(
  expand.grid(species = c('arugula'), rep = 1:3),
  expand.grid(species = c('basil', 'red_rom', 'green_rom'), rep = 1:4),
  cbind(species = 'butter', rep = 1)
) %>% 
  arrange(species) %>% 
  mutate(distance = Dist, 
         trayID = 1:nrow(.)+200) #keeping track of these trays by making them in the 200s 

control <- expand.grid(species = spp, rep = 'Control') %>% 
  arrange(species) %>% 
  mutate(distance = 1.5, 
         trayID = 1:nrow(.)+ 300) #controls are in the 300s

#MAKE TRAYS---------------------------------------
Ninoc <- 12
get_grid <- function(df, x, Ninoc){
  design_x <- df[x,]
  which_spp <- which(spp == design_x$species)
  grid <- with(design_x, sample_community(which_spp, planting_dist = distance, ninoc = Ninoc))
  grid$trayID <- design_x$trayID
  grid$rep <- design_x$rep
  grid$richness <- 1
  return(grid)
}

#simulate trays 
set.seed(2020)
#normal dimensions
extras <- list(NULL)
for(i in 1:nrow(design)){
  extras[[i]] <- get_grid(design, i, Ninoc)
}

controls <- list(NULL)
for(i in 1:nrow(control)){
  controls[[i]] <- get_grid(control, i, 0)
}

#extra trays compiled together
extra_trays <- c(extras, controls)




