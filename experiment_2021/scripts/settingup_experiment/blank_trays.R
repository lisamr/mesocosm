
rm(list=ls())
source('IBM/scripts/IBM_functions.R')
library(tidyverse)
theme_set(theme_classic())


#DESIGN PARAMETERS----
#distances wanted
Dist <- c(1.7)

#spp used (high, med, low competency species)
spp <- c('radish', 'arugula','basil', 'green_rom', 'red_rom', 'butter')

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

controls <- expand.grid(species = spp, rep = 'Control') %>% 
  arrange(species) %>% 
  mutate(distance = 1.5, 
         trayID = 1:nrow(.)+ max(design$trayID) )

#MAKE TRAYS----
Ninoc <- 12
get_grid <- function(x){
  design_x <- design[x,]
  which_spp <- which(spp == design_x$species)
  grid <- with(design_x, sample_community(which_spp, planting_dist = distance, ninoc = Ninoc))
  grid$trayID <- design_x$trayID
  grid$rep <- design_x$rep
  return(grid)
}

#simulate trays 
set.seed(2020)
#normal dimensions
grid_list <- list(NULL)
for(i in 1:nrow(design)){
  grid_list[[i]] <- get_grid(i)
}

plot_maps(grid_list[[2]])


#manually assign which get challenged
C12hi <- c(31, 36, 41, 87, 92, 97, 143, 148, 153, 199, 204, 209)
C10hi <- c(33,38,82, 87, 92, 145, 150, 194, 199, 204)
C10lo <- c(12, 15, 33, 36, 39, 54, 59, 70, 72, 75)

C12hitrays <- design %>% 
  filter(species != 'radish') %>% 
  pull(trayID)
C10hitrays <- design %>% 
  filter(species == 'radish', distance == 1.7) %>% 
  pull(trayID)
C10lotrays <- design %>% 
  filter(species == 'radish', distance == 2.8) %>% 
  pull(trayID)
 
f <- function(i, challenged){
  grid_list[[i]]$state0 <- as.character(grid_list[[i]]$state0)
  grid_list[[i]]$state0[challenged] <- 'C'
  grid_list[[i]]$state0 <- as.factor(grid_list[[i]]$state0)
  return(grid_list[[i]])
}

grid_list2 <- c(lapply(C10lotrays, function(x) f(x, C10lo)),
  lapply(C10hitrays, function(x) f(x, C10hi)),
  lapply(C12hitrays, function(x) f(x, C12hi)))



#print maps, 2 per page.
plot_list <- lapply(grid_list2, function(x){
  plot_maps(x) + scale_fill_manual(values = 'grey80')
})
pdf('post_quarantine_trials/figures/sandtrials_7indomes_Jan172021.pdf')
plot_list
dev.off()

#get dataframe
df <- do.call(rbind, lapply(grid_list, function(x) x@data))
