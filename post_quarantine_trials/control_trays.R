#printing trays for high density controls (6 trays)

rm(list=ls())
source(here::here('IBM/scripts/IBM_functions.R'))
library(cowplot)
library(lme4)
library(tidyverse)
theme_set(theme_classic())


#DESIGN PARAMETERS----
Dist <- c(1.5)
spp <- c('radish', 'arugula', 'basil', 'green romaine', 'red romaine', 'butter lettuce')

#design
#6 trays for high/med competency species, 4 for low
design <- data.frame(
  species = spp,
  trayID = 1:length(spp),
  rep = 1
)


#MAKE TRAYS----
#interplanting distance is just 1.7cm

get_grid <- function(x){
  design_x <- design[x,]
  which_spp <- which(spp == design_x$species)
  grid <- sample_community(which_spp, 0, Dist)
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

#print maps, 2 per page.
plot_list <- lapply(grid_list, function(x){
  plot_maps(x) + scale_fill_manual(values = 'grey80')
})
pdf(here::here('post_quarantine_trials/figures/maps_controls.pdf'))
plot_list
dev.off()
