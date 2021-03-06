#design 1010 tray densities
#need to figure out how many inds will be in a 1010 tray. look at both square and hexagonal grids
rm(list=ls())
library(raster)
library(tidyverse)
library(ggplot2)

#create a raster layer for tray. planting area for 1010 tray is 9.5in^2. 
width <- 9.5*2.54
tray <- raster(ncol=1, nrow=1, crs="+proj=utm +zone=1 +datum=WGS84", xmn=0, xmx=width, ymn=0, ymx=width) 
values(tray) <- 1
tray <- rasterToPolygons(tray) #convert to sp poly

#create square grid
make_grid_sq <- function(r){
  sq_raster <- raster(tray, resolution=r)
  values(sq_raster) <- 1
  sq_grid <- rasterToPolygons(sq_raster)
  return(sq_grid)
}
#create hexagon grid
make_grid_hex <- function(r){
  hex_points <- spsample(tray, type = "hexagonal", cellsize = r, offset=c(0.01,0.01)) #setting square offset value ensures the same grid is drawn every time. 
  hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = r)
  return(hex_grid)
}

#plot grids
plot_grid <- function(grid, axes=T, ...){
  
  #get title information
  x <- coordinates(grid)[,1]
  dist <- x[2]-x[1]
  title <- paste(length(grid), "cells, dist=", dist)
  
  #plot
  plot(tray, col = "grey90", axes = F, main=title, ...)
  plot(grid, border = "royalblue4", add = T)
  points(coordinates(grid), pch=16, cex=.5, col='black')
  
  #put in axes
  if(axes==T){
    xs <- unique(sort(round(x, 4)))
    xs <- xs[seq(1, length(xs), by=2)] #get odds
    ys <- coordinates(grid)[,2] %>% 
      round(4) %>% sort %>% unique 
    axis(1, at=xs, labels=c(1:length(xs)), pos=0, las=2) #bottom
    axis(3, at=xs, labels=c(1:length(xs)), pos=width, las=2) #top
    axis(2, at=ys, labels=rev(LETTERS[c(1:length(ys))]), pos=0, las=2) #left
    axis(4, at=ys, labels=rev(LETTERS[c(1:length(ys))]), pos=width, las=2) #right
  }else{}
  
}
#plot_grid(make_grid_sq(1))
#plot_grid(make_grid_hex(1.4), axes = F)


#plot distance against # cells
Dist <- seq(1,3,by=.01) %>% round(2) #seq makes imprecise numbers, %in% needs perfect equality. 
sq <- sapply(Dist, function(x) length(make_grid_sq(x)))
hex <- sapply(Dist, function(x) length(make_grid_hex(x)))
dat <- data.frame(Dist, sq, hex)
dat2 <- pivot_longer(dat, cols=c(sq, hex), "grid", values_to = "ncells") %>% arrange(grid)
#plot
p1 <- ggplot(dat2, aes(Dist, ncells, color=grid)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits=c(0,750)) 
#p1+scale_y_log10()
plotly::ggplotly(p1)

#envision how the densities may change with the richness levels in your experiment
richness <- c(1,2,4,6)
distances <- c(2.85, 2.25, 1.75, 1.55) #interplanting distances in cm #2.9, 2.4, 2, 1.75

#summarize in df
names(richness) <- distances #index key for recoding
dat3 <- dat2 %>% 
  filter(Dist %in% distances) %>% 
  mutate(richness = recode(Dist, !!!richness))

#plot
ggplot(dat3, aes(richness, ncells, color=grid)) + geom_point() +geom_line() + scale_y_continuous(limits = c(0, max(dat3$ncells)))

#how many seeds is that for 1 whole series of densities? gonna need to order more seeds. for this tinkering around experiment, just do 3 distances. 
distances2 <- c(3, 2, 1.4)
dat4 <- dat2 %>% 
  filter(Dist %in% distances2)
dat4 %>% group_by(grid) %>% 
  summarise(cells_per_series = sum(ncells),
            cells_per_sp = cells_per_series*3, 
            cells_total = cells_per_sp*5,
            ninoculated = .1*cells_total) 

#Make templates
pdf("GH_plots/grid_templates.pdf")
#lapply(distances2, function(x) plot_grid(make_grid_sq(x)))
lapply(distances, function(x) plot_grid(make_grid_hex(x)))
dev.off()

#1010 trials planted on 1/14/20----
#arugula: 1.4, 2, 3, 1.4, 2, 3
#mustard: 1.4, 2, 3
#PC_F1: 1.4, 2, 3
#romaine pellets: 1.4, 2, 3
#romaine: 1.4
#radish: 1.4

#count how many individuals need to be inoculated. 3ish hours.
ntrays <- c(7, 5, 5)
dat4hex <- dat4 %>% 
  filter(grid=="hex")
(rep(dat4hex$ncells, times=ntrays) * .1) %>% 
  round() %>% 
  sum

#create data recording sheets. create df saying what you planted and then create the grids.
species <- c('arugula', 'mustard', 'PC_F1', 'romaine_pel', 'romaine', 'radish')
ntrays_sp <- c(6, 3, 3, 3, 1, 1)
dists <- c(1.4, 2, 3, 1.4, 2, 3, 1.4, 2, 3, 1.4, 2, 3, 1.4, 2, 3, 1.4, 1.4)
dat5 <- data.frame(species = rep(species, times=ntrays_sp), 
           distances = dists)

plot_grid_inoc <- function(r, ...){
  #data
  gridh <- make_grid_hex(r)
  #select which individuals get inoculated
  n_inoc <- round(length(gridh)*.1)
  samp <- sample(1:length(gridh), n_inoc)
  #get title information
  x <- coordinates(gridh)[,1]
  dist <- x[2]-x[1]
  title <- paste(length(gridh), "cells, dist=", dist)
  #plot
  plot(tray, col = "grey90", main=title, ...)
  plot(gridh, border = "royalblue4", add = T)
  points(coordinates(gridh[samp]), pch=16, cex=.5, col='black')
  #put in axes
  xs <- unique(sort(round(x, 4)))
  xs <- xs[seq(1, length(xs), by=2)] #get odds
  ys <- coordinates(gridh)[,2] %>% 
    round(4) %>% sort %>% unique 
  axis(1, at=xs, labels=c(1:length(xs)), pos=0, las=2) #bottom
  axis(3, at=xs, labels=c(1:length(xs)), pos=width, las=2) #top
  axis(2, at=ys, labels=rev(LETTERS[c(1:length(ys))]), pos=0, las=2) #left
  axis(4, at=ys, labels=rev(LETTERS[c(1:length(ys))]), pos=width, las=2) #right
}
plot_grid_inoc(1.4)

#print maps for 1010 trial
head(dat5)
pdf("GH_plots/1010trial_inoc_maps.pdf")
set.seed(1)
for(i in 1:nrow(dat5)){
  plot_grid_inoc(r=dat5$distances[i], sub=dat5$species[i])
}
dev.off()

#print maps for mom's cricket die cut
plot_grid_cricket <- function(r, ...){
  #data
  gridh <- make_grid_hex(r)
  #select which individuals get inoculated
  n_inoc <- round(length(gridh)*1)
  samp <- sample(1:length(gridh), n_inoc)
  #get title information
  x <- coordinates(gridh)[,1]
  dist <- x[2]-x[1]
  title <- paste(length(gridh), "cells, dist between dots=", dist, 'cm')
  #plot
  plot(tray, main=title, col='white', ...)
  plot(gridh, border = "black", add = T)
  points(coordinates(gridh[samp]), pch=16, cex=.5, col='black')
  #put in axes
  #xs <- unique(sort(round(x, 4)))
  #xs <- xs[seq(1, length(xs), by=2)] #get odds
  #ys <- coordinates(gridh)[,2] %>% round(4) %>% sort %>% unique 
  #axis(1, at=xs, labels=c(1:length(xs)), pos=0, las=2) #bottom
  #axis(3, at=xs, labels=c(1:length(xs)), pos=width, las=2) #top
  #axis(2, at=ys, labels=rev(LETTERS[c(1:length(ys))]), pos=0, las=2) #left
  #axis(4, at=ys, labels=rev(LETTERS[c(1:length(ys))]), pos=width, las=2) #right
}
dat3.1 <- filter(dat3, grid=='hex')
head(dat3.1)
pdf("GH_plots/cricket_diecut_border.pdf")
for(i in 1:nrow(filter(dat3, grid=='hex'))){
  plot_grid_cricket(r=dat3.1$Dist[i])
}
dev.off()



#
