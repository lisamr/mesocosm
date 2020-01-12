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
plot(tray)

#create square grid
make_grid_sq <- function(r){
  sq_raster <- raster(tray, resolution=r)
  values(sq_raster) <- 1
  sq_grid <- rasterToPolygons(sq_raster)
  return(sq_grid)
}
#create hexagon grid
make_grid_hex <- function(r){
  hex_points <- spsample(tray, type = "hexagonal", cellsize = r)
  hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = r)
  return(hex_grid)
}

#plot grids
plot_grid <- function(grid){
  x <- coordinates(grid)[,1]
  dist <- x[2]-x[1]
  title <- paste(length(grid), "cells, dist=", dist)
  plot(tray, col = "grey90", axes = F, main=title)
  plot(grid, border = "royalblue4", add = T)
  points(coordinates(grid), pch=16, cex=.5, col='black')
}
plot_grid(make_grid_sq(1))
plot_grid(make_grid_hex(1))


#plot distance against # cells
Dist <- seq(1,3,by=.05) %>% round(2) #seq makes imprecise numbers, %in% needs perfect equality. 
sq <- sapply(Dist, function(x) length(make_grid_sq(x)))
hex <- sapply(Dist, function(x) length(make_grid_hex(x)))
dat <- data.frame(Dist, sq, hex)
dat2 <- pivot_longer(dat, cols=c(sq, hex), "grid", values_to = "ncells") %>% arrange(grid)
#plot
p1 <- ggplot(dat2, aes(Dist, ncells, color=grid)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits=c(0,750)) 
p1+scale_y_log10()
plotly::ggplotly(p1)

#envision how the densities may change with the richness levels in your experiment
richness <- c(1,2,4,6)
distances <- c(3, 2.4, 1.8, 1.5) #interplanting distances in cm

#summarize in df
names(richness) <- distances #index key for recoding
dat3 <- dat2 %>% 
  filter(Dist %in% distances) %>% 
  mutate(richness = recode(Dist, !!!richness))

#plot
ggplot(dat3, aes(richness, ncells, color=grid)) +
  geom_point() +
  geom_line()

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
lapply(distances2, function(x) plot_grid(make_grid_sq(x)))
lapply(distances2, function(x) plot_grid(make_grid_hex(x)))
dev.off()






#
