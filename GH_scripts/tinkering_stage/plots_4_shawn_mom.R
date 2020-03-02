#create plots for shawn and mom

#total numbers of plants
#1 1.det.add.1            85
#2 1.det.add.2           154
#3 1.det.add.4           238
#4 1.det.add.6           304
#5 1.det.sub.1           238

#load source code----
rm(list=ls())
source('IBM/scripts/IBM_functions.R')

#parameters----
#distances wanted
Dist <- c(2.8, 2.1, 1.7, 1.5)

#create tray
#tray dimensions
width <- 9.5*2.54
length <- width 
tray <- make_tray(width, length) #contraining the hexes to 9.5 inches.
tray_bigger <- make_tray(9.75*2.54, 9.75*2.54) #actual size of the tray. 

#plot----
#customize `plot_maps` for the cricket.
plot_grid_cricket <- function(r, ...){
  
  make_grid <- function(r, tray){
    #tray=output from `make_tray()`
    hex_points <- spsample(tray, type = "hexagonal", cellsize = r, offset=c(0.1,0.1)) #setting square offset value ensures the same grid is drawn every time. 
    hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = r)
    return(hex_grid)
  }
  #data
  gridh <- make_grid(r, tray)
  #select which individuals get inoculated
  n_inoc <- round(length(gridh)*1)
  samp <- sample(1:length(gridh), n_inoc)
  #get title information
  x <- coordinates(gridh)[,1]
  dist <- x[2]-x[1]
  title <- paste(length(gridh), "cells, dist between dots=", dist, 'cm')
  
  #plot
  plot(tray_bigger, main=title, col='white', ...)
  plot(gridh, border = "black", add = T)
  points(coordinates(gridh[samp]), pch=16, cex=1.5, col='black')
  
  #put in axes
  #xs <- unique(sort(round(x, 4)))
  #xs <- xs[seq(1, length(xs), by=2)] #get odds
  #ys <- coordinates(gridh)[,2] %>% round(4) %>% sort %>% unique 
  #axis(1, at=xs, labels=c(1:length(xs)), pos=0, las=2) #bottom
  #axis(3, at=xs, labels=c(1:length(xs)), pos=width, las=2) #top
  #axis(2, at=ys, labels=rev(LETTERS[c(1:length(ys))]), pos=0, las=2) #left
  #axis(4, at=ys, labels=rev(LETTERS[c(1:length(ys))]), pos=width, las=2) #right
}

pdf('GH_plots/cricket_diecut/grids_final_densities.pdf')
lapply(1:length(Dist), function(x) plot_grid_cricket(Dist[x]))
dev.off()

#vertical and horizontal distances between points
make_grid <- function(r, tray){
  #tray=output from `make_tray()`
  hex_points <- spsample(tray, type = "hexagonal", cellsize = r, offset=c(0.1,0.1)) #setting square offset value ensures the same grid is drawn every time. 
  hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = r)
  return(hex_grid)
}
spacing <- function(r){
  coords <- coordinates(make_grid(r, tray))
  horizontal <- coords[2,1] - coords[1,1] #horizontal spacing
  ys <- unique(round(coords[,2], 6))
  vertical <- ys[2] - ys[1] #vertical spacing
  results <- c(horizontal, vertical)
  names(results) <- c('horizontal', 'vertical')
  return(results)
}

sapply(1:length(Dist), function(x) spacing(Dist[x])) %>% t

#height should be vertical*sqrt(3)/2
height <- function(x) x*sqrt(3)/2
height(Dist) #agrees :)
 




