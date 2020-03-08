#create plots for shawn and mom

#total numbers of plants
#1 1.det.add.1            85
#2 1.det.add.2           154
#3 1.det.add.4           238
#4 1.det.add.6           304
#5 1.det.sub.1           238

#ID            total  sp_1  sp_2  sp_3  sp_4  sp_5  sp_6
#<chr>         <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#1 1.det.add.1      85    85     0     0     0     0     0
#2 1.det.add.2     154    89    65     0     0     0     0
#3 1.det.add.4     238    86    61    50    41     0     0
#4 1.det.add.6     304    87    63    51    42    35    26
#5 1.det.sub.1     238   238     0     0     0     0     0
#6 1.det.sub.2     238   138   100     0     0     0     0
#7 1.det.sub.4     238    86    61    50    41     0     0
#8 1.det.sub.6     238    68    50    40    33    27    20

#horizontal vertical
#[1,]        2.8 2.424871
#[2,]        2.1 1.818653
#[3,]        1.7 1.472243
#[4,]        1.5 1.299038

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
plot_grid_cricket <- function(r, axes=F, ...){
  
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
  
  if(axes==T){
    #put in axes
    xs <- unique(sort(round(x, 4)))
    xs <- xs[seq(1, length(xs), by=2)] #get odds
    ys <- coordinates(gridh)[,2] %>% round(4) %>% sort %>% unique 
    axis(1, at=xs, labels=c(1:length(xs)), pos=0, las=2) #bottom
    axis(3, at=xs, labels=c(1:length(xs)), pos=width, las=2) #top
    axis(2, at=ys, labels=rev(LETTERS[c(1:length(ys))]), pos=0, las=2) #left
    axis(4, at=ys, labels=rev(LETTERS[c(1:length(ys))]), pos=width, las=2) #right
  }

}

plot_grid_cricket(2.8, axes = F)
#pdf('GH_plots/cricket_diecut/grids_final_densities.pdf')
lapply(1:length(Dist), function(x) plot_grid_cricket(Dist[x]))
#dev.off()

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

#plot actual size map for Mom----
spdf_list <- readRDS('GH_output/species_distributions/spdf_list_tighterdens.RDS')
spp <- c('radish', 'arugula', 'pac_choy', 'romaine', 'basil', 'clover')

#export dimensions need to be different for each density
#plot 21, 64, 107, 150 showing det/sub series
#plot 1, 41, 84, 127 showing det/add series

#tray 21, 2.8cm
plot_maps(spdf_list[[21]], plotted_points = c('S', 'C'), point_cex = 2)
ggsave('GH_plots/cricket_diecut/maps_actualsize_2.8cm.pdf', width = 32.135, height=32.135, units = 'cm') #the dimensions should be actual size. check with illustrator. 

#tray 64, 2.1cm
plot_maps(spdf_list[[64]], plotted_points = c('S', 'C'), point_cex = 2)
ggsave('GH_plots/cricket_diecut/maps_actualsize_2.1cm.pdf', width = 32.108, height=32.108, units = 'cm') #the dimensions should be actual size. check with illustrator. 

#tray 107, 1.7cm
plot_maps(spdf_list[[107]], plotted_points = c('S', 'C'), point_cex = 2)
ggsave('GH_plots/cricket_diecut/maps_actualsize_1.7cm.pdf', width = 31.565, height=31.565, units = 'cm') #the dimensions should be actual size. check with illustrator. 

#tray 150, 1.5cm
plot_maps(spdf_list[[150]], plotted_points = c('S', 'C'), point_cex = 2)
ggsave('GH_plots/cricket_diecut/maps_actualsize_1.5cm.pdf', width = 31.678, height=31.678, units = 'cm') #the dimensions should be actual size. check with illustrator. 

