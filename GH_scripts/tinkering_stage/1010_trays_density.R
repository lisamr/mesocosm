#design 1010 tray densities
#need to figure out how many inds will be in a 1010 tray. look at both square and hexagonal grids
library(raster)
library(tidyverse)
library(ggplot2)

#create a raster layer for tray. 1010 tray is 25.4cm^2. 
width <- 25.4
tray <- raster(ncol=1, nrow=1, crs="+proj=utm +zone=1 +datum=WGS84", xmn=0, xmx=width, ymn=0, ymx=width) 
values(tray) <- 1
tray <- rasterToPolygons(tray) #convert to sp poly
plot(tray)

#create hexagon grid to overlay.
size <- 1
hex_points <- spsample(tray, type = "hexagonal", cellsize = size)
hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = size)
plot(tray, col = "grey50", axes = TRUE)
plot(hex_points, col = "black", pch = 20, cex = 0.5, add = T)
plot(hex_grid, border = "orange", add = T)

#square grid
sq_raster <- raster(tray, resolution=1)
values(sq_raster) <- 1
sq_grid <- rasterToPolygons(sq_raster)
#plot
plot(tray, col = "grey50", axes = TRUE)
plot(sq_grid, border = "orange", add = T)

#quantify number of points
length(hex_grid) 
length(sq_grid) 

#check out number of cells over a series of different interplanting distances
make_grid_sq <- function(r){
  sq_raster <- raster(tray, resolution=r)
  values(sq_raster) <- 1
  sq_grid <- rasterToPolygons(sq_raster)
  return(sq_grid)
}
make_grid_hex <- function(r){
  hex_points <- spsample(tray, type = "hexagonal", cellsize = r)
  hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = r)
  return(hex_grid)
}

#plot distance against # cells
Dist <- seq(1,2,by=.1)
sq <- sapply(Dist, function(x) length(make_grid_sq(x)))
hex <- sapply(Dist, function(x) length(make_grid_hex(x)))
dat <- data.frame(Dist, sq, hex)
dat2 <- pivot_longer(dat, cols=c(sq, hex), "grid", values_to = "ncells")
#plot
ggplot(dat2, aes(Dist, ncells, color=grid)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits=c(0,750)) 
