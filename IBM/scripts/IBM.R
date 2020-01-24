#code for individual based model simulating rhizoctonia spread

#assume community filled in a hexagonal grid with r interplanting distance with s species, each with p competency. A portion of the communities will be inoculated and the epidemics will spread as a funciton of the community composition (species identity and density). States are C (challenged), S (susceptible), and I (infected). Transmission depends on number, identity, and distance from infected neighbors, as well as susceptibility of recipient host. 

#for now, haven't added in time varying transition rates, real compositions, and have guessed abou the disease parameters.

#load packages----
rm(list=ls())
library(raster)
library(tidyverse)
library(ggplot2)


#parameters----
#tray dimensions
width <- 9.5*2.54
length <- width 

#DESIGN 
pinoc <- .1 #percent inoculated at start of experiemnt
r <- 2 #interplanting distance
s <- 6 #number of species
spp <- c(paste0('sp_', rep(1:s)))

#DISEASE 
## transition from C -> S ##
#inoculum decay rate invariant of species.
delta <- 1/5 #1/average number of days inoc stays around

## transition from S -> I ##
#create a matrix of the beta_ij values. these will come from empirical data, but will just make them up for now. 
comp <- c(1, .5, .3, .2, 0, 0) #vector of relative "competencies"
names(comp) <- spp #named vector
beta_ij <- matrix(comp, nrow=length(comp), ncol = length(comp))
beta_ij <- beta_ij * t(beta_ij)
rownames(beta_ij) <- colnames(beta_ij) <- spp

## transition from C -> I ##
#rate of infection from inoculum to plant. varies by species. 
a <- comp 

#for plotting
pal <- RColorBrewer::brewer.pal(n = length(spp), name = "RdYlBu")
theme_set(theme_bw())

#create grid----
#First create hexagonal grid to figure out constraints.
#planting area for 1010 tray is 9.5in^2. 
tray <- raster(ncol=1, nrow=1, crs="+proj=utm +zone=1 +datum=WGS84", xmn=0, xmx=width, ymn=0, ymx=length) 
values(tray) <- 1
tray <- rasterToPolygons(tray) #convert to sp poly

#create hexagon grid
make_grid_hex <- function(r){
  hex_points <- spsample(tray, type = "hexagonal", cellsize = r, offset=c(0.1,0.1)) #setting square offset value ensures the same grid is drawn every time. 
  hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = r)
  return(hex_grid)
}

grid <- make_grid_hex(r)

#create communities----

#spatial dataframe will be more streamlined. will come from the design from your data. 
ninoc <- round(.1*length(grid))
grid_df <- SpatialPolygonsDataFrame(
  grid, data.frame(
    ID=row.names(grid),
    x=coordinates(grid)[,1],
    y=coordinates(grid)[,2],
    spID=sample(c('sp_1', 'sp_2'), length(grid), T),
    inoculated=sample(c(rep(1, ninoc), rep(0, length(grid)-ninoc)))
  ))

#quick plot
#plot(grid_df) + text(coordinates(grid_df), cex=.5, col=as.numeric(grid_df$spID))

#define agents dataframe----

#create dataframe that defines ID, coordinates, species ID, state, # infected neighbors, # total neighbors, % infected neighbors, neighbor ID vector. This dataframe will stay stationary through time for the most part, except for state and # and % infected neighbors.
agents <- grid_df@data

#pairwise adjacency matrix
adj <- rgeos::gTouches(grid_df, byid = T) 

#pairwise transmission matrix
#fill in the community grid with the beta_ij values
trans <- beta_ij[agents$spID, agents$spID] #get pairwise betas for all of the individuals given their species identity.
#change the names to their IDs
colnames(trans) <- rownames(trans) <- row.names(grid_df)

#NOTE: if you ever need to, you can create a pairwise distance decay matrix to account for declining transmission with distance. 


#make the dataframe
agents$state <- ifelse(agents$inoculated==1, "C", "S") #initial states
agents$adj <- sapply(1:nrow(agents), function(i) {
  x <- which(adj[row.names(agents)[i],]==T)
  unname(x)
  }) #neighbors

#count how many neighbors
agents$n_adj <- unlist(lapply(1:nrow(agents), function(x) length(agents$adj[[x]])))

#inoculum decay rate (transition from C -> S)
agents$C_to_S <- 1 - exp(-delta) #1/average number of days inoc stays around

#incolum to plant transmission rate. will vary by time. 
#C -> I through primary transmission (alpha)
rate <- a[agents$spID]
agents$C_to_Ia <- 1 - exp(-(rate))

#plant to plant transmission. will vary by time. 
#C -> I through secondary transmission (beta)
#rate is the secondary transmission function for each ID. infections happen at the rate: e^(-sum_j(beta_ij * (# infected_j) / (total # neighbors)))
#can just sum up all of the beta_ij's for all of the IDs which are infected and divide by total # neighbors. 
plant_to_plant <- function(i){
  beta_adj <- trans[agents$ID[i], agents$adj[[i]]] #beta's of the adjacent plants for a given ID
  adj_is_I <- agents$state[agents$adj[[i]]]=="I" #ask if neighbors are infected
  1 - exp(-1*sum(beta_adj*adj_is_I)/agents$n_adj[i]) #rate of transmission from S -> I 
}
agents$C_to_Ib <- sapply(1:nrow(agents), plant_to_plant)

#rate that plant stays C
agents <- agents %>% mutate(C_to_C = 1 - C_to_S - C_to_Ia - C_to_Ib)

head(agents)

#rules of spread----

#run epidemic----