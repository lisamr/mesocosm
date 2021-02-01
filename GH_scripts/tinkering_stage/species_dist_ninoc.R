#create grids for experiment with 1010 trays and hexagonal layout.
#collapse all chunks: option command O

#LOAD PACKAGES----
rm(list=ls())
library(raster)
library(tidyverse)
library(ggplot2)

#DEFINE PARAMETERS----
set.seed(0) #same results every time you redo your simulation
#tray dimensions
width <- 9.5*2.54
length <- width 

#initial design pars
richness <- c(1,2,4,6)
#interplanting distances
Dist <- c(2.8, 2.1, 1.7, 1.5) #in between. Go with this.
#Dist <- c(2.55, 1.95, 1.55, 1.4) #tighter
#Dist <- c(2.85, 2.25, 1.75, 1.55) #original
Sub <- data.frame(dens="sub", Dist=1.7, ncells=238, richness)#values informed by hexagon grid below.
nreps <- 10
comp <- c(1, .3, .2, .1, .001, 0) #vector of relative "competencies". have to differentiate sp5 and 6 somehow.
nspp <- length(comp)
spp <- c('radish', 'arugula', 'basil', 'green', 'red', 'butter')
#spp <- paste0('sp_', 1:nspp) #use for the maps for mom.
ninoc <- 10 #number inoculated at start of experiemnt

#for plotting
pal <- RColorBrewer::brewer.pal(n = max(richness), name = "RdYlBu")
names(pal) <- spp

theme_set(theme_bw())

#CREATE DESIGN SKELETON----

#First create hexagonal grid to figure out constraints.
#planting area for 1010 tray is 9.5in^2. 
tray <- raster(ncol=1, nrow=1, crs="+proj=utm +zone=1 +datum=WGS84", xmn=0, xmx=width, ymn=0, ymx=length) 
values(tray) <- 1
tray <- rasterToPolygons(tray) #convert to sp poly

#create hexagon grid (can also use the already defined function `make_grid`)
make_grid_hex <- function(r){
  hex_points <- spsample(tray, type = "hexagonal", cellsize = r, offset=c(0.1,0.1)) #setting square offset value ensures the same grid is drawn every time. 
  hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = r)
  return(hex_grid)
}

#get number of plants with each interplanting distance
hex <- sapply(Dist, function(x) length(make_grid_hex(x)))
dat <- data.frame(Dist, ncells=hex, richness)

#add in density treatments. additive will increase linearly; substitutive will stay the same.
design <- rbind(cbind(dens="add", dat), Sub)



#STOCHASTIC SPECIES ORDER----
#instead of randomly creating different species orders for stochastic treatments, should do a random stratified order. that way, with only 20 (10 of each densities) replicates, the broadest range of orders are represented. 

#matrix of species orders
order_mat <- matrix(NA, nrow = nreps*2, ncol = length(comp))

#figure out first species for each replicate
y <- floor(2*nreps/nspp) #minimum number of times each species is the first species 
x <- rep(1:nspp, y)
extrax <- sample(1:nspp, nrow(order_mat)-length(x), F) #some species get to start an extra time. this is because nreps may not be a multiple of nspp. 
order_mat[,1] <- c(x, extrax)

#fill in the rest of the columns
for(i in 1:nrow(order_mat)){
  firstsp <- order_mat[i,1]
  nextsp <- sample(c(1:nspp)[-firstsp])
  order_mat[i,2:ncol(order_mat)] <- nextsp
}

#reshuffle the rows where first 2 species are the same.
reshuffle_rows <- 1 #initiate the vector to keep track of while loop. make sure its a number so it can be a T/F, i.e. no NA
while(sum(reshuffle_rows)!=0){ #repeats loop until condition met
  
  #all the pairwise combinations to compare
  pairs <- combn(nrow(order_mat), 2) 
  
  #do lots of pairwise comparisons of first 2 species
  reshuffle <- rep(NA, ncol(pairs))
  for(i in 1:ncol(pairs)){
    A <- pairs[1,i]
    B <- pairs[2,i]
    reshuffle[i] <- order_mat[A,1]==order_mat[B,1] & order_mat[A,2]==order_mat[B,2] #compare them
  }
  
  #reshuffle the bottom row. then check again.
  reshuffle_rows <- pairs[2,reshuffle] 
  for(i in reshuffle_rows){
    firstsp <- order_mat[i,1]
    nextsp <- sample(c(1:nspp)[-firstsp])
    order_mat[i,2:ncol(order_mat)] <- nextsp
  }
}

order_mat #species order for stochastic treatments :)

#orders to be added to design dataframe
stochastic_order <- as.vector(apply(t(order_mat), 2, rep, 4))
deterministic_order <- rep(1:6, length.out = length(stochastic_order))
competencies <- comp[c(deterministic_order, stochastic_order)]

#SIMULATE SPECIES DISTRIBUTION----

#assume relative abundances come from a log-normal distribution. Will randomly draw from a lognormal 10000 times and take the mean relative proportions.
nspp <- max(richness) #number of species
sims <- 10000 #number of simulations
relabund <- function(x, sd=.5) {
  y <- sort(rlnorm(1:nspp, meanlog = 1, sdlog = sd), decreasing = T)
  y/sum(y) #relative proportion
} 


#create experimental design given a log-normal distribution
get_design <- function(sd){
  
  #get average relative proportions
  y <- replicate(sims, relabund(sd = sd))
  rel_prop <- apply(y, 1, mean) 
  
  #calculate number of individuals for each treatment 
  set.seed(0)
  M <- matrix(0, nrow=nrow(design), ncol=length(rel_prop))
  for(i in 1:nrow(design)){
    sp <- design$rich[i]
    y <- c(rel_prop[1:sp], rep(0, length(rel_prop)-sp)) #relative abundances given a certain number of species 
    n1 <- design$n[i]+sp #buffering the number of species so there aren't errors due to rounding
    N <- round(n1 * (y/sum(y))) #slightly greater than I want
    x <- table(sample(1:sp, sum(N)-design$n[i], replace = T, prob = N[1:sp]))#randomly remove excess individuals
    M[i,] <- N-c(unname(x), rep(0,length(N)-length(x)))
  }
  
  #make one of the additive rows the same as the sub row. differences due to sampling. 
  sums <- apply(M, 1, sum)
  whichrow <- which(abs(sapply(1:length(richness), function(i) sums[length(richness)+1]-sums[i])) < 5)
  M[whichrow,] <- M[whichrow + length(richness), ]
  
  #bind M to design df
  M <- as.data.frame(M)
  design2 <- cbind(design, M) %>% 
    reshape2::melt(., id.vars = c("dens", "Dist",  "richness"), measure.vars = paste0("V", 1:nspp), variable.name = "order", value.name = "nind") %>% 
    arrange(dens, richness, Dist)
  
  #create design for each tray
  #4 treatments: add/det, sub/det, add/stoch, sub/stoch, each replicated 10 times. 
  #determ: competency decreases with order. 
  #stoch: competency varies independently of order
  rand <- c("det", "stoch")
  design3 <- cbind(
    rep=rep(1:nreps, each=nrow(design2)*length(rand)), 
    rand=rep(rand, each=nrow(design2)),
    design2)
  
  #assigning competency values that vary deterministically or randomly (stochastic) with order.
  #sp <- paste0("sp_", 1:nspp)
  sp <- spp
  names(sp) <- comp
  
  
  #create design with treatments, competencies, & species
  design4 <- design3 %>% 
    arrange(rand, dens, rep, richness) %>% 
    mutate(comp = competencies, 
           species = as.factor(recode(comp, !!!sp))) %>%
    mutate(comp = round(comp, 1), 
           SD = sd, 
           ID = interaction(rep, rand, dens, richness))
  
  design4$species <- fct_relevel(design4$species, spp)
  
  
  
  
  #create design with treatments, competencies, & species
  #  design4 <- suppressWarnings(
  #    design3 %>% 
  #     group_by(rep, rand, dens) %>% 
  #    mutate(
  #     comp = case_when( 
  #      rand=="det" ~ comp[order],
  #     rand=="stoch" ~ sample(comp, replace = F)[order]),
  #  species = as.factor(recode(comp, !!!sp)) ) %>% 
  #  mutate(comp = round(comp, 1), 
  #        SD = sd, 
  #       ID = interaction(rep, rand, dens, richness))) 
  
  #get design!
  return(design4)
}

design_orig <- get_design(sd=.5)

#figure of design with 4 treatments
ggplot(filter(design_orig, rep==4), aes(richness, nind, group=fct_rev(species), fill=species))+ geom_col()+ facet_wrap(~rand+dens)+ scale_fill_manual(values=pal)

#check out change in numbers
design_orig %>% group_by(ID) %>% pivot_wider(names_from = species, values_from = nind) %>% dplyr::select(-order, -comp, -SD, -rep, -dens)

design_orig %>% group_by(rep, rand, dens, richness) %>% 
  summarise(total = sum(nind), 
            sp_1 = nind[species==spp[1]], 
            sp_2 = nind[species==spp[2]],
            sp_3 = nind[species==spp[3]],
            sp_4 = nind[species==spp[4]],
            sp_5 = nind[species==spp[5]],
            sp_6 = nind[species==spp[6]]) 

#DECOUPLE RICHNESS AND COMPETENCY----

#redo what you did above, but with a different sd. then filter out just what you need.
filter_design <- function(sd){
  design <- get_design(sd)
  #find which ones are high in CC and richness. 
  keep <- design %>% 
    get_CC_data %>% 
    filter(CC>.5, richness>1, !duplicated(richness)) %>% 
    mutate(ID = interaction(rep, rand, dens, richness))
  #filter those ones out, which can be added to the final design.
  filtered <- design %>% filter(ID %in% keep$ID)
  
  return(filtered)
}

#additional trays
get_CC_data <- function(design){
  design %>% group_by(rep, rand, dens, richness) %>% summarise(CC = sum(nind*comp)/sum(nind)) %>% ungroup()
}
sd5 <- filter_design(5)
sd2 <- filter_design(2)
sd1 <- filter_design(1)

#look at their summaries
#get_CC_data(sd5)
#get_CC_data(sd2)
#get_CC_data(sd1)

#augment to the original design
design_augmented <- bind_rows(design_orig, sd5, sd2, sd1)
design_augmented$trayID <- as.numeric(interaction(design_augmented$SD, design_augmented$ID, drop = T))

#CREATE MAPS----
#create spatial maps to allow you to sow and monitor your trays.

#splits df by trays
dflist <- suppressWarnings(split(design_augmented, f = with(design_augmented, list(SD, ID)))) %>% 
  #creates empty tibbles when treatment doesnt exist. remove.
  discard(function(x) nrow(x) == 0)

dflist_to_sp <- function(i){
  
  #turn df to spatial map, adding attribute data along the way.
  tmp <- dflist[[i]] %>% arrange(species)
  SP <- make_grid_hex(tmp$Dist[1])
  SPdf <- SpatialPolygonsDataFrame(
    SP, data.frame(
      ID=row.names(SP),
      x=coordinates(SP)[,1],
      y=coordinates(SP)[,2],
      trayID = rep(tmp$trayID[1], length(SP)),
      rand = rep(tmp$rand[1], length(SP)),
      dens = rep(tmp$dens[1], length(SP)),
      richness = rep(tmp$richness[1], length(SP)),
      rep = rep(tmp$rep[1], length(SP)),
      SD = rep(tmp$SD[1], length(SP)),
      spID = rep(tmp$species, times=tmp$nind),
      comp = rep(tmp$comp, times=tmp$nind))) 
  
  #inoculate subset of individuals (strictly proportional to relative abundance)
  #discrete numbers only. ensure total inoculuated is the same every time.
  ninoc_sp <- floor(ninoc/length(SP)*tmp$nind)
  x <- ninoc - sum(ninoc_sp) #remainder
  x <- sample(tmp$species, x, replace = T, prob = tmp$nind) %>% table() 
  ninoc_sp <- data.frame(ninoc_sp + x)
  
  #get vector for inoculated or not (Challenged or Susceptible)
  state0 <- unlist(sapply(1:nrow(ninoc_sp), function(i) c(rep("C", ninoc_sp$Freq[i]), rep("S", tmp$nind[i]-ninoc_sp$Freq[i])) )) 
  SPdf$state0 <- state0
  
  #randomize locations of species. just randomize last three columns
  randomized <- SPdf@data[sample(1:length(SP)),(ncol(SPdf)-3):ncol(SPdf)]
  SPdf@data <- cbind(SPdf@data[,1:(ncol(SPdf)-4)], randomized)
  
  return(SPdf)
}

#create spatial polygons dataframe for each tray. takes a few minutes
spdf_list <- lapply(1:length(dflist), dflist_to_sp)

names(dflist) #ID of the maps

#EXPORT----

#design dataframe 
#write.csv(design_augmented, 'GH_output/species_distributions/design_augmented.csv', row.names = F)

#list of spatial dataframes
saveRDS(spdf_list, 'GH_output/species_distributions/spdf_list_ninoc.RDS')
