#create grids for experiment with 1010 trays and hexagonal layout.
#collapse all chunks: option command O

#load packages----
rm(list=ls())
library(raster)
library(tidyverse)
library(ggplot2)

#define functions----
#create hexagon grid
make_grid_hex <- function(r){
  hex_points <- spsample(tray, type = "hexagonal", cellsize = r, offset=c(0.1,0.1)) #setting square offset value ensures the same grid is drawn every time. 
  hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = r)
  return(hex_grid)
}

#create experimental design given a log-normal distribution
get_design <- function(sd, sims=10000){
  
  #assume relative abundances come from a log-normal distribution. Will randomly draw from a lognormal 10000 times and take the mean relative proportions.
  nspp <- max(richness) #number of species
  set.seed(0)
  relabund <- function(x, sd=.5) {
    y <- sort(rlnorm(1:nspp, meanlog = 1, sdlog = sd), decreasing = T)
    y/sum(y) #relative proportion
  } 
  
  #get average relative proportions
  y <- replicate(sims, relabund(sd = sd))
  rel_prop <- apply(y, 1, mean) 
  
  #calculate number of individuals for each treatment 
  M <- matrix(0, nrow=nrow(design), ncol=length(rel_prop))
  for(i in 1:nrow(design)){
    sp <- design$rich[i]
    y <- c(rel_prop[1:sp], rep(0, length(rel_prop)-sp)) #relative abundances given a certain number of species 
    n1 <- design$n[i]+sp #buffering the number of species so there aren't errors due to rounding
    N <- round(n1 * (y/sum(y))) #slightly greater than I want
    x <- table(sample(1:sp, sum(N)-design$n[i], replace = T, prob = N[1:sp]))#randomly remove excess individuals
    M[i,] <- N-c(unname(x), rep(0,length(N)-length(x)))
  }
  
  #if interplanting distances the same, make one of the additive rows the same as the sub row. differences due to sampling. 
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
  sp <- paste0("sp_", 1:nspp)
  names(sp) <- comp
  #create design with treatments, competencies, & species
  design4 <- suppressWarnings(
    design3 %>% 
      group_by(rep, rand, dens) %>% 
      mutate(
        comp = case_when( 
          rand=="det" ~ comp[order],
          rand=="stoch" ~ sample(comp, replace = F)[order]),
        species = as.factor(recode(comp, !!!sp)) ) %>% 
      mutate(comp = round(comp, 1), 
             SD = sd, 
             ID = interaction(rep, rand, dens, richness))
  ) 
  
  #get design!
  return(design4)
}

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
  ninoc_tot <- round(pinoc*length(SP))
  ninoc_sp <- floor(pinoc*tmp$nind)
  x <- ninoc_tot - sum(ninoc_sp) #remainder
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


#define parameters----
#tray dimensions
width <- 9.5*2.54
length <- width 

#initial UNCHANGING design pars
richness <- c(1,2,4,6)
nreps <- 1
comp <- c(1, .3, .2, .1, .001, 0) #vector of relative "competencies". have to differentiate sp5 and 6 somehow.
pinoc <- .1 #percent inoculated at start of experiemnt

#for plotting
pal <- RColorBrewer::brewer.pal(n = max(richness), name = "RdYlBu")

theme_set(theme_bw())

#create design skeleton----

#First create hexagonal grid to figure out constraints.
#planting area for 1010 tray is 9.5in^2. 
tray <- raster(ncol=1, nrow=1, crs="+proj=utm +zone=1 +datum=WGS84", xmn=0, xmx=width, ymn=0, ymx=length) 
values(tray) <- 1
tray <- rasterToPolygons(tray) #convert to sp poly

#simulate species distribution----

#not very streamlined, but oh well...
#need to play with interplanting distances. 

#iteration1
#interplanting distances
Dist <- c(2.4, 1.95, 1.55, 1.4)#additive assembly
Sub <- data.frame(dens="sub", Dist=1.75, ncells=224, richness)#sub assembly
#get number of plants with each interplanting distance
hex <- sapply(Dist, function(x) length(make_grid_hex(x)))
dat <- data.frame(Dist, ncells=hex, richness)
#add in density treatments. additive will increase linearly; substitutive will stay the same.
design <- rbind(cbind(dens="add", dat), Sub)
design1 <- get_design(sd=.5)
#filter only what you want
design1 <- design1 %>% 
  ungroup() %>% 
  filter(rand=='det', dens=='sub'|(dens=='add' & richness==1)) %>% 
  mutate(ID = droplevels.factor(ID),
         trayID = as.numeric(ID)) 

#iteration2
#interplanting distances
Dist2 <- c(2.7, 1.95, 1.55, 1.75) #additive assembly
Sub2 <- data.frame(dens="sub", Dist=1.4, ncells=340, richness)#sub assembly
#get number of plants with each interplanting distance
hex <- sapply(Dist, function(x) length(make_grid_hex(x)))
dat <- data.frame(Dist, ncells=hex, richness)
#add in density treatments. additive will increase linearly; substitutive will stay the same.
design <- rbind(cbind(dens="add", dat), Sub2)
design2 <- get_design(sd=.5)
#filter only what you want
design2 <- design2 %>% 
  ungroup() %>% 
  filter(rand=='det', dens=='sub'|(dens=='add' & richness==1)) %>% 
  mutate(ID = droplevels.factor(ID),
         trayID = as.numeric(ID)+max(design1$trayID)) 


#iteration3
#interplanting distances
Dist3 <- c(3, 1.95, 1.55, 1.4) #additive assembly
#get number of plants with each interplanting distance
hex <- sapply(Dist, function(x) length(make_grid_hex(x)))
dat <- data.frame(Dist, ncells=hex, richness)
#add in density treatments. additive will increase linearly; substitutive will stay the same.
design <- rbind(cbind(dens="add", dat), Sub2)
design3 <- get_design(sd=.5)
#filter only what you want
design3 <- design3 %>% 
  ungroup() %>% 
  filter(rand=='det', dens=='add' & richness==1) %>%
  mutate(ID = droplevels.factor(ID),
         trayID = as.numeric(ID)+max(design2$trayID)) 

#bind all of them together
finaldesign <- suppressWarnings(bind_rows(design1, design2, design3)) 

#create maps----
#create spatial maps to allow you to sow and monitor your trays.

#splits df by trays
dflist <- suppressWarnings(split(finaldesign, f = finaldesign$trayID )) %>% 
  #creates empty tibbles when treatment doesnt exist. remove.
  discard(function(x) nrow(x) == 0)

#create spatial polygons dataframe for each tray. takes some time.
spdf_list <- lapply(1:length(dflist), dflist_to_sp)

#plot maps----

#plot maps of trays with ggplot
names(dflist) #ID of the maps
plot_maps <- function(i){ #i is which tray
  tmp <- spdf_list[[i]] #ex/157 "0.5.9.det.add.6"   
  #make df readable to ggplot
  # add to data a "ID" column for each feature
  tmp$id <- rownames(tmp@data)
  # create a data.frame from our spatial object
  tmp_df <- fortify(tmp, region = "id") %>% 
    #merge the "fortified" data with attribute data
    merge(., tmp@data, by = "id")
  
  #add in color id for the species so its the same every time
  Colors <- pal
  names(Colors) <- levels(tmp_df$spID)
  
  #make a df of centroids to plot inoculated plants. 
  tmp_centroids <- data.frame(coordinates(tmp), state0=tmp$state0)
  
  #calculate axis labels
  xs <- unique(sort(round(tmp_centroids$X1, 4)))
  xs <- xs[seq(1, length(xs), by=2)] #get odds
  ys <- tmp_centroids$X2 %>% round(4) %>% sort %>% unique 
  planting_dist <- xs[2]-xs[1]
  
  #plot in ggplot!
  p1 <- ggplot(data = tmp_df, aes(x=long, y=lat)) +
    geom_polygon(aes(group = group, fill = spID))  +
    geom_path(aes(group = group), color = "white") +
    geom_point(data=tmp_centroids, aes(X1, X2), 
               alpha=ifelse(tmp_centroids$state0=="C", 1, 0)) +
    coord_equal() +
    scale_fill_manual(values=Colors) +
    labs(title = paste(paste0('Tray', tmp_df$trayID[1], "-"), tmp_df$rand[1], tmp_df$dens[1], paste0('replicate', tmp_df$rep[1]), sep = '-'),
         subtitle = paste("SD =", tmp_df$SD, ', nplants = ', length(tmp), ', spacing = ', planting_dist, 'cm')) +
    scale_x_continuous(name='', breaks=xs, labels=1:length(xs), sec.axis = dup_axis()) +
    scale_y_continuous(name='', breaks=ys, labels=rev(LETTERS[1:length(ys)]), sec.axis = dup_axis())
  
  return(p1)
}

#check out a tray
plot_maps(9)

#export----

#design dataframe 
write.csv(finaldesign, 'GH_output/tinkering_stage/1010_trial1.csv', row.names = F)

#spatial maps to physically print
pdf('GH_plots/tinkering_stage/1010_trial1_maps.pdf')
lapply(1:length(spdf_list), plot_maps) #run this to print the maps
dev.off()

