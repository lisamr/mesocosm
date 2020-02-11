#create grids for experiment with 1010 trays and hexagonal layout.
#collapse all chunks: option command O

#load packages----
rm(list=ls())
library(raster)
library(tidyverse)
library(ggplot2)

#define parameters----
#tray dimensions
width <- 9.5*2.54
length <- width 

#initial design pars
richness <- c(1,2,4,6)
Dist <- c(2.85, 2.25, 1.75, 1.55) #interplanting distances
Sub <- data.frame(dens="sub", Dist=1.75, ncells=224, richness)#values informed by hexagon grid below.
nreps <- 10
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

#create hexagon grid
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


#simulate species distribution----

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

design_orig <- get_design(sd=.5)

#decouple richness and Community competency----

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

#create maps----
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

#create spatial polygons dataframe for each tray. takes a few minutes
spdf_list <- lapply(1:length(dflist), dflist_to_sp)


#plot figures----

#plot density vs richness for both treatments.
ggplot(design, aes(richness, ncells, color=dens)) + 
  geom_point() + 
  geom_line() + 
  scale_y_continuous(limits=c(0,400)) 
ggsave('GH_plots/species_distributions/density_richness.pdf')

#change rand and dens names for plotting
rand.labs <- c('deterministic', 'stochastic')
names(rand.labs) <- c('det', 'stoch')
dens.labs <- c("additive", "substitutive")
names(dens.labs) <- c("add", "sub")

#figure of design with 4 treatments
ggplot(filter(design_orig, rep==2), aes(richness, nind, group=fct_rev(species), fill=species))+
  geom_col()+
  facet_wrap(~rand+dens, labeller = labeller(rand=rand.labs, dens=dens.labs))+
  scale_fill_manual(values=pal)
ggsave('GH_plots/species_distributions/design_4trtments.pdf', width = 7, height = 4.5)

#plot community competency against richness
#plot 1
design_augmented %>% 
  group_by(rep, rand, dens, richness, SD) %>% 
  summarise(CC = sum(nind*comp), 
            avgCC = sum(nind*comp)/sum(nind)) %>% 
  ggplot(., aes(richness, avgCC, color=as.factor(ifelse(SD==.5, "original", "augmented")))) +
  geom_point(alpha=.5) +
  labs(color="dataset", y='average host competency')
ggsave('GH_plots/species_distributions/avghostcomp_richness.pdf')
#plot2
design_augmented %>% 
  group_by(rep, rand, dens, richness, SD) %>% 
  summarise(CC = sum(nind*comp), 
            avgCC = sum(nind*comp)/sum(nind)) %>% 
  ggplot(., aes(richness, CC, color=as.factor(ifelse(SD==.5, "original", "augmented")))) +
  geom_point(alpha=.5) +
  labs(color="dataset", y='Total community competency')
ggsave('GH_plots/species_distributions/comcomp_richness.pdf')


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
plot_maps(17)

#plot 14 to 17 showing det/sub series
blank_theme <- theme(title = element_blank(), legend.position="none", axis.text = element_blank(), axis.ticks = element_blank())
tray14 <- plot_maps(14)+ blank_theme
tray15 <- plot_maps(15)+ blank_theme
tray16 <- plot_maps(16)+ blank_theme
tray17 <- plot_maps(17)+ blank_theme
cowplot::plot_grid(tray14, tray15, tray16, tray17)
ggsave('GH_plots/species_distributions/four_maps_detsub.pdf', width = 7, height = 7)

#plot 26 to 29 showing det/add series
tray26 <- plot_maps(26)+ blank_theme
tray27 <- plot_maps(27)+ blank_theme
tray28 <- plot_maps(28)+ blank_theme
tray29 <- plot_maps(29)+ blank_theme
cowplot::plot_grid(tray26, tray27, tray28, tray29)
ggsave('GH_plots/species_distributions/four_maps_detadd.pdf', width = 7, height = 7)


#export----

#design dataframe 
write.csv(design_augmented, 'GH_output/species_distributions/design_augmented.csv', row.names = F)


#list of spatial dataframes
saveRDS(spdf_list, 'GH_output/species_distributions/spdf_list.RDS')
#open with `readRDS('filename')`

#spatial maps to physically print
#lapply(1:length(spdf_list), plot_maps) #run this to print the maps


