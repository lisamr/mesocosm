#IBM tutorial. "IBM_functions.R" is the source code for all the functions. Here, I'll put those to use, but the code will be cleaner (won't have the functions defined here).

#using sample community composition

rm(list=ls())

#load source code----
source('IBM/scripts/IBM_functions.R')

#parameters----
#tray dimensions
width <- 9.5*2.54
length <- width 

#DESIGN (keep the names the same. simplifies downstream functions)
pinoc <- .1 #percent inoculated at start of experiemnt
s <- 6 #max number of species
spp <- c(paste0('sp_', rep(1:s))) #names of species
tfinal <- 20 #how many time steps
comp <- c(1, .5, .3, .2, 0, 0) #vector of relative "competencies"

#DISEASE 
#transmission curves
#mean values approximated from Otten et al. (2003)
beta_curve <- function(t) .2*exp(-3*(log(t/11))^2)
alpha_curve <- function(t) .4*(1-.3)^t

#distance decay in probability of secondary transmission from Kleczkowski et al. (1997)
dist_decay <- function(x, a=14.9, d=.36, tau=3.73, sigma=.00099){
  C <- a*d*exp(d*tau)
  y <- C*exp(-sigma*x^2)
  #standardize to 0-1, relative to a distance of 0
  y0 <- C*exp(-sigma*0^2)
  y/y0
}

### transmission from C -> S ###
#inoculum decay rate invariant of species.
delta <- 1/5 #1/average number of days inoc stays around

### transmission from S -> I ###
#create a matrix of the beta_ij values. Assume that pairwise values will control the amplitude of beta(t). these will come from empirical data, but will just make them up for now.

#create matrix of amplitudes of the beta_ij 
beta_ij_t <- make_beta_ij_t(comp)

### transmission from C -> I ###
#rate of infection from inoculum to plant. varies by species and time
alpha_i_t <- make_alpha_i_t(comp)

#create community----

#I created a sample community in the IBM source code. It's a simple 2 species community. However, you'll likely want to make your own. 
sample_grid <- sample_community(which_spp = 2, perc_inoc = .01, planting_dist = 1)



#run epidemic----

#for every individual, if state==C, then C->C, C->S, or C->I; if state==S, then S->S, or S->I; if state==I, then always stays I. Each transition has a probability and the fate of the transition determined by a value drawn from a uniform distribution between 0 and 1. Will need to loop through every individual every time step. At each time step, record the state of every individual. Output should be a matrix of states, with rows equalling # individuals and cols equalling # time steps. 

#try out transmission with nearest neighbor only infections and ones with distance decay (kernel) on two communities with different densities

#communities
samplegrid1 <- sample_community(1, .1, 2.55) #species 1
samplegrid2 <- sample_community(2, .1, 2.55) #species 2

#simulate!
testrunNN1 <- IBM(samplegrid1, Type = "NN")
testrunNN2 <- IBM(samplegrid2, Type = "NN")
testrunKernel1 <- IBM(samplegrid1, Type = "Kernel", spatialdecay = .002)
testrunKernel2 <- IBM(samplegrid2, Type = "Kernel", spatialdecay = .002)

#view
head(testrunNN1); head(testrunKernel1)
head(testrunNN2); head(testrunKernel2)

#visualize----

#first plot summary of S and I
p1 <- plotS_I(testrunNN1)[[2]]
p2 <- plotS_I(testrunNN2)[[2]]
p3 <- plotS_I(testrunKernel1)[[2]]
p4 <- plotS_I(testrunKernel2)[[2]]
cowplot::plot_grid(p1, p2, p3, p4, labels = c('NN1', 'NN2', 'Kernel1', 'Kernel2'))


#now plot spatial map of the spread (animation is ~3 minutes)
plot_spread_map(samplegrid1, testrunNN1, animate = F)
plot_spread_map(samplegrid1, testrunKernel1, animate = F)
plot_spread_map(samplegrid2, testrunKernel2, animate = F)
#anim_save('IBM/plots/spread_map.gif') #saves last animation


#plot maps
plot_maps <- function(sp_grid){ #i is which tray
  tmp <- sp_grid #ex/157 "0.5.9.det.add.6"   
  #make df readable to ggplot
  # add to data a "ID" column for each feature
  tmp$id <- rownames(tmp@data)
  # create a data.frame from our spatial object
  tmp_df <- fortify(tmp, region = "id") %>% 
    #merge the "fortified" data with attribute data
    merge(., tmp@data, by = "id")
  
  #add in color id for the species so its the same every time
  Colors <- pal(spp)
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
               alpha=ifelse(tmp_centroids$state0=="C", 1, 1)) +
    coord_equal() +
    scale_fill_manual(values=Colors) +
    labs(title = paste(paste0('Tray', tmp_df$trayID[1], "-"), tmp_df$rand[1], tmp_df$dens[1], paste0('replicate', tmp_df$rep[1]), sep = '-'),
         subtitle = paste("SD =", tmp_df$SD, ', nplants = ', length(tmp), ', spacing = ', planting_dist, 'cm')) +
    scale_x_continuous(name='', breaks=xs, labels=1:length(xs), sec.axis = dup_axis()) +
    scale_y_continuous(name='', breaks=ys, labels=rev(LETTERS[1:length(ys)]), sec.axis = dup_axis())
  
  return(p1)
}

#for the cricut machine
samplegrid <- sample_community(c(2, 3), .1, 2) #species 2
plot_maps(samplegrid)
ggsave('GH_plots/cricket_diecut/cricket_2cm.pdf',units = 'in', width = 14.9, height = 12.13) #print with these dimensions to get actual size spacing. 
samplegrid <- sample_community(c(2, 3), .1, 1.75) #species 2
plot_maps(samplegrid)
ggsave('GH_plots/cricket_diecut/cricket_1.75cm.pdf',units = 'in', width = 14.9, height = 12.13) #print with these dimensions to get actual size spacing. 
samplegrid <- sample_community(c(2, 3), .1, 2.4) #species 2
plot_maps(samplegrid)
ggsave('GH_plots/cricket_diecut/cricket_2.4cm.pdf',units = 'in', width = 14.9, height = 12.13) #print with these dimensions to get actual size spacing. 
