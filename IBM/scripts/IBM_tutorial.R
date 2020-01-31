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
sample_grid <- sample_community(nspecies = 1, planting_dist = 1)

#quick plot
plot(sample_grid) + text(coordinates(sample_grid), cex=.5, col=as.numeric(sample_grid$spID))


#run epidemic----

#for every individual, if state==C, then C->C, C->S, or C->I; if state==S, then S->S, or S->I; if state==I, then always stays I. Each transition has a probability and the fate of the transition determined by a value drawn from a uniform distribution between 0 and 1. Will need to loop through every individual every time step. At each time step, record the state of every individual. Output should be a matrix of states, with rows equalling # individuals and cols equalling # time steps. 

#try out transmission with nearest neighbor only infections and ones with distance decay (kernel) on two communities with different densities

#communities
samplegrid1 <- sample_community(1, 1) #planting dist=1cm
samplegrid2 <- sample_community(1, 2) #planting dist=2cm

#simulate!
testrunNN1 <- IBM(samplegrid1, Type = "NN")
testrunNN2 <- IBM(samplegrid2, Type = "NN")
testrunKernel1 <- IBM(samplegrid1, Type = "Kernel", spatialdecay = .002)
testrunKernel2 <- IBM(samplegrid2, Type = "Kernel", spatialdecay = .002)

#view
head(testrunNN1); head(testrunKernel1)

#visualize----

#first plot summary of S and I
p1 <- plotS_I(testrunNN1)[[2]]
p2 <- plotS_I(testrunNN2)[[2]]
p3 <- plotS_I(testrunKernel1)[[2]]
p4 <- plotS_I(testrunKernel2)[[2]]
cowplot::plot_grid(p1, p2, p3, p4, labels = c('NN1', 'NN2', 'Kernel1', 'Kernel2'))


#now plot spatial map of the spread (animation is ~3 minutes)
plot_spread_map(samplegrid2, testrunNN2, animate = F)
plot_spread_map(samplegrid2, testrunKernel2, animate = F)
#anim_save('IBM/plots/spread_map.gif') #saves last animation

