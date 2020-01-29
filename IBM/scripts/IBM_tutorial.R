#IBM tutorial. Will keep the code in "IBM.R" as the source code for all the functions. Here, I'll put those to use, but the code will be cleaner (won't have the functions defined here).

#use community compositions from experimental design (take output from 'species_distributions.R')

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
#values approximated from Otten et al. (2003)
beta_curve <- function(x) .2*exp(-3*(log(x/11))^2)
alpha_curve <- function(x) .4*(1-.3)^x

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
sample_grid <- sample_community(nspecies = 1)

#quick plot
plot(sample_grid) + text(coordinates(sample_grid), cex=.5, col=as.numeric(sample_grid$spID))

#define agents dataframe----

#create dataframe that defines ID, coordinates, species ID, state, # total neighbors, neighbor ID vector. This dataframe will stay stationary through time. `IBM` will keep track of the state changes. 
agents <- make_agents_df(sample_grid)
head(agents)

#run epidemic----

#for every individual, if state==C, then C->C, C->S, or C->I; if state==S, then S->S, or S->I; if state==I, then always stays I. Each transition has a probability and the fate of the transition determined by a value drawn from a uniform distribution between 0 and 1. Will need to loop through every individual every time step. At each time step, record the state of every individual. Output should be a matrix of states, with rows equalling # individuals and cols equalling # time steps. 

testrun <- IBM(agents)
head(testrun)

#visualize----

#first plot summary of S and I
plotS_I(testrun)

#now plot spatial map of the spread (takes about 3 minutes)
plot_spread_map(sample_grid, testrun)
#anim_save('GH_plots/spread_map.gif') #saves last animation

