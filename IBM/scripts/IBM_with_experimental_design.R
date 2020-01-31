#IBM of my experimental design. "IBM_functions.R" is the source code for all the functions. Here, I'll put those to use on the same compostions that will be in the greenhouse.

#use community compositions from experimental design (take output from 'species_distributions.R')

rm(list=ls())

#load source code----
library(beepr)
source('IBM/scripts/IBM_functions.R')


#functions----
#for plotting spdf_list objects
plot_maps <- function(spatialdataframe){ #i is which tray
  tmp <- spatialdataframe   
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
  tmp_centroids <- data.frame(coordinates(tmp), state=tmp$state)
  
  #calculate axis labels
  xs <- unique(sort(round(tmp_centroids$X1, 4)))
  xs <- xs[seq(1, length(xs), by=2)] #get odds
  ys <- tmp_centroids$X2 %>% round(4) %>% sort %>% unique 
  
  #plot in ggplot!
  p1 <- ggplot(data = tmp_df, aes(x=long, y=lat)) +
    geom_polygon(aes(group = group, fill = spID))  +
    geom_path(aes(group = group), color = "white") +
    geom_point(data=tmp_centroids, aes(X1, X2), 
               alpha=ifelse(tmp_centroids$state=="C", 1, 0)) +
    coord_equal() +
    scale_fill_manual(values=Colors) +
    labs(title = paste(paste0('ID', tmp_df$trayID[1], "-"), tmp_df$rand[1], tmp_df$dens[1], paste0('replicate', tmp_df$rep[1]), sep = '-'),
         subtitle = paste("SD =", tmp_df$SD)) +
    scale_x_continuous(name='', breaks=xs, labels=1:length(xs), sec.axis = dup_axis()) +
    scale_y_continuous(name='', breaks=ys, labels=rev(LETTERS[1:length(ys)]), sec.axis = dup_axis())
  
  return(p1)
}

#bind the IBM output (state changes) to the treatment attribute data. 
bind_treatment_to_states <- function(spatialdataframe, IBM_output){
  
  #summarize state change matrix 
  sum_states <- function(matrix) cbind(sum(matrix %in% c("S", "C")), sum(matrix =="I"))
  #put into tall dataframe
  df1 <- data.frame(time=1:ncol(IBM_output), t(apply(IBM_output, 2, sum_states)))
  names(df1) <- c('time', 'S', 'I')
  #df1 <- pivot_longer(df1, cols=c('S', 'I'), names_to = 'state', values_to = 'count')
  
  #bind treatment to state changes
  treatment <- spatialdataframe@data[1,] %>% select(trayID, rand, dens, richness, rep, SD)
  row.names(treatment) <- NULL
  
  df2 <- cbind(treatment, df1)
  df2 <- df2 %>% mutate(percI = I/(S+I))
  
  return(df2)
}

#parameters----
#tray dimensions
width <- 9.5*2.54
length <- width 

#DESIGN (keep the names the same. simplifies downstream functions)
pinoc <- .1 #percent inoculated at start of experiemnt
s <- 6 #max number of species
spp <- c(paste0('sp_', rep(1:s))) #names of species
tfinal <- 25 #how many time steps
comp <- c(1, .5, .3, .2, 0, 0) #vector of relative "competencies"

#DISEASE 
#transmission curves
#values approximated from Otten et al. (2003)
beta_curve <- function(x) .2*exp(-3*(log(x/11))^2)
alpha_curve <- function(x) .4*(1-.3)^x
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
#read in list of spatial polygons dataframes for each tray
spdf_list <- readRDS('GH_output/species_distributions/spdf_list.RDS')

#plot one of them
plot_maps(spdf_list[[10]]) #can take ~30sec. if error, try again.

#run epidemic----

#for every individual, if state==C, then C->C, C->S, or C->I; if state==S, then S->S, or S->I; if state==I, then always stays I. Each transition has a probability and the fate of the transition determined by a value drawn from a uniform distribution between 0 and 1. Will need to loop through every individual every time step. At each time step, record the state of every individual. Output should be a matrix of states, with rows equalling # individuals and cols equalling # time steps. 

ptm <- proc.time()# Start the clock!
IBM_list_NN <- lapply(spdf_list, function(x) IBM(x, "NN") )
howlongIBM_NN <- proc.time() - ptm# Stop the clock

ptm <- proc.time()# Start the clock!
IBM_list_Kernel <- lapply(spdf_list, function(x) IBM(x, "Kernel", spatialdecay = .002) )
howlongIBM_Kernel <- proc.time() - ptm# Stop the clock
beep("mario")

#save IBM output
saveRDS(IBM_list_NN, 'IBM/outputs/IBM_list_NN.RDS')
saveRDS(IBM_list_Kernel, 'IBM/outputs/IBM_list_Kernel.RDS')

#read IBM output
IBM_list_NN <- readRDS('IBM/outputs/IBM_list_NN.RDS')
IBM_list_Kernel <- readRDS('IBM/outputs/IBM_list_Kernel.RDS')

#visualize single tray----

#### check out a single tray ####

#first plot summary of S and I
plotS_I(IBM_list_NN[[1]])
plotS_I(IBM_list_Kernel[[1]])

#now plot spatial map of the spread (takes about 3 minutes)
plot_spread_map(spdf_list[[1]], IBM_list_NN[[1]], animate = F)
plot_spread_map(spdf_list[[1]], IBM_list_Kernel[[1]], animate = F)
#anim_save('GH_plots/spread_map.gif') #saves last animation

#how do treatments affect D-D relationship----

#need to associate the state changes to treatments, which are contained in the spatial polygons dataframes. Do it for each community. Then you can see how the treatments affect disease across all trays.

trt_states_list <- lapply(1:length(spdf_list), function(i) bind_treatment_to_states(spdf_list[[i]], IBM_list_Kernel[[i]]))

#turn it into one dataframe
trt_states_df <- bind_rows(trt_states_list)
head(trt_states_df)

#plot changes in I (#inf and %inf) over time
filter(trt_states_df, SD==.5) %>% 
  ggplot(., aes(time, I, group=trayID)) +
  geom_line() +
  facet_grid(cols = vars(richness),
             rows = vars(rand, dens))

filter(trt_states_df, SD==.5) %>% 
  ggplot(., aes(time, percI, group=trayID)) +
  geom_line() +
  facet_grid(cols = vars(richness),
             rows = vars(rand, dens))

#plot I (#inf and %inf) at final time step
trt_states_df %>% filter(time==tfinal, SD==.5) %>% 
  ggplot(., aes(richness, I, group = rep)) +
  geom_line(alpha=.5, color='dodgerblue4') +
  facet_grid(cols = vars(dens),
             rows = vars(rand))

trt_states_df %>% filter(time==tfinal, SD==.5) %>% 
  ggplot(., aes(richness, percI, group = rep)) +
  geom_line(alpha=.5, color='dodgerblue4') +
  facet_grid(cols = vars(dens),
             rows = vars(rand))

#Community competency----
#community competency can be calculated before the simulation runs
#turn spdf_list into dataframes, bind the rows, then summarize by tray and include treatment, what it's community competency is. can then compare to disaese too.
binded_spdf <- suppressWarnings(bind_rows(lapply(1:length(spdf_list), function(x)spdf_list[[x]]@data))) 
CCsummary <- binded_spdf %>% 
  group_by(trayID, rand, dens, richness) %>% 
  summarise(nind = length(spID), CC = sum(comp), relCC =CC/nind)

#merge Com Comp summary with infection data
CCsummary2 <- trt_states_df %>% 
  filter(time==tfinal) %>% 
  select(trayID, I, percI, SD) %>% 
  left_join(CCsummary, ., by = "trayID") 
head(CCsummary2)

#plot CC vs infections
ggplot(CCsummary2, aes(relCC, percI, color=nind)) +
  geom_point() +
  scale_color_viridis_c()
ggplot(CCsummary2, aes(CC, I, color=richness)) +
  geom_point() 

#plot CC vs richness
ggplot(CCsummary2, aes(richness, relCC, color=percI, shape=as.factor(SD))) +
  geom_point() 
ggplot(CCsummary2, aes(richness, CC, color=percI, shape=as.factor(SD))) +
  geom_point() 










