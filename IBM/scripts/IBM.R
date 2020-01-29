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
tfinal <- 20 #how many time steps

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
comp <- c(1, .5, .3, .2, 0, 0) #vector of relative "competencies"
names(comp) <- spp #named vector
beta_ij_amp <- matrix(comp, nrow=length(comp), ncol = length(comp)) #an array of pairwise beta_ij's across time. 
beta_ij_amp <- beta_ij_amp * t(beta_ij_amp)

#create an array of beta_ij, which varies by time. multiply the amplitude by beta(t)
beta_ij_t <- array(beta_ij_amp, c(nrow(beta_ij_amp),nrow(beta_ij_amp),tfinal))
beta_t <- beta_curve(1:tfinal) 
beta_ij_t <- beta_ij_t * rep(beta_t, each=length(beta_ij_amp))
rownames(beta_ij_t) <- colnames(beta_ij_t) <- spp #name dims

### transmission from C -> I ###
#rate of infection from inoculum to plant. varies by species and time. 
a <- alpha_curve(1:tfinal)
a_t <- t(matrix(comp, nrow=tfinal, ncol=length(comp), byrow = T) * a)
rownames(a_t) <- spp #rows are species, cols are time.


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
    spID=sample(c('sp_1', 'sp_2'), length(grid), T, prob = c(1.5,1)),
    state=sample(c(rep("C", ninoc), rep("S", length(grid)-ninoc)))
  ))

#quick plot
plot(grid_df) + text(coordinates(grid_df), cex=.5, col=as.numeric(grid_df$spID))

#define agents dataframe----

#create dataframe that defines ID, coordinates, species ID, state, # infected neighbors, # total neighbors, % infected neighbors, neighbor ID vector. This dataframe will stay stationary through time for the most part, except for state and # and % infected neighbors.

make_agents_df <- function(grid_df){
  agents <- grid_df@data
  
  #pairwise adjacency matrix
  adj <- rgeos::gTouches(grid_df, byid = T) 
  
  #NOTE: if you ever need to, you can create a pairwise distance decay matrix to account for declining transmission with distance. 
  
  
  #make the dataframe
  #which ones are neighbors
  agents$adj <- sapply(1:nrow(agents), function(i) {
    x <- which(adj[row.names(agents)[i],]==T)
    unname(x)
  }) 
  #count how many neighbors
  agents$n_adj <- unlist(lapply(1:nrow(agents), function(x) length(agents$adj[[x]])))
  
  return(agents)
}

agents <- make_agents_df(grid_df)
head(agents)

#rules of spread----

#for every individual, if state==C, then C->C, C->S, or C->I; if state==S, then S->S, or S->I; if state==I, then always stays I. Each transition has a probability and the fate of the transition determined by a value drawn from a uniform distribution between 0 and 1. Will need to loop through every individual every time step. At each time step, record the state of every individual. Output should be a matrix of states, with rows equalling # individuals and cols equalling # time steps. 

IBM <- function(agents){
  #agents is the data frame that gives info on the individuals
  
  #define matrix to fill in all the states. nrow=n_IDs, ncol=times
  states_matrix <- matrix(NA, nrow = nrow(agents), ncol = tfinal)
  #initiate first time step
  states_matrix[,1] <- as.character(agents$state)
  
  #inoculum decay
  C_to_S <- 1 - exp(-delta) 
  
  #define function for prob of inoculum to plant transmission
  C_to_Ia <- function(t){
    #alpha_i at time t
    rate <- a_t[agents$spID,t]
    #rate of transmission C -> I via inoculum
    unname(1 - exp(-(rate))) 
  }
  
  #an array of pairwise transmission matrix
  #fill in the community grid with the beta_ij_t values
  trans <- beta_ij_t[agents$spID, agents$spID,] #get pairwise betas for all of the individuals given their species identity.
  #change the species names to their IDs
  colnames(trans) <- rownames(trans) <- agents$ID
  
  #define function for prob of plant to plant transmission
  C_or_S_to_Ib <- function(t, i){
    #infections happen at the rate: 
    #e^(-sum_j(beta_ij*(# infected_j)/(total # neighbors)))
    #sum all of the beta_ij's for all of the infected IDs and divide by total # neighbors. 
    #beta's of the adjacent plants for a given ID
    beta_adj <- trans[,,t][agents$ID[i], agents$adj[[i]]]
    #are neighbors are infected?
    adj_is_I <- states_matrix[,t][agents$adj[[i]]]=="I" 
    #rate of transmission from S -> I 
    1 - exp(-1*sum(beta_adj*adj_is_I)) #number infected
    #1 - exp(-1*sum(beta_adj*adj_is_I)/agents$n_adj[i]) #proporiton infection
  }
  
  for(t in 1:(tfinal-1)){ #at every time step
    #get alpha_i(t) (inoculum to plant transmission coef)
    alpha_i <- C_to_Ia(t)

    for(i in 1:nrow(agents)){ #for every individual
      P <- runif(1, 0, 1)#draw random number
      
      if(states_matrix[i,t]=="C"){
        #calculate the probability of each transition
        CS <- C_to_S #inoc decay 
        CI <- C_to_Ia(t)[i] + C_or_S_to_Ib(t, i) #infection
        CC <- 1 - CS - CI #stasis
        
        #decide what the fate (state change) will be
        outcome_num <- 1 + (P < CS) + (P < (CS + CI))
        #C becomes... S if 3, I if 2, C if 1
        new_state <- switch(outcome_num, 'C', 'I', 'S')
        
      }else{
        if(states_matrix[i,t]=="S"){
          #calculate the probability of each transition
          SI <- C_or_S_to_Ib(t, i)
          SS <- 1 - SI
          
          #decide what the fate (state change) will be
          outcome_num <- 1 + (P < SI)
          #S becomes... I if 2, S if 1
          new_state <- switch(outcome_num, 'S', 'I')
          
        }else{
          #if state is I, stays I always
          new_state <- 'I'
        }
      }
      #fill states matrix for all individuals
      states_matrix[i,t+1] <- new_state
    }
  }
  return(states_matrix)
}

#run epidemic----
testrun <- IBM(agents)
head(testrun)

#visualize----

#first plot summary of S and I
sum_states <- function(data) cbind(sum(data %in% c("S", "C")), sum(data =="I"))
df1 <- data.frame(time=1:ncol(testrun), t(apply(testrun, 2, sum_states)))
names(df1) <- c('time', 'S', 'I')
df1 <- pivot_longer(df1, cols=c('S', 'I'), names_to = 'state', values_to = 'count')
ggplot(df1, aes(time, count, color=state)) +
  geom_line()

#now plot spatial map of the spread

#dataframe matrix of states
states_dfm <- data.frame(ID = grid_df$ID, testrun) 

# create a ggplot-readable df from our spatial object
tmp <- fortify(grid_df, region = "ID") %>% 
  rename('ID'='id') %>% 
  #merge with attribute data and matrix of states
  merge(., grid_df@data, by= "ID") %>% 
  merge(., states_dfm, by = "ID") %>% 
  #remove state column (it's redundant)
  dplyr::select(-state) %>% 
  #pivot longer by states/time
  pivot_longer(cols = c(paste0("X", 1:tfinal)), 
               names_to = 'time', values_to = 'state') %>% 
  #rename time to be numeric
  mutate(time = as.integer(gsub('X', '', time)) )

#create a df of the centroids so infections can be plotted
tmp_centroids <- tmp %>% 
  group_by(ID) %>% 
  top_n(tfinal, order) #takes 1 point per time step


#add in color id for the species so its the same every time
Colors <- pal
names(Colors) <- levels(tmp$spID)

#map it!
require(gganimate)
staticplot <- ggplot(data=tmp, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = spID)) +
  geom_path(color = "white") +
  coord_equal() +
  geom_point(data=tmp_centroids, aes(x, y), size = 3, alpha = ifelse(tmp_centroids$state=="I", .5, 0)) +
  scale_fill_manual(values=Colors)

#takes about 3 minutes
animplot <- staticplot + 
  transition_states(time) +
  ggtitle('time step {closest_state} of {tfinal}')
anim_save('GH_plots/spread_map.gif', animplot)

