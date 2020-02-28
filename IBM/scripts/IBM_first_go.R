#source code for individual based model simulating rhizoctonia spread

#assume community filled in a hexagonal grid with r interplanting distance with s species, each with p competency. A portion of the communities will be inoculated and the epidemics will spread as a funciton of the community composition (species identity and density). States are C (challenged), S (susceptible), and I (infected). Transmission depends on number, identity, and distance from infected neighbors, as well as susceptibility of recipient host. 

#load packages----
rm(list=ls())
require(raster)
require(tidyverse)
require(ggplot2)
require(gganimate)


#parameters----
#tray dimensions
width <- 9.5*2.54
length <- width 

#DESIGN (keep the names the same. simplifies downstream functions)
pinoc <- .1 #percent inoculated at start of experiemnt
r <- 2 #interplanting distance
s <- 6 #number of species
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
  
  #standardizing to be between 0 and 1, relative to a distance of 0. 
  y0 <- C*exp(-sigma*0^2)
  y/y0
  
  #standardize value so that at 20 mm, value is 1. reason is that beta is derived from interplanting distances of 20mm and I want everything to be relative to that. 
  #y20 <- C*exp(-sigma*20^2)
  #y/y20
}


### transmission from C -> S ###
#inoculum decay rate invariant of species.
delta <- 1/5 #1/average number of days inoc stays around

### transmission from S -> I ###
#create a matrix of the beta_ij values. Assume that pairwise values will control the amplitude of beta(t). these will come from empirical data, but will just make them up for now.

#create matrix of amplitudes of the beta_ij 
make_beta_ij_t <- function(comp){
  beta_ij_amp <- matrix(comp, nrow=length(comp), ncol = length(comp)) #an array of pairwise beta_ij's across time. 
  beta_ij_amp <- beta_ij_amp * t(beta_ij_amp)
  
  #create an array of beta_ij, which varies by time. multiply the amplitude by beta(t)
  beta_ij_t <- array(beta_ij_amp, c(nrow(beta_ij_amp),nrow(beta_ij_amp),tfinal))
  beta_t <- beta_curve(1:tfinal) 
  beta_ij_t <- beta_ij_t * rep(beta_t, each=length(beta_ij_amp))
  rownames(beta_ij_t) <- colnames(beta_ij_t) <- spp #name dims
  
  return(beta_ij_t)
}

beta_ij_t <- make_beta_ij_t(comp)

### transmission from C -> I ###
#rate of infection from inoculum to plant. varies by species and time. 
make_alpha_i_t <- function(comp){
  a <- alpha_curve(1:tfinal)
  a_t <- t(matrix(comp, nrow=tfinal, ncol=length(comp), byrow = T) * a)
  rownames(a_t) <- spp #rows are species, cols are time.
  return(a_t)
}
alpha_i_t <- make_alpha_i_t(comp)


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
    spID=sample(c('sp_1', 'sp_6'), length(grid), T, prob = c(1.5,1)),
    state=sample(c(rep("C", ninoc), rep("S", length(grid)-ninoc)))
  ))

grid_df$state[1] <- NA #introduce NA to test function. might happen when seed doesn't germinate.

#quick plot
plot(grid_df) + text(coordinates(grid_df), cex=.5, col=as.numeric(grid_df$state))

#define agents dataframe----

#create dataframe that defines ID, coordinates, species ID, state. This dataframe will stay stationary through time. 
#FOR NEAREST NEIGHBOR: as long as you are directly adjacent to infected, interplanting distance doesn't matter. 
#FOR KERNEL MODEL: you need to know the distance of all indidividuals. Those neighbors will have betas that decay with distance.

get_NN <- function(grid_df){
  agents <- grid_df@data
  
  #pairwise adjacency matrix
  adj <- rgeos::gTouches(grid_df, byid = T) 
  
  return(list(agents=agents, adjacent=adj))
}

get_pairwise_dist <- function(grid_df, sigma){
  
  agents <- grid_df@data
  
  #translate a pairwise distance matrix (between centroids!!) into how probabilty of transmission decays with distance. make sure to convert cm to mm
  pair_dist <- as.matrix(dist(cbind(agents$x, agents$y), diag = T, upper=T)*10)
  pair_dist2 <- dist_decay(pair_dist, sigma=sigma) #this gets multiplied by its specific beta_ij. low value will decrease tranmsission probability. 
  diag(pair_dist2) <- 0 #set diagonals to zero. Self can't infect self. 
  
  return(list(agents=agents, pair_dist=pair_dist2))
}

#agents_NN <- get_NN(grid_df)
#agents_kernel <- get_pairwise_dist(grid_df)
#head(agents_NN)
#head(agents_kernel)

#rules of spread----

#for every individual, if state==C, then C->C, C->S, or C->I; if state==S, then S->S, or S->I; if state==I, then always stays I. Each transition has a probability and the fate of the transition determined by a value drawn from a uniform distribution between 0 and 1. Will need to loop through every individual every time step. At each time step, record the state of every individual. Output should be a matrix of states, with rows equalling # individuals and cols equalling # time steps. 

IBM <- function(grid_df, Type, spatialdecay=.001){

  #Type = type of transmission--"NN" or "Kernel". affects which function you run to get agents matrices
  if(Type=="NN"){
    agents <- get_NN(grid_df)
  }else{ #must be kernel
    agents <- get_pairwise_dist(grid_df, sigma=spatialdecay)
  }
  
  #split up the agents list
  agentsdf <- agents[[1]]
  agents_neigbors <- agents[[2]]
  
  #define matrix to fill in all the states. nrow=n_IDs, ncol=times
  states_matrix <- matrix(NA, nrow = nrow(agentsdf), ncol = tfinal)
  #initiate first time step
  states_matrix[,1] <- as.character(agentsdf$state)
  
  #inoculum decay
  C_to_S <- 1 - exp(-delta) 
  
  #define function for prob of inoculum to plant transmission
  C_to_Ia <- function(t){
    #alpha_i at time t
    rate <- alpha_i_t[as.character(agentsdf$spID),t]
    #rate of transmission C -> I via inoculum
    unname(1 - exp(-(rate))) 
  }
  
  #an array of pairwise transmission matrix
  #fill in the community grid with the beta_ij_t values
  #get pairwise betas for all of the individuals given their species identity.
  trans <- beta_ij_t[as.character(agentsdf$spID), as.character(agentsdf$spID),] 
  #account for distance decay or nearest neigbor
  trans2 <- trans*rep(agents_neigbors, times=dim(trans)[3])

  #change the species names to their IDs
  colnames(trans2) <- rownames(trans2) <- agentsdf$ID
  
  #define function for prob of plant to plant transmission
  C_or_S_to_Ib <- function(t, i){

    #infections happen at the rate: 
    #e^(-sum_j(beta_ij*(# infected_j)))
    #sum all of the beta_ij's for all of the infected IDs
    
    #are neighbors are infected?
    is_I <- states_matrix[,t] %in% "I" 
    #beta's of the infected neighbors
    beta_neighbors <- trans2[i,,t][is_I]
    
    #rate of transmission from S -> I 
    1 - exp(-1*sum(beta_neighbors)) #number infected
  }
  
  #which agents are not NA? NA usually because plant didn't germinate.
  germinated <- which(!is.na(agentsdf$state))
  
  for(t in 1:(tfinal-1)){ #at every time step
    #get alpha_i(t) (inoculum to plant transmission coef)
    alpha_i <- C_to_Ia(t)
    
    for(i in germinated){ #for every GERMINATED individual
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

testrunNN <- IBM(grid_df, Type = "NN")
testrunKernel <- IBM(grid_df, Type = "Kernel", spatialdecay = .002)
head(testrunNN)
head(testrunKernel)

#visualize----

#first plot summary of S and I
plotS_I <- function(IBM_output){
  sum_states <- function(data) cbind(sum(data %in% c("S", "C")), sum(data %in%"I"))
  df1 <- data.frame(time=1:ncol(IBM_output), t(apply(IBM_output, 2, sum_states)))
  names(df1) <- c('time', 'S', 'I')
  df1 <- pivot_longer(df1, cols=c('S', 'I'), names_to = 'state', values_to = 'count')
  p <- ggplot(df1, aes(time, count, color=state)) +
    geom_line()
  return(list(df1, p))
}

plotS_I(testrunNN)
plotS_I(testrunKernel)

#now plot spatial map of the spread

plot_spread_map <- function(spatialgrid_df, IBMoutput, animate=T){
  
  #for plotting
  pal <- function(spp) RColorBrewer::brewer.pal(n = length(spp), name = "RdYlBu")
  
  #dataframe matrix of states
  states_dfm <- data.frame(ID = spatialgrid_df$ID, IBMoutput) 
  
  # create a ggplot-readable df from our spatial object
  tmp <- fortify(spatialgrid_df, region = "ID") %>% 
    rename('ID'='id') %>% 
    #merge with attribute data and matrix of states
    merge(., spatialgrid_df@data, by= "ID") %>% 
    merge(., states_dfm, by = "ID") %>% 
    #remove state column (it's redundant)
    dplyr::select(-state) %>% 
    #pivot longer by states/time
    pivot_longer(cols = c(paste0("X", 1:ncol(IBMoutput))), 
                 names_to = 'time', values_to = 'state') %>% 
    #rename time to be numeric
    mutate(time = as.integer(gsub('X', '', time)) )
  
  #create a df of the centroids so infections can be plotted
  tmp_centroids <- tmp %>% 
    group_by(ID) %>% 
    top_n(tfinal, order) #takes 1 point per time step
  
  
  #add in color id for the species so its the same every time
  Colors <- pal(spp)
  names(Colors) <- spp
  
  #map it!
  staticplot <- ggplot(data=tmp, aes(long, lat, group = group)) +
    geom_polygon(aes(fill = spID)) +
    geom_path(color = "white") +
    coord_equal() +
    scale_fill_manual(values=Colors) +
    labs(x="", y="")
  
  if(animate==T){
    #takes about 3 minutes
    plot <- staticplot + 
      geom_point(data=tmp_centroids, aes(x, y), size = 3, alpha = ifelse(tmp_centroids$state=="I", .5, 0)) +
      transition_states(time) +
      ggtitle('time step {closest_state} of {tfinal}')
    
  }else{#plot static plot
    plot <- staticplot + 
      geom_point(data=tmp_centroids, aes(x, y), size=4, alpha = ifelse(tmp_centroids$state=="I", .05, 0)) #reduce alpha so it's easier to see
  }
  
  return(plot)
}

plot_spread_map(grid_df, testrunKernel, animate = T)
#anim_save('GH_plots/spread_map.gif') #saves last animation

