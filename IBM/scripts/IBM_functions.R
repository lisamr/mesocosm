#source code for individual based model simulating rhizoctonia spread

#assume community filled in a hexagonal grid with r interplanting distance with s species, each with p competency. A portion of the communities will be inoculated and the epidemics will spread as a funciton of the community composition (species identity and density). States are C (challenged), S (susceptible), and I (infected). Transmission depends on number, identity, and distance from infected neighbors, as well as susceptibility of recipient host. 

#load packages----
require(raster)
require(tidyverse)
require(ggplot2)
require(gganimate)

#transmission coefficients----

#transmission curves for Rhizoctonia solani
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

### transmission from C -> I ###
#rate of infection from inoculum to plant. varies by species and time. 
make_alpha_i_t <- function(comp){
  a <- alpha_curve(1:tfinal)
  a_t <- t(matrix(comp, nrow=tfinal, ncol=length(comp), byrow = T) * a)
  rownames(a_t) <- spp #rows are species, cols are time.
  return(a_t)
}

#for plotting
pal <- function(spp) RColorBrewer::brewer.pal(n = length(spp), name = "RdYlBu")
theme_set(theme_bw())

#create grid----
#First create hexagonal grid to figure out constraints.
#planting area for 1010 tray is 9.5in^2. 
make_tray <- function(width, length){
  tray <- raster(ncol=1, nrow=1, crs="+proj=utm +zone=1 +datum=WGS84", xmn=0, xmx=width, ymn=0, ymx=length) 
  values(tray) <- 1
  tray <- rasterToPolygons(tray) #convert to sp poly
  return(tray)
}

#create hexagon grid
make_grid <- function(r, tray){
  #tray=output from `make_tray()`
  hex_points <- spsample(tray, type = "hexagonal", cellsize = r, offset=c(0.1,0.1)) #setting square offset value ensures the same grid is drawn every time. 
  hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = r)
  return(hex_grid)
}



#create communities----

#spatial dataframe will be more streamlined. will come from the design from your data. 
sample_community <- function(which_spp, perc_inoc=.1, planting_dist=2){
  #create size of tray and grid
  tray <- make_tray(24.13, 24.13)
  grid <- make_grid(r=planting_dist, tray) #2 cm interplanting distance
  
  #how many species
  spp <- spp[which_spp]
  
  #fill community
  ninoc <- round(perc_inoc*length(grid))
  grid_df <- SpatialPolygonsDataFrame(
    grid, data.frame(
      ID=row.names(grid),
      x=coordinates(grid)[,1],
      y=coordinates(grid)[,2],
      spID=sample(spp, length(grid), T),
      state0=sample(c(rep("C", ninoc), rep("S", length(grid)-ninoc)))
    ))
  
  return(grid_df)
}

#sample_grid <- sample_community(3)
#quick plot
#plot(sample_grid) + text(coordinates(sample_grid), cex=.5, col=as.numeric(sample_grid$spID))

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
  
  
  #define function for prob of inoculum to plant transmission
  alpha_t <- function(t){
    #alpha_i at time t
    alpha_i_t[as.character(agentsdf$spID),t]
  }
  
  #an array of pairwise transmission matrix
  #fill in the community grid with the beta_ij_t values
  #get pairwise probabilities of infection for all of the individuals given their species identity.
  pij_t <- beta_ij_t[as.character(agentsdf$spID), as.character(agentsdf$spID),] 
  #account for distance decay or nearest neigbor
  beta_t <- pij_t*rep(agents_neigbors, times=dim(pij_t)[3])
  
  #change the species names to their IDs
  colnames(beta_t) <- rownames(beta_t) <- agentsdf$ID
  
  #define function for prob of plant to plant transmission
  beta_I <- function(t, i){
    
    #(sum_j(beta_ij*(# infected_j)))
    #sum all of the beta_ij's for all of the infected IDs
    
    #are neighbors are infected?
    is_I <- states_matrix[,t] %in% "I" 
    #beta's of the infected neighbors
    beta_infecteds <- beta_t[i,,t][is_I]
    
    #summed beta for all infected neighbors
    return(sum(beta_infecteds))
  }
  
  #which agents are not NA? NA usually because plant didn't germinate.
  germinated <- which(!is.na(agentsdf$state))
  
  for(t in 1:(tfinal-1)){ #at every time step
    #get alpha_i(t) (inoculum to plant transmission coef)
    alpha_i <- alpha_t(t)
    for(i in germinated){ #for every GERMINATED individual
      if(states_matrix[i,t]=="C"){
        #calculate the probability of infection
        CI <- 1 - exp(-(alpha_i[i] + beta_I(t, i)))
        #decide if individual gets infected
        new_state <- ifelse(rbinom(1, 1, prob = CI)==1, 'I', 'C')
      }else{
        if(states_matrix[i,t]=="S"){
          #calculate the probability of infection
          SI <- 1 - exp(-beta_I(t, i))
          #decide if individual gets infected
          new_state <- ifelse(rbinom(1, 1, prob = SI)==1, 'I', 'S')
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
#IBM(sample_grid, Type="NN") #OR...
#IBM(sample_grid, Type="Kernel")


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

#now plot spatial map of the spread
plot_spread_map <- function(spatialgrid_df, IBMoutput, animate=T, alpha_intensity=.05){
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
    dplyr::select(-state0) %>% 
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
    #takes about a minute
    plot <- staticplot + 
      geom_point(data=tmp_centroids, aes(x, y), shape = ifelse(is.na(tmp_centroids$state), 4, 16),size = 3, alpha = ifelse(tmp_centroids$state %in% c("I", NA), .9, 0)) +
      transition_states(time, .5, 1) +
      ggtitle('time step {closest_state} of {tfinal}')
    plot <- animate(plot, nframes=tfinal, fps=5, duration=10)
    
  }else{#plot static plot
    plot <- staticplot + 
      geom_point(data=tmp_centroids, aes(x, y), shape = ifelse(is.na(tmp_centroids$state), 4, 16), size=4, alpha = ifelse(tmp_centroids$state %in% c("I", NA), alpha_intensity, 0)) #reduce alpha so it's easier to see
  }
  
  return(plot)
}

#plotting functions after simulation----
#for plotting spdf_list objects
plot_maps <- function(spatialdataframe, plotted_points=c("C"), point_cex = 1, point_shape = 16){ #i is which tray
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
  #names(Colors) <- levels(tmp_df$spID)
  names(Colors) <- spp
  
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
               cex = point_cex,
               shape = point_shape,
               alpha=ifelse(tmp_centroids$state0 %in% plotted_points, 1, 0)) +
    coord_equal() +
    scale_fill_manual(values=Colors) +
    labs(title = paste(paste0('Tray', tmp_df$trayID[1], "-"), tmp_df$rand[1], tmp_df$dens[1], paste0('replicate', tmp_df$rep[1]), sep = '-'),
         subtitle = paste("SD =", tmp_df$SD, ', nplants = ', length(tmp), ', spacing = ', planting_dist, 'cm')) +
    scale_x_continuous(name='', breaks=xs, labels=1:length(xs), sec.axis = dup_axis()) +
    scale_y_continuous(name='', breaks=ys, labels=rev(LETTERS[1:length(ys)]), sec.axis = dup_axis())
  
  return(p1)
}

#bind the IBM output (state changes) to the treatment attribute data. 
bind_treatment_to_states <- function(spatialdataframe, IBM_output){
  
  #summarize state change matrix 
  sum_states <- function(matrix) cbind(sum(matrix %in% c("S", "C")), sum(matrix %in% "I"))
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


#for statistics----
#need to keep track of individuals in order to do a bernoulli regression. going to assume we're keeping 6 species.
track_individuals <- function(spatialdataframe, IBM_output){
  
  #get final state
  spatialdataframe$state_tf <- IBM_output[,ncol(IBM_output)]
  
  #summarize things about the tray: density, density of species, avg host competency
  trayinfo <- spatialdataframe %>% 
    as.data.frame() %>% 
    summarise(density = length(spID),
              nsp1 = sum(spID=="sp_1"),
              nsp2 = sum(spID=="sp_2"),
              nsp3 = sum(spID=="sp_3"),
              nsp4 = sum(spID=="sp_4"),
              nsp5 = sum(spID=="sp_5"),
              nsp6 = sum(spID=="sp_6"),
              avgCC = sum(comp)/density)
  
  #bind trayinfo to individual info
  output <- cbind(as.data.frame(spatialdataframe), trayinfo) 
  return(output)
}


