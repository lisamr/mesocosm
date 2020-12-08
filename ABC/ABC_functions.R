#functions for ABC fitting

#convenience----
headmat <- function(mat, n=10) mat[1:n,1:n]
dens <- function(x, ...) plot(density(x), ...)


#creating community grid----

sample_community <- function(spp, perc_inoc=.1, planting_dist=2){
  #create size of tray and grid
  tray <- make_tray(24.13, 24.13)
  grid <- make_grid(r=planting_dist, tray) #2 cm interplanting distance
  
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


make_tray <- function(width, length){
  tray <- raster(ncol=1, nrow=1, crs="+proj=utm +zone=1 +datum=WGS84", xmn=0, xmx=width, ymn=0, ymx=length) 
  values(tray) <- 1
  tray <- rasterToPolygons(tray) #convert to sp poly
  return(tray)
}


make_grid <- function(r, tray){
  #tray=output from `make_tray()`
  hex_points <- spsample(tray, type = "hexagonal", cellsize = r, offset=c(0.1,0.1)) #setting square offset value ensures the same grid is drawn every time. 
  hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = r)
  return(hex_grid)
}




#run simulation----
K <- function(D, sigmak) exp(-sigmak*D^2) #distance kernel

f_sim1 <- function(beta, alpha, sigmak){
  
  #keep track of all individuals' states
  states_matrix <- matrix(NA, nrow=N, ncol = tfinal)
  states_matrix[,1] <- as.character(agents$state0)
  
  #vectors denoting specific states
  challenged <- as.numeric(agents$state0=='C')
  susceptible <- 1-challenged
  infected <- rep(0, N)
  
  #define distance kernel
  K1 <- K(d, sigmak)
  diag(K1) <- 0
  
  #Run the simulation
  for(t in 1:(tfinal-1)){ 
    
    #calculate probabilities of infections
    P_ci <- (1 - exp(-(alpha + (beta*K1)%*%infected)))*challenged
    P_si <- (1 - exp(- (beta*K1)%*%infected))*susceptible
    
    #generate new infections
    CI <- rbinom(N, 1, P_ci)
    SI <- rbinom(N, 1, P_si)
    new_infections <- CI + SI
    
    #update states
    challenged <- challenged - CI
    susceptible <- susceptible - SI
    infected <- infected + new_infections
    states_matrix[,t+1][challenged==1] <- 'C'
    states_matrix[,t+1][susceptible==1] <- 'S'
    states_matrix[,t+1][infected==1] <- 'I'
  }
  
  return(states_matrix)
}






#ABC algorithm----
#  Test if prior is non zero and within designated bounds
prior.non.zero<-function(par){
  H <- function(x) as.numeric(x>0)#Identity function: H(x)= 1 if x=T
  prod(sapply(1:length(par), function(a) H(par[a]-lower[a])* H(upper[a]-par[a])))
}

#SSE normalized by range of the constrained distribution
Norm.Eucl.dist<-function(p1,p2){ 
  sqrt(sum(((p1-p2)/(upper-lower))^2)) }

#  Covariance based on M neighbours
getSigmaNeighbours<-function(M, theta, Theta){
  #calculate top M parameter sets (closest to center of mass?)
  dist <- sapply(1:P, function(a) Norm.Eucl.dist(as.numeric(theta), as.numeric(Theta[a,])))
  temp <- data.frame(no=seq(1,P), dist)
  temp<-temp[order(temp$dist),]
  sigma <- cov(Theta[temp$no[1:(M+1)],])
  return(sigma)
}








#difference functions----
SSE <- function(S, Sstar){
  sqrt(sum(S - Sstar)^2)
}
chisq <- function(S, Sstar){ #pg 11 McKinley et al. 2009
  x <- ((S - Sstar)^2)/S
  x[is.na(x)] <- 0
  sum(x)
}


#summary statistics----

# 1. number of infecteds
nI <- function(x){
  sum(x=='I')
}


# 2. nearest infected neighbor

uniq.dist <- function(distmat) sort(unique(round(c(distmat), 8)))


nobs <- function(D, uniq.dist){
  n.obs <- matrix(NA, nrow = length(uniq.dist), ncol = ncol(D))
  for(i in 1:ncol(D)){
    inf <- which(D[,i]=="I")
    if(length(inf)>0){
      NinfN <- apply(d[inf, inf], 1, function(x) round(min(x[x!=0]), 8))
    }else NinfN <- 0
    n.obs[,i] <- sapply(1:length(uniq.dist), function(a) sum(NinfN==uniq.dist[a]))
  }
  return(n.obs)
}


# 3. O-ring
f_Oring <- function(ppp, Data, t, statei='I', R=NULL, Spar=.7){
  if(sum(Data[,t]==statei)>1 & sum(Data[,t]=="I")>1){
    marks(ppp) <- as.factor(Data[,t]) #update time
    K12 <- Kcross(ppp, statei, 'I') #calculate o-ring stat 
    g12 <- pcf(K12, method="a", spar=Spar)
    lambda2 <- summary(ppp)$marks['I',"intensity"]
    Oring <- eval.fv(lambda2*g12)
    if(is.null(R)){
      res <- Oring$pcf
    }else{
      res <- Oring$pcf[which.min(abs(Oring$r - R))] 
    }
  }else{
    res <- 0
  }
  return(res)
}


# 4. K function

#cumulative number of events across a range of distances. 
fK <- function(ppp, Data, t, state='I', R=NULL, C=c("border", "isotropic", "Ripley", "translate")){
  marks(ppp) <- as.factor(Data[,t])
  if(is.null(R)){
    tmp <- Kest(subset(ppp, marks == state), correction=C)
    ifelse(is.na(tmp[[3]]), 0, tmp[[3]])
  }else{
    tmp <- Kest(subset(ppp, marks == state), r = c(0,R), correction=C)
    ifelse(is.na(tmp[[3]][2]), 0, tmp[[3]][2])
    }
}

#marks(npp) <- as.factor(D[,5])
#Kest(subset(npp, marks == 'I')) %>% plot











#map plotting functions----
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
      geom_point(data=tmp_centroids, aes(x, y), shape = ifelse(is.na(tmp_centroids$state), 4, 16), size=4, alpha = ifelse(tmp_centroids$state %in% c("I", NA), .05, 0)) #reduce alpha so it's easier to see
  }
  
  return(plot)
}

plot_maps <- function(spatialdataframe, plotted_points=c("C"), point_cex = 1, point_shape = 16){ #i is which tray
  tmp <- spatialdataframe 
  pal <- function(spp) RColorBrewer::brewer.pal(n = length(spp), name = "RdYlBu")
  
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


