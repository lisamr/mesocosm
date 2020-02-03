#pipeline for processing 1010 tray inoculation trial data
#assumed that I've already created the designs previously, or just using some of the ones from the final experimental design. 

#load stuff----
rm(list=ls())
source('IBM/scripts/IBM_functions.R')
spdf_list <- readRDS('GH_output/species_distributions/spdf_list.RDS') #this is the final experimental design. can be any design though. 

#functions----
#for plotting spdf_list objects
plot_maps <- function(spatialdataframe, s=6){ #i is which tray
  tmp <- spatialdataframe   
  #make df readable to ggplot
  # add to data a "ID" column for each feature
  tmp$id <- rownames(tmp@data)
  # create a data.frame from our spatial object
  tmp_df <- fortify(tmp, region = "id") %>% 
    #merge the "fortified" data with attribute data
    merge(., tmp@data, by = "id")
  
  #add in color id for the species so its the same every time
  spp <- c(paste0('sp_', rep(1:s))) #names of species
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
    labs(title = paste(paste0('Tray', tmp_df$trayID[1], "-"), tmp_df$rand[1], tmp_df$dens[1], paste0('replicate', tmp_df$rep[1]), sep = '-'),
         subtitle = paste("SD =", tmp_df$SD)) +
    scale_x_continuous(name='', breaks=xs, labels=1:length(xs), sec.axis = dup_axis()) +
    scale_y_continuous(name='', breaks=ys, labels=rev(LETTERS[1:length(ys)]), sec.axis = dup_axis())
  
  return(p1)
}

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

#generate dataframe to fill out----
#ID numbers start in bottom left corner and move up L to R

#choose the trays you used
whichones <- c(1, 2)
#turn spdf_list into dataframe to record data
traydf <- suppressWarnings(bind_rows(lapply(whichones, function(i) spdf_list[[i]]@data))) 
traydf <- traydf %>% mutate(state_final = NA, day = NA)
write.csv(traydf, 'GH_output/tinkering_stage/sample_1010_datasheet.csv', row.names = F)

#make sense of data (part 1)----
data <- read.csv('GH_data/sample_1010_datasheet.csv')
head(data)

#put into the form of the IBM simulation. 
#generate a matrix of state changes. That can be plotted. Also bind state matrix to the original attribute data. 
tfinal <- max(data$day, na.rm = T) #number of days of monitoring
spp <- c(paste0('sp_', rep(1:6))) #names of species
state_mat <- matrix(data$state0, ncol=tfinal, nrow=nrow(data))

#could probably do something more elegant, but here's a 3-part if/else for-loop.
for(i in 1:nrow(state_mat)){
  if(is.na(data$day[i])){
    state_mat[i,] <- NA
  }else{
    if(data$state_final[i]=="S"){
      state_mat[i,] <- "S"
    }else{
      state_mat[i, data$day[i]:ncol(state_mat)] <- "I"
    }
  }
}

#split state matrix by tray so you can visualize individually
rows <- data %>%
  mutate(row=row_number()) %>%  
  group_by(trayID) %>% 
  summarise(firstrow=first(row), lastrow=last(row))

#choose a tray
whichtray <- rows %>% filter(trayID == 1)
state_mat_ind <- state_mat[c(whichtray$firstrow:whichtray$lastrow),]

#plot it
plotS_I(state_mat_ind)
plot_spread_map(spdf_list[[1]], state_mat_ind, animate = T)

#make sense of data (part 2)----

#this is not the fastest way to do this, but using the same functions as from the IBM simulation. 
bind_data <- function(i){
  whichtray <- rows %>% filter(trayID == i)
  state_mat_ind <- state_mat[c(whichtray$firstrow:whichtray$lastrow),]
  bind_treatment_to_states(spdf_list[[whichtray$trayID]], state_mat_ind)
}
trt_states_df <- bind_rows(lapply(whichones, bind_data))

#plot changes in I (#inf and %inf) over time
filter(trt_states_df, SD==.5) %>% 
  ggplot(., aes(time, I, group=trayID)) +
  geom_line(alpha=.5, color='dodgerblue4') +
  facet_grid(cols = vars(richness),
             rows = vars(rand, dens))

#plot I (#inf and %inf) at final time step
trt_states_df %>% filter(time==tfinal, SD==.5) %>% 
  ggplot(., aes(richness, I, group = rep)) +
  geom_line(alpha=.5, color='dodgerblue4') +
  facet_grid(cols = vars(dens),
             rows = vars(rand))




