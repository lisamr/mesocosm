#code to generate spreadsheet to enter data for big experiment. Also to process the data.


#LOAD SOURCE CODE----

rm(list=ls())
source('IBM/scripts/IBM_functions.R')

#functions to be added to source code...
test_track_individuals <- function(spatialdataframe, IBM_output){
  
  #get final state
  spatialdataframe$state_tf <- IBM_output[,ncol(IBM_output)]
  
  #summarize things about the tray: density, density of species, avg host competency
  spdf <- as.data.frame(spatialdataframe)
  trayinfo1 <- c('density'=nrow(spdf), 'avgCC'=mean(spdf$comp))
  trayinfo2 <- sapply(1:length(spp), function(i) sum(spdf$spID==spp[i]))
  names(trayinfo2) <- spp
  trayinfo <- c(trayinfo1, trayinfo2)
  trayinfo_mat <- matrix(rep(trayinfo, nrow(spdf)), ncol = length(trayinfo), byrow = T)
  colnames(trayinfo_mat) <- names(trayinfo)
  
  #bind trayinfo to individual info
  output <- cbind(as.data.frame(spatialdataframe), trayinfo_mat) 
  return(output)
}


#IMPORT DESIGN AND SPATIAL DATAFRAME----

design <- read.csv('GH_output/real_experiment/big_experiment_design_03302020.csv')
spdf_list <- readRDS('GH_output/real_experiment/big_experiment_spdf_list_03302020.RDS')


#GENERATE DATAFRAME----

#ID numbers start in bottom left corner and move up L to R
#trays missing last row. fill in with NAs.

add_axes <- function(spatial_object){
  #convert spatial object into dataframe
  spdf <- spatial_object@data
  
  #create axis columns that appear on maps
  #get x and y coordinates
  xs <- unique(sort(round(spdf$x, 4)))
  ys <- spdf$y %>% round(4) %>% sort %>% unique
  
  #translate coordinates into letters or integers
  horiz.ax <- rep(1:ceiling(length(xs)/2), each=2)
  #have to correct for rows that have odd numbers.
  horiz.ax <- if(length(xs)/2 != ceiling(length(xs)/2)){
    horiz.ax[-length(horiz.ax)]
  }else{
    horiz.ax
  }
  horiz.ax <- horiz.ax + rep(c(0,.5), length.out=length(horiz.ax))
  vert.ax <- LETTERS[rev(1:length(ys))]
  names(horiz.ax) <- xs
  names(vert.ax) <- ys
  
  #add axis columns
  spdf <- spdf %>% 
    mutate(x = round(x, 4),
           y = round(y, 4),
           horiz = recode(x, !!!horiz.ax),
           vert = recode(y, !!!vert.ax))
  
  return(spdf)
}

#choose the trays you used. Didn't plant trays 10, 24, 25, 29, 30, 34, 35, 37, 38, 39, 40
design

omit <- c(107:116)
whichones <- c(1:length(spdf_list))[-omit]


#turn grid_list into dataframe to record data
traydf <- suppressWarnings(bind_rows(lapply(whichones, function(i) add_axes(spdf_list[[i]])))) 

traydf <- traydf %>% mutate(day_infected = NA) %>% select(trayID, horiz, vert, everything(), -comp)

write.csv(traydf, 'GH_output/real_experiment/big_experiment_03302020_datasheet.csv', row.names = F)

#GET STATE MATRIX----

#after doing your experiment, enter data and save.
data <- read_csv('GH_data/real_experiment/big_experiment_03302020_datasheet.csv')
print(data, width=Inf)

#determine final state of each individual. 
data <- data %>%
  mutate(state_final = case_when(
    is.na(state0) ~ NA_character_, #can't use NA. must designate which NA to use. (NA_real_, NA_character_, NA_integer_, NA_complex_)
    is.na(day_infected) ~ "S",
    day_infected>0 ~ "I"
  )) 

#put into the form of the IBM simulation. 
#generate a matrix of state changes. That can be plotted. Also bind state matrix to the original attribute data. 
tfinal <- max(data$day_infected, na.rm = T) #number of days of monitoring
spp <- unique(data$spID) #names of species

#fill in the state matrix
state_mat <- matrix("S", ncol=tfinal, nrow=nrow(data))
for(i in 2:(tfinal)){
  state_mat[,i] <- state_mat[,i-1]#current time same as last time
  I <- which(data$state_final=='I' & data$day_infected==i)#find infecteds
  state_mat[I,i] <- "I"#update state if infected
}
#replace first column with challenged status and rows as NA
state_mat[,1] <- data$state0
state_mat[is.na(data$state_final),] <- NA

#split state matrix by tray
rows <- data %>%
  mutate(row=row_number()) %>%  
  group_by(trayID) %>% 
  summarise(firstrow=first(row), lastrow=last(row))

state_mat_list <- list(NULL)
for(i in rows$trayID){
  whichtray <- rows %>% filter(trayID == i)
  state_mat_list[[i]] <- state_mat[c(whichtray$firstrow:whichtray$lastrow),]
}

#update the spatial dataframes with any changes in NAs and states (if tray wasn't planted, inoculations moved, species identity changed, or no germination)
for(i in 1:length(state_mat_list)){
  spdf_list[[i]]$state0 <- state_mat_list[[i]][,1]
}

#plot a tray
plotS_I(state_mat_list[[1]])
plot_spread_map(spdf_list[[1]], state_mat_list[[1]], animate = F)

#MAKE SENSE OF DATA (PART 1)----

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

filter(trt_states_df, SD==.5) %>% 
  ggplot(., aes(time, percI, group=trayID)) +
  geom_line(alpha=.5, color='dodgerblue4') +
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

#QUESTION 3 ANALYSIS----
library(brms)
library(tidybayes)
dens <- function(x, ...) plot(density(x), ...)
inv_logit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))
Scale <- function(x){
  (x-mean(x))/(2*sd(x))
}

#do a logistic regression to predict the probability of a secondary infection. do on invididuals. bernoulli trial. 
ind_trials_list <- lapply(1:length(state_mat_list), function(x) test_track_individuals(spdf_list[[x]], state_mat_list[[x]]))
ind_trials <- bind_rows(ind_trials_list)
ind_trials$infected <- ifelse(ind_trials$state_tf=="I", 1, 0)
head(ind_trials)
tail(ind_trials)




