#processing 1010 trays. single-species at 2 cm apart.

#load stuff----
rm(list=ls())
source('IBM/scripts/IBM_functions.R')

#parameters----
#tray dimensions
width <- 9.5*2.54
length <- width 

#DESIGN (keep the names the same. simplifies downstream functions)
pinoc <- .1 #percent inoculated at start of experiemnt
spp <- c('kale', 'PC', 'arugula', 'mustard', 'radish') #names of species

#create trays----
spdf_list <- list(NULL)
for(i in 1:length(spp)){
  spdf_list[[i]] <- sample_community(which_spp = i, perc_inoc = pinoc, planting_dist = 2)
  spdf_list[[i]]$trayID <- i
}

#generate dataframe to fill out----

#ID numbers start in bottom left corner and move up L to R
#trays missing last row. fill in with NAs.

#choose the trays you used
whichones <- c(1:length(spdf_list))
#turn spdf_list into dataframe to record data
traydf <- suppressWarnings(bind_rows(lapply(whichones, function(i) spdf_list[[i]]@data))) 
traydf <- traydf %>% mutate(day = NA) %>% select(trayID, everything())
write.csv(traydf, 'GH_output/tinkering_stage/1010trial020720_datasheet.csv', row.names = F)

#get state matrix----

#after doing your experiment, enter data and save.
data <- readxl::read_xls('GH_data/1010trial020720_datasheet.xls')
is.na(data) <- data=="NA" #counterintuitivley assigns character NAs as true NAs.
data <- data %>% mutate(day=as.numeric(day))
head(data)

#put into the form of the IBM simulation. 
#generate a matrix of state changes. That can be plotted. Also bind state matrix to the original attribute data. 
tfinal <- max(data$day, na.rm = T) #number of days of monitoring
spp <- unique(data$spID) #names of species

#fill in the state matrix
state_mat <- matrix("S", ncol=tfinal, nrow=nrow(data))
for(i in 2:(tfinal)){
  state_mat[,i] <- state_mat[,i-1]#current time same as last time
  I <- which(data$state_final=='I' & data$day==i)#find infecteds
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

#update the spatial dataframes with any changes in NAs and states (if inoculations moved)
for(i in 1:length(spdf_list)){
  spdf_list[[i]]$state0 <- state_mat_list[[i]][,1]
}

#plot a tray
plotS_I(state_mat_list[[5]])
plot_spread_map(spdf_list[[5]], state_mat_list[[5]], animate = F)

#match data to model----

#calculate proportion of secondary infections. "relative competencies"
comp <- rep(NA, length=length(state_mat_list))
names(comp) <- spp
for(i in 1:length(state_mat_list)){
  states <- state_mat_list[[i]]
  S0 <- sum(states[,1]=="S", na.rm = T)
  Itfinal <- sum(states[,tfinal]=="I", na.rm = T)
  comp[i] <- Itfinal/S0 #proportion infected from secondary transmission
}
comp #use this to feed your model. 

#parameters for model
delta <- 1/5 #1/average number of days inoc stays around
beta_ij_t <- make_beta_ij_t(comp)#matrix of amplitudes of the beta_ij 
alpha_i_t <- make_alpha_i_t(comp)#transmission from C -> I 

#simulate
testrunKernel1 <- IBM(samplegrid1, Type = "Kernel", spatialdecay = .002)


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




