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
#comp <- comp/comp['radish'] #comp standardized to radish

#parameters for model
delta <- 1/5 #1/average number of days inoc stays around
beta_ij_t <- make_beta_ij_t(comp)#matrix of amplitudes of the beta_ij 
alpha_i_t <- make_alpha_i_t(comp)#transmission from C -> I 

#simulate and return percent infected
sim_pI <- function(grid_df, spatialdecay=.001){
  simulation <- IBM(grid_df, Type = "Kernel", spatialdecay)
  sum(simulation[,ncol(simulation)] %in% "I")/sum(simulation[,ncol(simulation)] %in% c("I", "S"))
}

nreps <- 100
simulation_matrix <- matrix(NA, nrow = length(state_mat_list), ncol=nreps)
for(i in 1:length(state_mat_list)){
  simulation_matrix[i,] <- replicate(nreps, sim_pI(spdf_list[[i]])) 
}
#calculate deviation from 'observed'
apply(simulation_matrix, 1, mean) - comp #good for mustard and radish, underestimating disease for other species.


pdf('GH_plots/tinkering_stage/histogram_simvsobs_020720.pdf')
lapply(1:5, function(sp) 
  {hist(simulation_matrix[sp,], xlim=c(0,1), main='100 simulations vs. observed', sub = spp[sp], xlab='prop. infected')
  abline(v=comp[sp], col='red', lwd=2)
} )
dev.off()

#simulate single plot 
x <- 5 #which trayID
testrunKernel1 <- IBM(spdf_list[[x]], Type = "Kernel", spatialdecay = .001)

pdf('GH_plots/tinkering_stage/real_vs_simulated_radish.pdf')
#plot
plotS_I(state_mat_list[[x]]) #real
plotS_I(testrunKernel1) #simulated
plot_spread_map(spdf_list[[x]], state_mat_list[[x]], animate = F) #real
plot_spread_map(spdf_list[[x]], testrunKernel1, animate = F)#simulated
dev.off()

plot_spread_map(spdf_list[[x]], state_mat_list[[x]], animate = T) #real
anim_save('GH_plots/tinkering_stage/real_radish_animation.gif') #saves last animation

plot_spread_map(spdf_list[[x]], testrunKernel1, animate = T)#simulated
anim_save('GH_plots/tinkering_stage/simulated_radish_animation.gif') #saves last animation
