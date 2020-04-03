#testing out how I might fit model to data. starting with RMSE. 

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

#get state matrix for observed data----

#after doing your experiment, enter data and save.
data <- readxl::read_xls('GH_data/1010trial020720_datasheet.xls')
is.na(data) <- data=="NA" #counterintuitivley assigns character NAs as true NAs.
data <- data %>% mutate(day=as.numeric(day))
head(data)

#put into the form of the IBM simulation. 
#generate a matrix of state changes. That can be plotted. Also bind state matrix to the original attribute data. 
tmonitored <- max(data$day, na.rm = T) #number of days plants monitored
tfinal <- tmonitored - 2 #for simulation. offset time by 2 days since I inoculate on day3, instead of day1
spp <- unique(data$spID) #names of species

#fill in the state matrix
state_mat <- matrix("S", ncol=tmonitored, nrow=nrow(data))
for(i in 2:(tmonitored)){
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
plotS_I(state_mat_list[[3]])
plot_spread_map(spdf_list[[3]], state_mat_list[[3]], animate = F)

#simulate data----
#PI <- c(0.3931624, 0.1617647, 0.3161765, 0.5895522, 0.6335878)#Percent infected for the 5 trays

#start with just 1 tray 
tray <- 4
obs_data <- state_mat_list[[tray]]
obs_tray <- spdf_list[[tray]]

####################################
#COARSE GRID APPROXIMATION
####################################

#parameters for model
#make a grid of parameter values. length of ab_combos is how many parameter vectors you'll be testing 
daysmonitored <- (2:tfinal)
daysmonitored <- daysmonitored[daysmonitored %% 2 != 0] 

betas <- seq(0,1, by=.1)
alphas <- seq(0,1, by=.1)
ab_combos <- expand.grid(beta = betas, alpha = alphas) #parameters tested
beta_mat <- matrix(rep(ab_combos[,1], each = length(spp)), ncol = 5, byrow = T)
alpha_mat <- matrix(rep(ab_combos[,2], each = length(spp)), ncol = 5, byrow = T)
delta <- 1/5 #guessed 1/average number of days inoc stays around
sp_decay <- .001 #spatial decay parameter
nsims <- 2 #how many simulations performed per parameter value 

#assess fit (with least squares right now)
#make sure comparisons are made only on the days you monitor
f_SSE <- function(pred, obs){
  I_pred <- apply(pred, 2, function(x) sum(x %in% "I"))
  I_obs <- apply(obs, 2, function(x) sum(x %in% "I"))
  I_obs <- I_obs[-c(1:2)] #remove first 2 days so obs and pred can be directly compared relative to days 'inoculated'
  RMSE <- sqrt(sum((I_pred[daysmonitored] - I_obs[daysmonitored])^2) / length(daysmonitored)) #sum of sq errors
  return(RMSE)
}

#calculate RMSE for all combos of parameters
error_coarse <- matrix(NA, ncol = nsims, nrow = nrow(beta_mat))
for(i in 1:nrow(beta_mat)){
  beta_ij_t <- make_beta_ij_t(beta_mat[i,])
  alpha_i_t <- make_alpha_i_t(alpha_mat[i,])

  #simulate model with above parameters many times
  sim_test <- replicate(nsims, IBM(grid_df = obs_tray, Type = 'Kernel', spatialdecay = sp_decay))
  #dim(sim_test) #array of output matrices
  #[1] 168  22  50
  error_coarse[i,] <- sapply(1:nsims, function(x) f_SSE(sim_test[,,x], obs_data))

  }

meanRMSE <- apply(error_coarse, 1, mean)
medianRMSE <- apply(error_coarse, 1, median)
sdRMSE <- apply(error_coarse, 1, sd)
#put error into data frame with parameters
dat_coarse <- data.frame(ab_combos, error_coarse, meanRMSE, medianRMSE, sdRMSE)
dat_coarse2 <- dat_coarse %>% 
  pivot_longer(
    starts_with("X"),
    names_prefix = 'X',
    names_to = 'sim_rep',
    values_to = 'RMSE')
  
#viz 1, surface plot
dat_coarse2 %>% filter(sim_rep == 1) %>% 
  ggplot(., aes(alpha, beta, z = medianRMSE)) +
  geom_raster(aes(fill = medianRMSE)) + 
  geom_contour(color = 'white', bins = 20) +
  scale_fill_viridis_c()

#viz 2, scatter plot  
ggplot(dat_coarse2, aes(beta, RMSE, color = alpha, group = alpha)) +
  geom_point() +
  geom_smooth() +
  scale_color_viridis_c() +
  facet_wrap(~alpha)


#radish: beta=.6 alpha=.75
#arugula: beta=.3 alpha=.9
#mustard: beta=.52 alpha=.52

####################################
#FINE GRID APPROXIMATION
####################################
#parameters for model
betas <- seq(.4,.8, by=.05)
alphas <- seq(.25,1, by=.05)
ab_combos <- expand.grid(beta = betas, alpha = alphas) #parameters tested
beta_mat <- matrix(rep(ab_combos[,1], each = length(spp)), ncol = 5, byrow = T)
alpha_mat <- matrix(rep(ab_combos[,2], each = length(spp)), ncol = 5, byrow = T)
delta <- 1/5 #guessed 1/average number of days inoc stays around
sp_decay <- .001 #spatial decay parameter
nsims <- 20 #how many simulations performed per parameter value 

#calculate SSE for all combos of parameters
error_fine <- matrix(NA, ncol = nsims, nrow = nrow(beta_mat))
for(i in 1:nrow(beta_mat)){
  beta_ij_t <- make_beta_ij_t(beta_mat[i,])
  alpha_i_t <- make_alpha_i_t(alpha_mat[i,])
  
  #simulate model with above parameters many times
  sim_test <- replicate(nsims, IBM(grid_df = obs_tray, Type = 'Kernel', spatialdecay = sp_decay))
  #dim(sim_test) #array of output matrices
  #[1] 168  22  50
  error_fine[i,] <- sapply(1:nsims, function(x) f_SSE(sim_test[,,x], obs_data))
}

meanSSE <- apply(error_fine, 1, mean)
medianSSE <- apply(error_fine, 1, median)
sdSSE <- apply(error_fine, 1, sd)
#put error into data frame with parameters
dat_fine <- data.frame(ab_combos, error_fine, meanSSE, medianSSE, sdSSE)
dat_fine2 <- dat_fine %>% 
  pivot_longer(
    starts_with("X"),
    names_prefix = 'X',
    names_to = 'sim_rep',
    values_to = 'SSE')
  

#viz 1, surface plot
dat_fine2 %>% filter(sim_rep==1) %>% 
  ggplot(., aes(alpha, beta, z = medianSSE)) +
  geom_tile(aes(fill = medianSSE)) + 
  #geom_contour(color = 'white', bins = 20) +
  scale_fill_viridis_c()

#viz 2, scatter plot  
ggplot(dat_fine2, aes(beta, medianSSE, color = alpha, group = alpha)) +
  geom_point() +
  geom_line() +
  scale_color_viridis_c()
ggplot(dat_fine2, aes(beta, medianSSE, color = alpha, group = alpha)) +
  geom_pointrange(aes(ymin=medianSSE-sdSSE, ymax=medianSSE+sdSSE)) +
  scale_color_viridis_c() 
ggplot(dat_fine2, aes(beta, medianSSE)) +
  geom_pointrange(aes(ymin=medianSSE-sdSSE, ymax=medianSSE+sdSSE)) +
  facet_wrap(~alpha)

#best fitting model...
bestfit <- dat_fine2[which.min(dat_fine2$meanSSE),]


#predict from best fitting model
beta_ij_t <- make_beta_ij_t(rep(bestfit$beta, 5))
alpha_i_t <- make_alpha_i_t(rep(bestfit$alpha, 5))

nsims=30
predicted <- replicate(nsims, IBM(grid_df = obs_tray, Type = 'Kernel', spatialdecay = sp_decay))
nI <- function(x) apply(x, 2, function(x) sum(x %in% "I"))
predI <- sapply(1:nsims, function(x) nI(predicted[,,x]))
obsI <- nI(obs_data[,-c(1:2)])
plot(NULL, ylim=c(0, 150), xlim=c(1, length(obsI)), xlab = 'days after inoculation', ylab = '# infected')
lapply(1:nsims, function(x) lines(predI[,x], col=scales::alpha('grey20', alpha=.3), lwd=1.5))
lines(obsI, lwd=2, col='firebrick')



