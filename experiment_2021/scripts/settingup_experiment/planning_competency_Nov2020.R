#Planning out host competency trial in Nov 2020

#QUESTIONS
#1. Should I plant more trays at lower capacity, or fewer trays at full capacity?  
#2. How many trays do I need for the low competency hosts vs. high competency hosts? 
#3. Are results the same if I quantify final secondary infection prevalence vs. rate of spread?

#APPROACH
#1. simulate two sets of data and crunch the numbers. Compare the distributions of mean infection rates. 
#2. Vary up the number of trays for each of the species and see how it changes the results (sensitivity analysis).
#3. Crunch the numbers for two different analyses--a)binomial model, b) non-linear curve fitting?
#4. fit the curves with a SIR model and estimate beta and R0. 

#ANSWERS
#1. PLANT AT FULL CAPACITY. SMALLER TRAYS LEADS TO CONSISTENTLY LOWER COMPETENCY VALUES.
#2. VARIATION AFTER APPROXIMATELY 4 AND 6 SEEMS TO LEVEL OFF for low and high, respectively. 5 for high might be okay too, but just depends on time and seed quantities.
#3. I'm not sure how i can use the parameters from a logistic curve to judge transmission rates. the b and e pars in the low competency trays are the same as high competency trays and only differ in the upper assymptote. If I wanted to know # of infecteds, I would just do a poisson model. Go with the binomial models. 


rm(list=ls())
source('IBM/scripts/IBM_functions.R')
library(cowplot)
library(brms)
library(tidybayes)
library(lme4)
library(drc)

#Functions-----
dens <- function(x, ...) plot(density(x), ...)
inv_logit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

#allow alteration of tray dimensions
sample_community2 <- function(which_spp, perc_inoc=.1, planting_dist=2, height=24.13, width = 24.13){
  #create size of tray and grid
  tray <- make_tray(width, height)
  grid <- make_grid(r=planting_dist, tray) #2 cm interplanting distance
  spp <- spp[which_spp]#how many species
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



#DESIGN PARAMETERS----
#distances wanted
Dist <- c(1.7)

#spp used (high, med, low competency species)
spp <- c('high', 'med', 'low')
comp <- c(.6, .35, .15) #vector of relative "competencies"
names(comp) <- spp

#MAKE TRAYS----
#interplanting distance is just 1.7cm

#design
#6 species, max of 100 different simulations per species, 1 interplanting distance
design <- expand.grid(species = spp, rep = 1:100, distance = Dist) %>% arrange(species)
design$trayID <- row.names(design)

get_grid <- function(x, W, H){
  design_x <- design[x,]
  which_spp <- which(spp == design_x$species)
  grid <- sample_community2(which_spp, .1, Dist, H, W)
  grid$comp <- comp[design_x$species]
  grid$trayID <- design_x$trayID
  grid$rep <- design_x$rep
  return(grid)
}

#simulate trays 
set.seed(2020)
#normal dimensions
grid_list <- list(NULL)
width = height = 24.13
for(i in 1:nrow(design)){
  grid_list[[i]] <- get_grid(i, width, height)
}

#smaller dimensions
grid_list_S <- list(NULL)
width_S =17
height_S = 16
for(i in 1:nrow(design)){
  grid_list_S[[i]] <- get_grid(i, width_S, height_S)
}

plot_maps(grid_list[[1]]) #238 plants
plot_maps(grid_list_S[[1]]) #110 plants


#SIMULATE SPREAD-----

#model parameters
tfinal <- 25 #how many time steps
beta_curve <- function(x) .2*exp(-3*(log(x/11))^2)
alpha_curve <- function(x) .4*(1-.3)^x
delta <- 1/5 #1/average number of days inoc stays around
beta_ij_t <- make_beta_ij_t(comp) #matrix of amplitudes of the beta_ij 
alpha_i_t <- make_alpha_i_t(comp) #rate of infection from inoculum to plant

#simulate
set.seed(2020)
IBM_list <- lapply(1:length(grid_list), function(x) IBM(grid_list[[x]], 'Kernel', .001))
IBM_list_S <- lapply(1:length(grid_list_S), function(x) IBM(grid_list_S[[x]], 'Kernel', .001))

#plot
ID <- 14
plot_spread_map(grid_list[[ID]], IBM_list[[ID]], animate = T)
plot_grid(plotS_I(IBM_list[[ID]])[[2]],
          plot_spread_map(grid_list[[ID]], IBM_list[[ID]], animate = F))
plot_grid(plotS_I(IBM_list_S[[ID]])[[2]],
          plot_spread_map(grid_list_S[[ID]], IBM_list_S[[ID]], animate = F))



#prep data----
#get dataframe for log-logistic curve fitting
LL_data <- function(IBM_L, grid_L){
  attributes <- grid_L@data[1,c('spID', 'trayID', 'rep')]
  row.names(attributes) <- NULL
  tmp <- IBM_L[!IBM_L[,1]=='C',]#rm primary infections
  data.frame(attributes,
             time = 1:ncol(tmp),
             plants = nrow(tmp), 
             I = colSums(tmp=='I')
  )
}
LLdf <- lapply(1:length(IBM_list), function(x) LL_data(IBM_list[[x]], grid_list[[x]])) %>% bind_rows()
LLdf_S <- lapply(1:length(IBM_list_S), function(x) LL_data(IBM_list_S[[x]], grid_list_S[[x]])) %>% bind_rows()

#get dataframe for binomial model (just last time point)
binomdf <- LLdf %>% filter(time == max(time)) %>% mutate(rep = as.factor(rep))
binomdf_S <- LLdf_S %>% filter(time == max(time))%>% mutate(rep = as.factor(rep))


#Run binomial models-----

#dataframes with different # of reps, split by species
nreps <- rep(3:10, each = 15)
tmp <- split(binomdf, binomdf$spID)
binomdf_list <- list(NULL)#split by species
for(i in 1:length(tmp)){
  binomdf_list[[i]] <- lapply(nreps, function(x) filter(tmp[[i]], rep %in% sample(100, x))) 
}


#repeat for smaller trays
tmp <- split(binomdf_S, binomdf_S$spID)
binomdf_list_S <- list(NULL)#split by species
for(i in 1:length(tmp)){
  binomdf_list_S[[i]] <- lapply(nreps, function(x) filter(tmp[[i]], rep %in% sample(100, x))) 
}



#run loop for lists of data

binom_sensitivity <- function(data_list, Title){
  tmp_mod <- list(NULL)
  tmp_mod <- sapply(data_list, function(x) {
    fit <- glm(cbind(I, plants-I) ~ 1, family = binomial, data = x)
    y <- c(coef(fit), confint(fit))
    names(y) <- c('mean', 'lower', 'upper')
    return(y)
  }) 
  tmp_mod %>% t %>% 
    as.data.frame() %>% 
    mutate(nreps = nreps) %>% 
    ggplot(., aes(nreps, mean)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge2(width = .3), alpha = .5) +
    labs(title = Title) 
}

#visualize
p1 <- binom_sensitivity(binomdf_list[[1]], 'High competency, normal tray')
p2 <- binom_sensitivity(binomdf_list[[3]], 'Med competency, normal tray')
p3 <- binom_sensitivity(binomdf_list[[2]], 'Low competency, normal tray')

p4 <- binom_sensitivity(binomdf_list_S[[1]], 'High competency, small tray')
p5 <- binom_sensitivity(binomdf_list_S[[3]], 'Med competency, small tray')
p6 <- binom_sensitivity(binomdf_list_S[[2]], 'Low competency, small tray')

plot_grid(p1, p2, p3, p4, p5, p6,
          nrow = 2)



#Fit logistic curves-----

#b = slope at inflection point (care most about?)
#c = fixed at zero, lower asymptote
#d = upper asymptote
#e = inflection point

high <- LLdf %>% filter(trayID %in% 1:5)
modhigh <- drm(I ~ time, fct = L.3(), data = high)
summary(modhigh)
#b:(Intercept)  -0.479393   0.037862 -12.661 < 2.2e-16 ***
#d:(Intercept) 179.191965   3.769016  47.543 < 2.2e-16 ***
#e:(Intercept)  13.772216   0.205205  67.114 < 2.2e-16 ***

med <- LLdf %>% filter(trayID %in% 101:105)
modmed <- drm(I ~ time, fct = L.3(), data = med)
summary(modmed)
#b:(Intercept) -0.351294   0.033117 -10.608 < 2.2e-16 ***
#d:(Intercept) 50.831930   2.290480  22.193 < 2.2e-16 ***
#e:(Intercept) 16.562847   0.413777  40.029 < 2.2e-16 ***


low <- LLdf %>% filter(trayID %in% 201:205)
modlow <- drm(I ~ time, fct = L.3(), data = low)
summary(modlow)
#b:(Intercept) -0.47305    0.20256 -2.3354   0.02115 *  
#d:(Intercept)  2.65779    0.30913  8.5975 3.281e-14 ***
#e:(Intercept) 13.12453    1.20233 10.9159 < 2.2e-16 ***

pal <- rev(RColorBrewer::brewer.pal(3, 'YlOrRd'))
plot(high$time, high$I, col=pal[1], pch=16)
lines(predict(modhigh, data.frame(time = 1:25)))
points(med$time, med$I, col=pal[2], pch=16)
lines(predict(modmed, data.frame(time = 1:25)))
points(low$time, low$I, col=pal[3], pch=16)
lines(predict(modlow, data.frame(time = 1:25)))


#try with nls to make sure I"m getting the same answer...well different parameter estimates, but probabaly cuz the funciton is defined differently. Still cant use a single parameter to summarize the speed of transmission. 

## using a selfStart model
nlsHi <- nls(I ~ SSlogis(time, Asym, xmid, scal), high)
summary(nlsHi)
#Asym 179.1970     3.6445   49.17   <2e-16 ***
#xmid  13.7724     0.1964   70.12   <2e-16 ***
#scal   2.0863     0.1625   12.84   <2e-16 ***
  
nlsmed <- nls(I ~ SSlogis(time, Asym, xmid, scal), med)
summary(nlsmed)
#Asym  50.8295     2.3073   22.03   <2e-16 ***
#xmid  16.5623     0.4119   40.21   <2e-16 ***
#scal   2.8463     0.2793   10.19   <2e-16 ***

nlslow <- nls(I ~ SSlogis(time, Asym, xmid, scal), low)
summary(nlslow)
#Asym   2.6578     0.3191   8.329  1.4e-13 ***
#xmid  13.1245     1.2094  10.852  < 2e-16 ***
#scal   2.1139     1.0066   2.100   0.0378 *  



#conclusion: I'm not sure how i can use these parameters to judge transmission rate. the b and e pars in the low competency trays are the same as high competency trays and only differ in the upper assymptote. If I wanted to know # of infecteds, I would just do a poisson model. Go with the binomial models. 

#If I want to know speed, I probably need to use an SIR model. 




#try fitting the logistic curve with the analytical solution?
#dN/dt = rN(1-N/K) #continuous time
#N(t) = (N0*exp(rt)) / (1+N0*(exp(rt) - 1)/K) #analytical solution

#custom function with SSE as loss function
fit.logistic <- function(par,y){
  r <- par[1]; k <- par[2]; n0 <- par[3]
  t <- y[,2]
  n <- y[,1]
  tmp <- k*n0 / (n0 + (k-n0)*exp(-r*t))
  #tmp <- n0 *exp(r*t)/(1 + n0 * (exp(r*t)-1)/k)
  sumsq <- sum((n - tmp)^2)
}

#data (US population census)
head(high)
r.guess <- .2
k.guess <- 200
n0.guess <- 1
par <- c(r.guess,k.guess, n0.guess) #guesses for parameters. necessary to start `optim()`

tmp <- filter(low, time>=0) %>% transmute(I, time)
#tmp <- cbind(high$I, high$time)
tmpfit <- optim(par,fit.logistic,y=tmp)
tmpfit$par #r=0.4789777 K=179.2322471   N0=0.2441312

#plot predictions against observed
logistic.int <- expression(n0 * exp(r*t)/(1 + n0 * (exp(r*t) - 1)/k))
r <- tmpfit$par[1]
k <- tmpfit$par[2]
n0 <- tmpfit$par[3]
t <- 1:25
plot(NULL, xlim=c(0,25), ylim=c(0,214), type="n", xlab="days", 
     ylab="I")
lines(t, eval(logistic.int), lwd=2, col=grey(0.85))
points(tmp$time,tmp$I, pch=16, col="black")

