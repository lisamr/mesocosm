#basic foundation of cellular simulation
library(raster)
library(dplyr)
library(ggplot2)
library(reshape2)
Rename <-  function(values_to, index_from, from_column){
  result <- (values_to[match(from_column, index_from)])
  return(result)
}

#1. DEFINE HOST AND PATHOGEN RASTER
host <- path <- raster(ncol=10, nrow=10, crs="+proj=utm +zone=1 +datum=WGS84", xmn=0, xmx=10, ymn=0, ymx=10) #create a raster layer for the host and pathogen

#host
values(host) <- rep(1:4, 25) # competence scores
values(host) <- sample(host, replace = F) #randomize distribution
h <- values(host) #could see values by plot(host)

#pathogen
values(path) <- 0
path[5,5] <- 1  # one inoculation

#2. LOAD FUNCTIONS
#for visualizing
animate <- function(r, v, pause=0.25) {
  for (i in 1:ncol(v)) { #for every time step
    values(r) <- v[,i] #assign the values of matrix v to raster r at each time step
    plot(r, legend=FALSE, asp=NA, main=i) #plot it
    dev.flush() #not sure what this does
    Sys.sleep(pause) #give me .25 sec to see each frame
  }
}

#for running simulation
simulate <- function(host, path, steps, fun) {
  h <- values(host)
  
  v <- matrix(0, ncell(path), steps)
  v[,1] <- values(path)
  a <- adjacent(host, 1:ncell(host), 4, sorted=TRUE)
  
  for (time in 1:(steps-1)) {
    v[,time+1] <- v[,time]
    for (i in 1:nrow(v)) {
      if (v[i, time+1] == 1 ) next
      adj <- a[a[,1] == i, 2]	
      if (any(v[adj, time] > 0)) {		#if any of the neighbors are infected
        v[i, time+1] <- fun(v[i, time], h[i]) #become infected at some function of current health status and host competency
      }
    }
  }
  v
}

#how to determine if cell is infected
f <- function(p, h) {
  if ( runif(1, 0.1, 1) <= h) {#youre infected if value is less than or equal your probability of getting infected (competency value). h=1 always infected. h=0 never infected.
    1
  } else {
    p #if not infected, remain the same as last time step, p
  }
}

#3. RUN THE FUNCTIONS
vv <- simulate(host, path, steps=20, f)
animate(path, vv, pause=0.1)

###########################
#change host composition
#load abundances
abund <- read.csv("simulation/outputs/hostabund.csv")

#make a competency table to relate to the abund df
c <- data.frame(species=1:6, max.comp=c(1, .7, .3, .1, 0, 0))

#assign new competency values to hosts
val <- apply(abund, 2, function(x) Rename(c$max.comp, c$species, x))
values(host) <- val[,4] # competence scores
values(host) <- sample(host, replace = F) #randomize distribution

plot(host, col=rev(heat.colors(6)))
#3. RUN THE FUNCTIONS
vv <- simulate(host, path, steps=20, f)
animate(path, vv, pause=0.1)

###########################
#change number of inoculations
values(path) <- sample(c(rep(1, 10), rep(0, 90)), replace = F)

#3. RUN THE FUNCTIONS
vv <- simulate(host, path, steps=20, f)
animate(path, vv, pause=0.05)

#########################################################################
#Infection probability changes after time from exposed
#########################################################################
#change host competency as a function of time
#exact distribution can be determined emperically: (P(I) vs. time since neighbor infected)
#load abundances
abund <- read.csv("simulation/outputs/hostabund.csv")
counts3 <- melt(abund)
spp <- as.factor(counts3$value)
spp <- factor(spp, levels = sort(levels(spp), T))
ggplot(counts3, aes( variable, group=spp, fill=spp)) +
  geom_bar()

#read in estimated competency distributions. these are really infection probabilities.
comp2 <- read.csv("simulation/outputs/competencyvalues.csv")
head(comp2)
ggplot(comp2, aes(x=time, y=comp, group=species, col=as.factor(species)))+
  geom_point()+
  geom_line()

#make function to get host values and matrix of infection probabilities
format.comp <- function(compdf=comp2, plants=abund$R2, steps=20){
  #plants=vector of plant identies
  #steps=total time of simulation
  #compdf=dataframe of competency values over time. eg. comp2. column values are time, species, comp. length(time)=steps.
  
  Rename <-  function(values_to, index_from, from_column){
    result <- (values_to[match(from_column, index_from)])
    return(result)
  }
  
  #assign values to host
  s <- matrix(rep(plants, steps), ncol=steps) #species matrix across time
  s <- s[sample(nrow(s)),] #randomize cells
  
  #get matrix of competency values for each time step after exposure
  #label it "h"
  val <- lapply(1:steps, function(x) Rename(compdf[compdf$time==x,]$comp, compdf[compdf$time==x,]$species, from_column = s[,x]))
  h <- matrix(unlist(val), ncol = length(val), byrow = F)
  return(list("infection_prob"=h, "species_ident"=s))
}

#make function to set up pathogen and host raster grids

#1. DEFINE HOST AND PATHOGEN RASTER
HPraster <- function(h, n=10, ncol=10, nrow=10){
  #h=output from format.comp()
  #n=number of inoculated plants
  #ncol/nrow dimensions of raster grid
  
  host <- path <- raster(ncol=ncol, nrow=nrow, crs="+proj=utm +zone=1 +datum=WGS84", xmn=0, xmx=10, ymn=0, ymx=10) #create a raster layer for the host and pathogen
  
  #host
  values(host) <- h1[[2]][,1] # species identity, already randomized
  
  #inoculate with pathogen
  values(path) <- 0
  values(path) <- sample(c(rep(1, n), rep(0, ncell(path)-n)), replace = F)
  
  return(list("path"=path, "host"=host))
}

#need to add: 
#1. time since neighbor infected; 
#2. h is drawn from matrix above at the correct time.

simulate <- function(h, HP, steps, fun) {
#setting up the raster grid
  #h=h1
  #HP=HP
  #steps=20
  #fun=f
    
  h <- h[[1]]  
  path <- HP$path
  host <- HP$host
  v <- matrix(0, ncell(path), steps)
  v
  v[,1] <- values(path) #values of infection status
  a <- adjacent(host, 1:ncell(host), 4, sorted=TRUE)
  com.rates <- matrix(rep(h[,1], steps), ncol = steps) #probability of infection rates over time, aka competency. change after exposure
  t.exposed <- rep(NA, ncell(host))
  
  for (time in 1:(steps-1)) { #for every time interval
    v[,time+1] <- v[,time] #previous status continues to next time step
    
    #chunk for getting time since exposure for each cell
    for (i in 1:nrow(v)) { #for every cell
      if (v[i, time+1] == 1 ) next #skipping those already infected
      adj <- a[a[,1] == i, 2]	#find all your neighbors
      if (any(v[adj, time] > 0)) { #if any of the neighbors are infected
        if ( !is.na(t.exposed[i]) ) next #skipping those already exposed
        t.exposed[i] <- (time) #record first time exposed
      }
    }
    
    #chunk for turning susceptibles into infected
    for (i in 1:nrow(v)) { #for every cell
     if (v[i, time+1] == 1 ) next #skipping those already infected
      adj <- a[a[,1] == i, 2]	#find all your neighbors
      if (any(v[adj, time] > 0)) { #if any neighbors are infected
        com.rates[i, time] <- h[i, time-t.exposed[i]+1] #update probability of infection as a function of time after exposure
        v[i, time+1] <- fun(v[i, time], com.rates[i, time]) #become infected at some function of current health status and host competency
      } 
    }
  }
  list("values"=v, "Inf.prob"=com.rates, "t.exposed"=t.exposed)
  }

#how to determine if cell is infected
f <- function(p, h) {
  if ( runif(1, 0.1, 1) <= h) {#youre infected if value is less than or equal your probability of getting infected (competency value). h=1 always infected. h=0 never infected.
    1
  } else {
    p #if not infected, remain the same as last time step, p
  }
}

#for visualizing
animate <- function(HP, v, pause=0.25, col=palette()) {
  h <- HP$host
  p <- HP$path
  
  for (i in 1:ncol(v)) { #for every time step
    values(p) <- v[,i] #assign the values of matrix v to raster p at each time step
    plot(h, asp=NA, col=col, main=i)
    plot(p, legend=FALSE, asp=NA, col=scales::alpha(c(NA, "black"), .8), add=T) #plot it
    dev.flush() #not sure what this does
    Sys.sleep(pause) #give me .25 sec to see each frame
  }
}

#run simulation
#library(RColorBrewer)
display.brewer.all()
pal <-rev(brewer.pal(6, "YlOrRd"))  

set.seed(2)
h1 <- format.comp(plants = abund$R4)
HP <- HPraster(h1, n = 10)
vv <- simulate(h1, HP, steps=20, f)

animate(HP, vv[[1]], pause=0.2, col=pal)
plot(HP$host, asp=NA, col=pal)
plot(HP$path, asp=NA, col=c("white", "black"))

