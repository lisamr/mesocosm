#Functions for mesocosm simulation
#Run in this order: 
#1. format.comp(), 2. HPraster(),  3. simulate(h=format.comp, HP = HPraster, fun = f), 4. f(), 5. animate(HP = HPraster, v = simulate)

#libraries and display
library(raster)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
#display.brewer.all()
pal <-rev(brewer.pal(6, "YlOrRd"))  

#data management functions
extract_after <-  function(char, vector){
    pattern <- paste0(".*", char)
    sub(pattern, '', vector)#after chacter
  }

############################################################
#1. format.comp()
############################################################
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
  val <- lapply(1:steps, function(x) Rename(compdf[compdf$time==x,]$comp, compdf[compdf$time==x,]$species, from_column = s[,x]))
  h <- matrix(unlist(val), ncol = length(val), byrow = F)
  return(list("infection_prob"=h, "species_ident"=s))
}

############################################################
#2. HPraster()
############################################################
#set up pathogen and host raster grids
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

############################################################
#3. simulate()
############################################################
#simulation
simulate <- function(h, HP, steps, fun) {
  #setting up the raster grid
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

############################################################
#4. f()
############################################################
#how to determine if cell is infected
f <- function(p, h) {
  if ( runif(1, 0.1, 1) <= h) {#youre infected if value is less than or equal your probability of getting infected (competency value). h=1 always infected. h=0 never infected.
    1
  } else {
    p #if not infected, remain the same as last time step, p
  }
}

############################################################
#5. animate()
############################################################
#for visualizing
animate <- function(HP, v, pause=0.25, col=palette()) {
  v <- v[[1]]
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
