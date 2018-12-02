#basic foundation of cellular simulation

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
  if (sample(h)[1] == h) { #fate decides your infected inverse to your competency value
    1
  } else {
    p #if not infected, remain the same as last time step, p
  }
}

#3. RUN THE FUNCTIONS
vv <- simulate(host, path, steps=20, f)
animate(path, vv, pause=0.25)

###########################
#change host composition
values(host) <- rep(1:4, 25) # competence scores
values(host) <- sample(host, replace = F) #randomize distribution
h <- values(host) #could see values by plot(host)
