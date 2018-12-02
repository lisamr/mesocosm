
library(raster)

# function for display
animate <- function(r, v, pause=0.25) {
	for (i in 1:ncol(v)) { #for every time step
		values(r) <- v[,i] #assign the values of matrix v to raster r at each time step
		plot(r, legend=FALSE, asp=NA, main=i) #plot it
		dev.flush() #not sure what this does
		Sys.sleep(pause) #give me .25 sec to see each frame
	}
}

host <- path <- raster(ncol=10, nrow=10, crs="+proj=utm +zone=1 +datum=WGS84", xmn=0, xmx=10, ymn=0, ymx=10) #create a raster layer for the host and pathogen

# adjacent cells
a <- adjacent(x=host, cells=1:ncell(host), directions=4, sorted=TRUE) #identify cells that are adjacent to a set of cells on the raster. 2 columns: from and to shows from which cell it is a neighbor to.

values(host) <- rep(1:4, 25) # competence scores
h <- values(host) #could see values by plot(host)

values(path) <- 0
path[5,5] <- 1  # one inoculation

steps <- 20 # time steps

# results
v <- matrix(0, ncell(path), steps) #v is 100 x 20. 100 cells, 20 time steps
v[,1] <- values(path) #first time step shows 1 inoculation

# run 1
# ignoring host competence
for (time in 1:(steps-1)) { #for each time step...
	v[,time+1] <- v[,time] #penultimate step assigned as the last time step
	for (i in 1:nrow(v)) { #for each cell in the matrix
		adj <- a[a[,1] == i, 2]	#for each cell in the matrix, get their neighbors
		if (any(v[adj, time] > 0)) { #if any of the cells has a neighbor with a value over 1 at t=time...
			v[i, time+1] <- 1 #turn that cell into value 1 at the next time step
		}
	}
}

animate(path, v, pause=0.25)

# run 2
# host 4 is immune 
for (time in 1:(steps-1)) {
	v[,time+1] <- v[,time]
	for (i in 1:nrow(v)) {
		adj <- a[a[,1] == i, 2]	
		if (any(v[adj, time] > 0)) { #if any neighbor is infected
			if (h[i] < 4) { #AND if the host have a value under 4
				v[i, time+1] <- 1 #it will become infected next time step
			}
		}
	}
}

animate(path, v, pause=0.25)

#change host values
values(host) <- rep(c(1, 7, 8, 10), 25) # competence scores
values(host) <- sample(values(host), replace = F) #randomize the positions of the hosts
h <- values(host) 
plot(host)

# run 3
# infection probability is stochastic; probability is high (1 for host = 1) to low (1/4 for host = 4) 
for (time in 1:(steps-1)) { #for each time step in each cell...
	v[,time+1] <- v[,time]
	for (i in 1:nrow(v)) {
		if (v[i, time+1] == 1 ) next  # skip over the already infected
		adj <- a[a[,1] == i, 2]	#get neighbors of each cell
		if (any(v[adj, time] > 0)) { #if any of the neightbors are infected,
			if (sample(h[i])[1] == h[i]) { #and if fate decides youre infected, (cell becomes infected based on a distribution. distribution is where competency value is equal to odds of getting infected)
				v[i, time+1] <- 1 #turn that cell into an infected at that time.
			}
		}
	}
}

animate(path, v, pause=0.25)



## an attempt to generalize

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

#must define that function
f <- function(p, h) { 
	1 #in this case, always get infected if your neighbor is infected
}
vv <- simulate(host, path, steps=20, f)
animate(path, vv, pause=0.25)

#redefining the function
f <- function(p, h) {
	if (sample(h)[1] == h) { #fate decides your infected inverse to your competency value
		1
	} else {
		p #if not infected, remain the same as last time step, p
	}
}

vv <- simulate(host, path, steps=20, f)
animate(path, vv, pause=0.25)

##################################################
#condensing vital script
##################################################
#DEFINE HOST AND PATHOGEN RASTER
host <- path <- raster(ncol=10, nrow=10, crs="+proj=utm +zone=1 +datum=WGS84", xmn=0, xmx=10, ymn=0, ymx=10) #create a raster layer for the host and pathogen

#host
values(host) <- rep(1:4, 25) # competence scores
values(host) <- sample(host, replace = F) #randomize distribution
h <- values(host) #could see values by plot(host)

#pathogen
values(path) <- 0
path[5,5] <- 1  # one inoculation

#LOAD FUNCTIONS
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

#RUN THE FUNCTIONS
vv <- simulate(host, path, steps=20, f)
animate(path, vv, pause=0.25)
