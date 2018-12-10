#0. load data

#load betas, design df, and species abundances
#B <- read.csv("simulation/outputs/betas.csv")
B <- read.csv("simulation/outputs/rates.csv")
design <- read.csv("simulation/outputs/design.csv")
design[2,2] <- 2.29 #for some reason R sucks and hates 2.28. changing it slightly.
D.add <-design %>% filter(d=="add") %>% pull(r)#(planting distances for addititve)
A <- getabund4(D.add, F, 1) 
A %>% 
  group_by(R) %>% 
  table()
A %>% 
  group_by(R) %>% 
  tally()

C=4
dens="add"
D <-design %>% filter(d==dens) %>% pull(r)#(planting distances)
A <- getabund4(D, F, 1) 
A1 <- A %>% 
  filter(R==unique(R)[C]) %>% 
  pull(A)
A1


######################################################################
#1. format.comp2()
######################################################################

#make function to get host values and matrix of infection probabilities
format.comp2 <- function(rich, beta, distances, plants=abund$R2, steps=20){
  plants=plants %>% as.character() %>% as.numeric()
  distances = (distances*10) #convert them into mm from cm and reverse order so first one is for R1
  #plants=vector of plant identies
  #steps=total time of simulation
  #beta=dataframe of competency values over time. eg. comp2. column values are time, species, comp. length(time)=steps.
  
  #assign values to host
  s <- matrix(rep(plants, steps), ncol=steps) #species matrix across time
  s <- s[sample(nrow(s)),] #randomize cells
  
  #get matrix of beta values for each time step after exposure
  betaP <- beta %>% 
    filter(trans=="prim") %>% 
    filter(dist==2)
  betaS <- beta %>% 
    filter(dist==distances[rich]) %>% 
    filter(trans=="sec") 
 
  valP <- lapply(1:steps, function(x) Rename(betaP[betaP$time==x,]$beta, betaP[betaP$time==x,]$species, from_column = s[,x]))
  h.P <- matrix(unlist(valP), ncol = length(valP), byrow = F)
  valS <- lapply(1:steps, function(x) Rename(betaS[betaS$time==x,]$beta, betaS[betaS$time==x,]$species, from_column = s[,x]))
  h.S <- matrix(unlist(valS), ncol = length(valS), byrow = F)
  
  return(list("betaP"=h.P,"betaS"=h.S, "species_ident"=s))
}


test <- format.comp2(rich = C, beta = B, distances = D, plants = A1, steps = 20)
summary(test$betaS)

######################################################################
#2. #HPraster2()
######################################################################
#set up pathogen and host raster grids

HPraster2 <- function(design, density, rich, h, n=10){
  #design=design
  #density="add"
  #rich=C
  #h=test
  #n=10
  
  design <- design %>% filter(d==density) %>% slice(rich)
  n=round(n*.01*design$n)
  ncol <- design$l
  nrow <- design$w
  #design=data frame with the number of plants in each treatment
  #density="add" or "sub"
  #rich=richness level (1-4)
  #h=output from format.comp()
  #n=number of inoculated plants

  
  host <- path <- raster(ncol=ncol, nrow=nrow, crs="+proj=utm +zone=1 +datum=WGS84", xmn=0, xmx=ncol, ymn=0, ymx=nrow) #create a raster layer for the host and pathogen
  
  #host
  values(host) <- h$species_ident[,1] # species identity, already randomized
  
  #inoculate with pathogen
  values(path) <- 0
  values(path) <- sample(c(rep(1, n), rep(0, ncell(path)-n)), replace = F)

  return(list("path"=path, "host"=host))
}
HP <- HPraster2(design = design, density = dens,rich =  C, h = test, 10)
#plot(HP$path)
#plot(HP$host)
############################################################
#3. simulate()
############################################################
#simulation
simulate2 <- function(h, HP, steps) {
  #h=test
  #HP=HP
  #steps=20
  #fun=f
  #time=1
  #i=88
  
  #setting up the raster grid
  path <- HP$path
  host <- HP$host
  v <- matrix(0, ncell(path), steps)
  v[,1] <- values(path) #values of infection status
  a <- adjacent(host, 1:ncell(host), 4, sorted=TRUE)

  #setting up objects that will be reported
  betaP <- h$betaP
  betaS <- h$betaS
  PK4 <- PK3 <- PK2 <- PK1 <- matrix(rep(betaS[,1], steps), ncol = steps) #"beta" at the current time for each neighbor. change after exposure
  t.exposed <- matrix(NA, nrow=ncell(host), ncol=4) 
  
  for (time in 1:(steps-1)) { #for every time interval
    v[,time+1] <- v[,time] #previous status continues to next time step
    
    #chunk for getting time since exposure for each cell-neighbor combo
    for (i in 1:nrow(v)) { #for every cell i
      if (v[i, time+1] == 1 ) next #skipping those already infected
      adj <- a[a[,1] == i, 2]	#find all your neighbors
      if (any(v[adj, time] > 0)) { #if any of the neighbors are infected
        adjI <- which((v[adj, time]==1)) #identify which ones
        for (k in 1:length(adjI)){ #for each of the infected neighbors
          if (!is.na(t.exposed[i, adjI[k] ])) next#skipping the neighbors that have already been recorded
          t.exposed[i, adjI[k]] <- (time) #record first time exposed for that neighbor
        }
        #if (!is.na(t.exposed[i,k])) next 
        #t.exposed[i,k] <- (time) #record first time exposed for that neighbor
      }
    }
    
    #chunk for turning susceptibles into infected
    for (i in 1:nrow(v)) { #for every cell
      if (v[i, time+1] == 1 ) next #skipping those already infected
      adj <- a[a[,1] == i, 2]	#find all your neighbors
      if (any(v[adj, time] > 0)) { #if any neighbors are infected
        #update probability of infection as a function of time after exposure from each neighbor (k1, k2, k3, k4)
        PK1[i, time] <- betaS[i, time-t.exposed[i, 1]+1] 
        PK2[i, time] <- betaS[i, time-t.exposed[i, 2]+1] 
        PK3[i, time] <- betaS[i, time-t.exposed[i, 3]+1] 
        PK4[i, time] <- betaS[i, time-t.exposed[i, 4]+1] 

        #function for defining susceptible becoming infected  
        f <- function(p, h) {
          if(is.na(h)) { #give NA values 0
            0
          } else{
            if ( runif(1, 0.1, 1) <= h) {#youre infected if value is less than or equal your probability of getting infected (competency value). h=1 always infected. h=0 never infected.
              1
            } else {
              p #if not infected, remain the same as last time step, p
            }
          }
        }
        
        #become infected at some function of current health status and host competency. value of 1 in P leads to infection.
        PK <- list(PK1, PK2, PK3, PK4)
        P <- sapply(1:4, function(k) f(v[i, time], PK[[k]][i, time]))
        v[i, time+1] <- ifelse(1 %in% P, 1, 0)
      } 
    }
  }

  list("values"=v, "betaS"=PK, "t.exposed"=t.exposed)
}
simtest <- simulate2(test, HP, 20)

############################################################
#5. animate()
############################################################
#for visualizing

animate <- function(HP, v, pause=0.1, col=pal) {
  v <- v$values
  h <- HP$host
  p <- HP$path
  pal<-(RColorBrewer::brewer.pal(6, "Spectral")) 
    
  which.color <- unique(values(h)) %>% sort()
  for (i in 1:ncol(v)) { #for every time step
    values(p) <- v[,i] #assign the values of matrix v to raster p at each time step
    plot(h, asp=NA, col=col[which.color], main=i)
    plot(p, legend=FALSE, asp=NA, col=scales::alpha(c(NA, "black"), .8), add=T) #plot it
    dev.flush() #not sure what this does
    Sys.sleep(pause) #give me .25 sec to see each frame
  }
}
animate(HP, simtest, pause = .2)

############################################################
#6. s1()
############################################################
#put it all together

#prep data
prepdata <- function(dens, rand){
  #dens="add" or "sub", rand=F or T
  D <- design %>% filter(d==dens) %>% pull(r)#(planting distances)
  A <- getabund4(D, rand, 1) 
  list(A=A, D=D, dens=dens, rand=rand)
}

#run simulation functions
s1 <- function(Rich, A, beta, design){
  #Rich=richness level, a=abundance dataframe from `prepdata`; beta=transition probs(csv), design=experimental design (csv)
  t <- 20
  n <- 10 #percent inoculated
  
  #load data
  data <- A
  A <- data$A #abundances
  D <- data$D #planting distances
  dens <- data$dens #density treatment
  #rand <- data$rand #randomness
  

  #do simulation for this community at richness level Rich
  A1 <- A %>% 
    filter(R==unique(R)[Rich]) %>% 
    pull(A)
  
  h <- format.comp2(Rich, beta, D, A1, t)
  HP <- HPraster2(design = design, density = dens, Rich, h, n)
  sim <- simulate2(h, HP, t)
  #animate(HP, simtest, pause = .2)
  
  #get response variables
  #RV1 <- response(sim, HP)
  #RV2 <- response.sum(RV1)
  
  #report
  return(list("format.comp"=h, "HPraster"=HP, "simulate"=sim))
 
  
}

#load betas, design df, and species abundances
B <- read.csv("simulation/outputs/rates.csv")
design <- read.csv("simulation/outputs/design.csv")
design[2,2] <- 2.29 #for some reason R sucks and hates 2.28. changing it slightly.

#load data. include density and random treatment
data <- prepdata("sub", F)
#data$A  %>% group_by(R) %>% table()

S1 <- s1(4, data, B, design)
HP <- S1$HPraster
sim <- S1$simulate
animate(HP, sim)
