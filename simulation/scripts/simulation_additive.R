#0. load data

#load betas, design df, and species abundances
#B <- read.csv("simulation/outputs/betas.csv")
B1 <- read.csv("simulation/outputs/rates111111.csv")
B2 <- read.csv("simulation/outputs/rates1.7.5.4.2.1.csv")
B3 <- read.csv("simulation/outputs/rates1.3.2.100.csv")

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
  t <- 25
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
data <- prepdata("add", T)
#data$A  %>% group_by(R) %>% table()

S1 <- s1(2, data, B, design)
HP <- S1$HPraster
sim <- S1$simulate
animate(HP, sim)

############################################################
#7. response()
############################################################
response <- function(sim, HP){
  sim=sim
  HP=HP
  #data=output from simulate; HP=output from HPraster
  data <- sim
  #get data
  dat <- data$values #state of each cell by time: N x time
  spp <- as.vector(values(HP$host)) #vector of species present in correct
  
  #get table of frequency of each species 
  fr <- as.data.frame(table(spp))
  N <- cbind.data.frame("tot", as.numeric(sum(fr$Freq)))
  names(N) <- names(fr) <- c("species", "n")
  fr <- rbind(fr, N)
  
  ############################################################
  #make tables of infecteds
  ############################################################
  #for each time step, get sum of infecteds and cbind onto existing f dataframe
  sp <- fr$species[1:length(fr$species)-1]
  I <- (dat * spp)
  frI <- fr
  for(t in 1:ncol(I)){
    tmp <- sapply(sp, function(i) {
      sum(I[,t]==sp[i])
    })
    tmp[length(tmp)+1] <- sum(tmp)
    frI <- cbind(frI, tmp)
  }
  names(frI) <- c(names(fr), 1:ncol(I))
  
  ############################################################
  #make tables of exposed
  ############################################################
  #for each time step, get sum of exposed and cbind onto existing f dataframe
  
  #find time infected for every cell
  t.inf <- sapply(1:nrow(dat), function(i) which(dat[i, ]==1)[1])
  
  #table of exposure
  #first time exposed for each cell
  data$t.exposed[is.na(data$t.exposed)] <- 99 #change all NAs to 99 so it doesn't have problems finding min
  t.exposed1 <- apply(data$t.exposed, 1, min)
  t.exposed1[t.exposed1==99] <- NA #change them back to NA
  exposure <- data.frame(cell=1:length(t.exposed1), exp.i=t.exposed1, exp.f=t.inf-1)
  #$exp.f==NA when never infected. $exp.i==NA when never exposed. NA's introduced when 1. never exposed, 2. inoculated, 3. exposed, but not infected
  #1. never exposed
  S <- is.na(exposure$exp.i) & is.na(exposure$exp.f)
  #2. inoculated
  I <- is.na(exposure$exp.i) & exposure$exp.f==0
  #3. exposed, never infected
  E <- exposure$exp.i>0 & is.na(exposure$exp.f)
  #change values
  exposure$exp.i[which(S)] <- 0
  exposure$exp.i[which(I)] <- 0
  exposure$exp.f[which(S)] <- 0
  exposure$exp.f[which(I)] <- 0
  exposure$exp.f[which(E)] <- ncol(dat)
  
  for(i in 1:nrow(dat)){
    #change these values to NA, corresponding to time exposed
    dat[i, exposure$exp.i[i]:exposure$exp.f[i]] <- NA 
  }
  
  #make matrix for just exposed plants
  e <- dat
  e[!is.na(e)] <- 0
  e[is.na(e)] <- 1
  
  #get frequencies of exposed plants
  E <- e * spp
  
  frE <- fr
  for(t in 1:ncol(E)){
    tmp <- sapply(sp, function(i) {
      sum(E[,t]==sp[i])
    })
    tmp[length(tmp)+1] <- sum(tmp)
    frE <- cbind(frE, E=tmp)
  }
  names(frE) <- c(names(fr), 1:ncol(E))
  
  
  ############################################################
  #make table tall
  ############################################################
  frI2 <- melt(frI, c("species", "n"), 3:ncol(frI), variable.name = "time", value.name = "n.I")
  frE2 <- melt(frE, c("species", "n"), 3:ncol(frE), variable.name = "time", value.name = "n.E")
  fr2 <- left_join(frI2, frE2, by=c("species", "time", "n")) %>% arrange( species) 
  
  ############################################################
  #get n.S and %inf
  ############################################################
  fr2$n.S <- fr2$n-fr2$n.I
  fr2$pI <- fr2$n.I/fr2$n
  
  ############################################################
  #get dI/dt
  ############################################################
  #dI/dt = I[t+1]-I[t]
  di <- function(t) fr2[t,]$n.I-fr2[t-1,]$n.I
  fr2$dI <- c(NA, di(2:length(fr2$n.I))) 
  
  #change all dI where time=1 to NA
  fr2$dI[fr2$time==1] <- NA
  
  #time needs to be numerical
  fr2$time <- as.numeric(fr2$time)
  
  return(fr2)
}
RV <- response(sim, HP)

#############################################################
#8. response.sum()
#############################################################
#Summary Response variables
response.sum <- function(resp){
  #n.species, n.infected, perc.infected for each species
  RV1 <- resp %>% 
    filter(time==max(time)) %>% 
    group_by(species) %>% 
    summarize(n, n.I, pI)
  
  #dIdt reports
  RV2 <- resp %>% 
    group_by(species) %>% 
    filter(dI==max(dI, na.rm = T)) %>% 
    filter(!duplicated(dI)) %>% #remove rows with duplicated dI
    summarise(max.dI=dI, tmax.dI=time)
  
  left_join(RV1, RV2, by="species")
}
RV2 <- response.sum(RV)

#############################################################
#8. s2()
#############################################################
#simulation function that packages together simulation functions and responses
s2 <- function(Rich, A, beta, design){
  #Rich=richness level, a=abundance dataframe from `prepdata`; beta=transition probs(csv), design=experimental design (csv)
  t <- 15
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
  RV1 <- response(sim, HP)
  RV2 <- response.sum(RV1)
  
  #report
  return(list("format.comp"=h, "HPraster"=HP, "simulate"=sim, "resp1"=RV1, "resp2"=RV2))
}

#load data. include density and random treatment
data <- prepdata("add", T)
#data$A  %>% group_by(R) %>% table()

S2 <- s2(3, data, B, design)
RV <- S2$resp1

#plot
animate(S2$HPraster, S2$simulate, .1)
ggplot(RV, aes(time, n.I, group=species, col=species))+
  geom_point()+
  geom_line()
ggplot(RV, aes(time, dI, group=species, col=species))+
  geom_point()+
  geom_line()
