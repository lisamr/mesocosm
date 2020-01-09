#functions in a r script 
library(raster)
library(dplyr)
library(ggplot2)
library(reshape2)
library(viridis)
library(cowplot)


#### load data ####
#load betas(sort of ~transition rates), design df
#see "oldcode_transitionprobs.R" and "species_distribution.R" for how I got the data below.
#B <- read.csv("simulation/outputs/betas.csv")
B <- read.csv("simulation/outputs/rates.csv")
design <- read.csv("simulation/outputs/design.csv")
design[2,2] <- 2.4 #for some reason R sucks and hates 2.39. changing it slightly.
dat <- read.csv("simulation/outputs/moretime_responses.csv")#substitutve has 392 individuals and is run for 30 time steps. See 'moretime_fitnls.R' for code

#Data managment ####
Rename <-  function(values_to, index_from, from_column){
  result <- (values_to[match(from_column, index_from)])
  return(result)
}


#indivual components for simulation #####
#getting abundances for all 4 treatments: sub/det, sub/rand, add/det, add/rand
#default planting distances (cm)
D.add <- rev(c(1.75, 1.9, 2.28, 3))
D.sub <- 2
getabund4 <- function(dist, L=50, rand=F, SD=.5){
  #load function
  getabund <- function(r, n){ 
    #r=richness, n=total individuals
    
    #get counts of each species
    n1 = n+max(r) #buffering the number of species so there aren't errors due to rounding
    
    #get distribution of each species. comes from a log distribution (mu=6, sd=.4) and 6 species are pulled. replicated 1000 times and the mean values taken.
    N <- 6
    sims <- 1000
    m2 <- matrix(NA, nrow = sims, ncol = N)
    for(i in 1:sims){
      y2 <- rlnorm(1:N, meanlog = 1, sdlog = SD) %>% 
        sort(decreasing = T)
      m2[i,] <- y2
    }
    pool <- apply(m2, 2, mean)/sum(apply(m2, 2, mean))
    
    #get abundances of each species for a richness level
    abund <- pool[1:r]*n1/sum(pool[1:r])
    
    #get actual population from the counts above 
    r2=1:length(r)
    pop <- rep(1:length(abund), abund)
    pop <- sort(sample(pop, n)) #bring the population back down to desired n. had issues with rounding.
    return(pop)
  }
  
  #design is based on dimensions of available propogation trays
  design <- tibble(r=dist, w=floor(25/r), l=floor(L/r), n=w*l, sp=c( 1,2,4,6))
  tmp <- lapply(c(1,2,4,6), function(i) {
    A <- getabund(i, design$n[design$sp==i])
    R <- rep(i, length(A))
    cbind(R, A)
  }
  )
  tmp <- do.call("rbind.data.frame", tmp)
  
  #setup randomization option
  sp <- unique(tmp$A)
  if (rand==T) {
    randsp <- sample(sp, length(sp))
    tmp$A <- Rename(randsp, sp, tmp$A)
  }
  tmp$R <- as.factor(as.numeric(tmp$R))
  tmp$A <- as.factor(as.numeric(tmp$A))
  tmp$A <- factor(tmp$A, levels = rev(levels(tmp$A)))
  tmp
}

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

#for visualizing
animate <- function(HP, v, pause=0.1, col=pal, saveplot=F, Name=NA) {
  v <- v$values
  h <- HP$host
  p <- HP$path
  #pal <- viridis::viridis(6,end = .95, direction = -1)
  pal<-(RColorBrewer::brewer.pal(6, "Spectral")) 
  
  which.color <- unique(values(h)) %>% sort()
  for (i in 1:ncol(v)) { #for every time step
    values(p) <- v[,i] #assign the values of matrix v to raster p at each time step
    plot(h, asp=1, col=col[which.color], main=i)#asp keeps 1:1 aspect
    
    if(saveplot==T){
      # 2. Create a plot
      plot(p, legend=FALSE, asp=NA, col=scales::alpha(c(NA, "black"), .8), add=T) #plot it
      # Open a png file of current dimensions and assign name
      dims <- dev.size(units = "in")
      name1 <- paste0(Name, i)
      name2 <- paste0(paste0("simulation/images/animation/", name1), ".png")
      dev.print(png, name2, units='in', res=150, width = dims[1], height = dims[2] )
      
    }else{
      plot(p, legend=FALSE, asp=NA, col=scales::alpha(c(NA, "black"), .8), add=T) #plot it
    }
    
    dev.flush() #not sure what this does
    Sys.sleep(pause) #give me .25 sec to see each frame
  }
}

#getting response variables out
response <- function(sim, HP){
  
  #data=output from simulate; HP=output from HPraster
  data <- sim
  #get data
  dat <- data[[1]] #state of each cell by time: N x time
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
  exposure <- data.frame(cell=1:length(data$t.exposed), exp.i=data$t.exposed, exp.f=t.inf-1)
  #collapse table into 3 columns: cell, exp.i, exp.f
  e.i <- apply(exposure[,2:5], 1, function(x) {min(x, na.rm = T)})
  e.i[is.infinite(e.i)] <- NA
  exposure <- data.frame(cell=1:length(data$t.exposed), exp.i=e.i, exp.f=t.inf-1)
  
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
  
  return((fr2))
}
#getting Force of infection
FOI <- function(sim, HP){
  
  #data=output from simulate; HP=output from HPraster
  data <- sim
  #get data
  dat <- data[[1]] #state of each cell by time: N x time
  spp <- as.vector(values(HP$host)) #vector of species present in correct
  
  #find time infected for every cell
  t.inf <- sapply(1:nrow(dat), function(i) which(dat[i, ]==1)[1])
  
  #table of exposure
  exposure <- data.frame(cell=1:nrow(data$t.exposed), exp.i=data$t.exposed, exp.f=t.inf-1)
  #collapse table into 3 columns: cell, exp.i, exp.f
  e.i <- apply(exposure[,2:5], 1, function(x) {min(x, na.rm = T)})
  e.i[is.infinite(e.i)] <- NA
  exposure <- data.frame(cell=1:nrow(data$t.exposed), exp.i=e.i, exp.f=t.inf-1)
  
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
  
  #get time exposed for each species
  exposure$exp.time <- exposure$exp.f - exposure$exp.i
  exposure$exposed <- ifelse(exposure$exp.time>0, 1, 0)
  exposure$infected <- ifelse(exposure$exp.f +1 <= ncol(dat), 1, 0)
  exposure$species <- spp
  head(exposure)
  
  #FOI = new infections/(# exposed * duration of exposure)
  FOI <- exposure %>% 
    group_by(species) %>% 
    summarise(n.E=sum(exposed),  FOI = sum(infected)/(sum(exposed) * mean(exp.time))) %>% 
    ungroup() %>% 
    bind_rows(summarise(.data = exposure, n.E=sum(exposed), FOI = sum(infected)/(sum(exposed) * mean(exp.time))))
  FOI$species[nrow(FOI)] <- "tot"
  return(FOI)
}

#summary of response variables
response.sum <- function(resp, foi){
  
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
  
  sum <- left_join(RV1, RV2, by="species")
  
  #FOI reports
  left_join(sum, foi, by = "species") 
  
}


#simulate! ####
#simulation function that packages together simulation functions and responses
s2 <- function(Rich, A, beta, design, t=20){
  #Rich=richness level, a=abundance dataframe from `prepdata`; beta=transition probs(csv), design=experimental design (csv)
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
  foi <- FOI(sim, HP)
  RV2 <- response.sum(RV1, foi)
  
  #report
  return(list("format.comp"=h, "HPraster"=HP, "simulate"=sim, "resp1"=RV1, "resp2"=RV2))
}

#functions for binding different designs. see moretreatments.R for example
prepdata2 <- function(dens, rand, L=50, design, SD=.5){
  #dens="add" or "sub", rand=F or T
  D <- design %>% filter(d==dens) %>% pull(r)#(planting distances)
  A <- getabund4(D,L, rand, SD) 
  list(A=A, D=D, dens=dens, rand=rand)
}

#simulation and results functions ("f.sim" below)
results3 <- function(f.sim=s2, t=20, B, L=50, nreps=10, density="sub", Random=F, design, SD=.5){
  #density="add" or "sub", Rand=T or F, f.sim=s2
  R <- c(1,2,4,6)
  
  #add binding functions
  bind.response <- function(r){
    #response
    tmp1 <- lapply(1:4, function(x) do.call("rbind", r[[x]][4]))
    tmp2 <- lapply(1:4, function(x) cbind(ident[x,], tmp1[[x]]))
    plyr::ldply(tmp2, data.frame)
  }
  bind.respsum <- function(r){
    #response summary
    tmp4 <- lapply(1:4, function(x) do.call("rbind", r[[x]][5]))
    tmp5 <- lapply(1:4, function(x) cbind(ident[x,], tmp4[[x]]))
    plyr::ldply(tmp5, data.frame)
  }
  
  #1. response
  bres <- list(NULL)
  bres2 <- list(NULL)
  for(i in 1:nreps){
    data <- prepdata2(density, Random,L, design, SD)
    res1 <- lapply(1:4, function(R) s2(R, data, B, design, t))
    ident <- data.frame(rep=i, dens=density, rand=Random, rich=R)
    bres[[i]] <- bind.response(res1)
    bres2[[i]] <- bind.respsum(res1)
  }
  bres <- plyr::ldply(bres, data.frame)
  bres2 <- plyr::ldply(bres2, data.frame)
  
  return(list("responses"= bres, "response.summary"= bres2))
}

#multiple replicates ####
#Getting results out for multiple replicates. Provided below are functions for the original experimental design, and one opened up to be more flexible for binding together different designs.

#simulate lots of replicates
fb3 <- function(Beta, t, L=50, n=10, design, sd=.5){
  SD <-  results3(s2,t,L, B=Beta, nreps=n, "sub", F, design, sd) #substitutive, deterministic
  SR <-  results3(s2,t,L, B=Beta, nreps=n, "sub", T, design, sd) #substitutive, random
  AD <-  results3(s2,t,L, B=Beta, nreps=n, "add", F, design, sd) #additive, deterministic
  AR <-  results3(s2,t,L, B=Beta, nreps=n, "add", T, design, sd) #additive, random
  response <- rbind(SD$responses, SR$responses, AD$responses, AR$responses)
  response.summary <- rbind(SD$response.summary, SR$response.summary, AD$response.summary, AR$response.summary)
  
  return(list(responses=response, response.summary=response.summary))
}


