#libraries and display
library(raster)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
#display.brewer.all()
pal <-rev(brewer.pal(6, "YlOrRd"))  

#simulation evaluation
#I need reports on response variables for each simulation set
#2 functions: 1. `response`=RV for each time step and species within community, 2. `response.sum` = summary response variables to compare across treatments

#Response variables
#max dI/dt; time to max dI/dt; %infected; all variables reported for each species and all together. 
#G matrix (communitity R0) = dominant eiganvalue of the N x N matrix. a=Bij*Pij/di. components are transmission rates and average duration of infection period. B=transmission rates, P=#infecting species, d=death rate, which I assume is 1 for this system. 

#load data
comp2 <- read.csv("simulation/outputs/competencyvalues.csv")
abund <- read.csv("simulation/outputs/hostabund.csv")

########################################################################
#percent infected
########################################################################
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

##############################################################
#Summary Response variables
##############################################################
head(RV)

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

##############################################################
#test it!
##############################################################

#run simulation
h1 <- format.comp(plants = abund$R4)
HP <- HPraster(h1, n = 10)
vv <- simulate(h1, HP, steps=20, f)
#animate(HP, vv,pause = .1, col = pal)

#get results
RV <- response(vv, HP)
RVsum <- response.sum(RV)
head(RVsum)
head(RV)

#plot
ggplot(RV, aes(time, n.I, group=species, col=species))+
  geom_point()+
  geom_line()
ggplot(RV, aes(time, dI, group=species, col=species))+
  geom_point()+
  geom_line()










