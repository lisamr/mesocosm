library(dplyr); library(reshape2)
library(ggplot2)

#change color palette
library(RColorBrewer)
pal <- brewer.pal(8, "BrBG")
ggplot <- function(...) ggplot2::ggplot(...) + scale_fill_manual(values=pal) + scale_color_manual(values=pal) 

#generating needed distributions for the mesocosm simulation

#################################
#SPECIES DISTRIBUTIONS
#################################
#lognormal distribution of species abundances
#important for knowing abundances of species at each richness level.

#you have 6 species. get proportion of individuals per species and standardize to a certain number of total individuals.


#adding in argument "rand" so I can specify if assembly is random
getabund3 <- function(n=100, mu=1, rand=F) {
  
  #load nested functions
  getabund <- function(r, n, mu){ 
    #r=richness, n=total individuals, mu=mean of log distribution
    
    #get counts of each species
    n1 = n+max(r) #buffering the number of species so there aren't errors due to rounding
    pool=dlnorm(1:6, mu) #distribution of each species
    abund <- pool[1:r]*n1/sum(pool[1:r]) #get abundances of each species for a richness level
    
    #get actual population from the counts above 
    r2=1:length(r)
    pop <- rep(1:length(abund), abund)
    pop <- sort(sample(pop, n)) #bring the population back down to desired n. had issues with rounding.
    return(pop)
  }
  getabund2 <- function(n, mu){
    x=c(1,2,4,6)
    counts <- lapply(x, function(r) getabund(r, n, mu))
    counts2 <- unlist(counts)
    counts2 <- as.data.frame(matrix(counts2, ncol=length(x))) 
    names(counts2) <- c("R1", "R2", "R4", "R6")
    counts2
  }
  
  #setup randomization option
  A <- getabund2(n, mu)
  sp <- unique(A$R6)
  if (rand==T) {
    randsp <- sample(sp, length(sp))
    A <- apply(A, 2, function(x) Rename(randsp, sp, x))
  }
    A
}

getabund3(200, 1, rand = T)  




melt(counts2)

#plot counts
counts3 <- melt(counts2)
spp <- as.factor(counts3$value)
spp <- factor(spp, levels = sort(levels(spp), T))
ggplot(counts3, aes( variable, group=spp, fill=spp)) +
  geom_bar()
#plot distribution of richest community
tmp <- counts3 %>% 
  group_by(variable) %>% 
  count(value) %>% 
  filter(variable=="R6")
plot(tmp$value, tmp$n)
hist(counts2$R4)

#export host abundances 
write.csv(counts2, "simulation/outputs/hostabund.csv", row.names = F)

##################################
#COMPETENCY CURVES
##################################
#competency is a function of susceptibility and infectivity. if we assume infectivity is constant over time and across species, it's just a function of susceptibility. I expect a plant to be susceptible highest a day or two after their neighbor is infected and a gradual decay over time. The distribution is from Otten 2003, describing the decay of the secondary transmission rates. I'm eyeballing the distribution to match something like the decay in primary transmission, same paper. formula is as follows: 

trans <- function(B=.182, g=3.06, t=1:10, tq){
  #B determines amplitude; gamma is slope of decay; t is time duration; tq is time of peak
  B_t = B*exp(-g*(log(t/tq))^2)
  return(B_t)
}

#explore the parameters
#Beta: amplitude
tmp <- lapply(seq(.1, .5, length.out = 10), function(x) trans(tq=2, B=x))
plot(1:10, tmp[[10]])
for(i in 1:length(tmp)){
  points(1:10, tmp[[i]], type='o')
}

#Gamma: steepness
tmp <- lapply(seq(.1, .5, length.out = 10), function(x) trans(tq=2, g=x))
plot(1:10, tmp[[10]])
for(i in 1:length(tmp)){
  points(1:10, tmp[[i]], type='o')
}

#tq: time of peak
tmp <- lapply(1:5, function(x) trans(tq=x))
plot(1:10, tmp[[1]])
for(i in 1:length(tmp)){
  points(1:10, tmp[[i]], type='o')
}

#it looks like tq=2 seems most sensible. B and g are from Otten 2003. Will keep g as is, but B will vary depending on plant species. 
tmp <- lapply(c(c(.9, .4, .3, .1, 0, 0)), function(x) trans(tq=2, B=x, t=1:20))
plot(1:20, tmp[[1]])
for(i in 1:length(tmp)){
  points(1:20, tmp[[i]], type='o')
}
#create dataframe
comp <- data.frame(time=1:20, do.call("cbind", tmp))
names(comp) <- c("time", 1:6)
comp2 <- melt(comp, measure.vars = 2:7, id.vars = "time", variable.name = "species", value.name = "comp")
comp2

#export "competency" values
write.csv(comp2, "simulation/outputs/competencyvalues.csv", row.names = F)

