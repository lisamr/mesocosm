library(dplyr); library(reshape2)
library(ggplot2)

#change color palette
library(RColorBrewer)
pal <- brewer.pal(8, "BrBG")
ggplot <- function(...) ggplot2::ggplot(...) + scale_fill_manual(values=pal) + scale_color_manual(values=pal) 

#generating needed distributions for the mesocosm simulation

#lognormal distribution of species abundances
#important for knowing abundances of species at each richness level.

#you have 6 species. get proportion of individuals per species and standardize to a certain number of total individuals.
getabund <- function(r, n=100, mu=1){ 
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

x=c(1,2,4,6)
counts <- lapply(x, function(r) getabund(r, n=100, mu=1.2))
counts2 <- unlist(counts)
counts2 <- as.data.frame(matrix(counts2, ncol=length(x))) 
names(counts2) <- c("R1", "R2", "R3", "R4")
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
  filter(variable=="R4")
plot(tmp$value, tmp$n)
hist(counts2$R4)

#export host abundances 
write.csv(counts2, "simulation/outputs/hostabund.csv")


