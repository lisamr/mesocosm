library(dplyr); library(reshape2)
library(ggplot2)

#change color palette
library(RColorBrewer)
library(viridis)
pal<-(RColorBrewer::brewer.pal(6, "Spectral")) 
pal <- rev(pal)


#generating needed distributions for the mesocosm simulation

#################################
#SPECIES DISTRIBUTIONS
#################################
#lognormal distribution of species abundances
#important for knowing abundances of species at each richness level.

#you have 6 species. get proportion of individuals per species and standardize to a certain number of total individuals.

getabund <- function(r, n){ 
  #r=richness, n=total individuals, mu=mean of log distribution
  
  #get counts of each species
  n1 = n+max(r) #buffering the number of species so there aren't errors due to rounding
  
  #get distribution of each species. comes from a log distribution (mu=6, sd=.4) and 6 species are pulled. replicated 1000 times and the mean values taken.
  N <- 6
  sims <- 1000
  m2 <- matrix(NA, nrow = sims, ncol = N)
  for(i in 1:sims){
    y2 <- rlnorm(1:N, meanlog = 1, sdlog = .5) %>% 
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

#adding in argument "rand" so I can specify if assembly is random
getabund3 <- function(n=100, rand=F) {
  
  #load nested functions
  getabund2 <- function(n){
    x=c(1,2,4,6)
    counts <- lapply(x, function(r) getabund(r, n))
    counts2 <- unlist(counts)
    counts2 <- as.data.frame(matrix(counts2, ncol=length(x))) 
    names(counts2) <- c("R1", "R2", "R4", "R6")
    counts2
  }
  
  #setup randomization option
  A <- getabund2(n)
  sp <- unique(A$R6)
  if (rand==T) {
    randsp <- sample(sp, length(sp))
    A <- apply(A, 2, function(x) Rename(randsp, sp, x))
  }
    A
}

counts2 <- getabund3(200, rand = F)  
counts2
#pool comes from a lognormal distribution
getabund(6, 200) %>% table %>% as.data.frame() %>% plot
getabund(6, 200) %>% hist
#################################
#Abundances with additive treatment
#################################
#trays are 10 (25) x 20 (50) inches, optimumum distance is 2 cm apart...12x25=300plants
#design abundances so abund vs. richness increases linearly. Unfortunately, number of species n doesn't remain the same throughout treatments. may need to be a saturating curve. 
design <- tibble(r=c(1.8, 2, 2.4, 2.8), w=round(25/r), l=round(50/r), n=w*l, sp=c( 6, 4, 2, 1))
plot(design$sp, design$n, ylim=c(0, max(design$n)))
abline(lm(design$n~design$sp))

test <- sapply(c(1,2,4,6), function(i) {
  A <- getabund(i, design$n[design$sp==i])
  R <- rep(i, length(A))
  cbind(R, A)
}
  )
test <- do.call("rbind.data.frame", test)
ggplot(test, aes(A, group=R))+
  geom_bar()+
  facet_grid(~R)

#Saturating curve. define the curves based on the species abundances of the richest community. I'm using 392, a somewhat arbitrary number that seems fine for the most abundant community. don't want to plant too much more than that. 
A <- getabund(6, 392) %>% table()
R <- c(1,2,4,6)
N <- sapply(R, function(i) sum(A[1:i]))
tmp <- sapply(1:length(N), function(i) {
  A <- getabund(R[i], N[i])
  R <- rep(i, length(A))
  cbind(R, A)
  })
tmp <- do.call("rbind.data.frame", tmp)
tmp$A <- as.factor(tmp$A)
tmp$A <- factor(tmp$A, levels = rev(levels(tmp$A)))

ggplot(tmp, aes(R, group=A, fill=as.factor(A)))+
  scale_fill_manual(values = pal)+
  geom_bar()
plot(c(1,2,4,6), N, type='o')
plot(1:4, N, type='o')

#redesigning. Probably the best I'll do.
N # 115 198 319 392 = optimum number to keep constant species abundaces across treatments
#design <- tibble(r=c(1.75, 1.9, 2.28, 3), w=floor(25/r), l=floor(50/r), n=w*l, sp=c( 6, 4, 2, 1)) #previous design using the wrong log-dist
design <- tibble(r=c(1.75, 1.92, 2.5, 3.1), w=floor(25/r), l=floor(50/r), n=w*l, sp=c( 6, 4, 2, 1)) 
design
tmp <- sapply(c(1,2,4,6), function(i) {
  A <- getabund(i, design$n[design$sp==i])
  R <- rep(i, length(A))
  cbind(R, A)
}
)
tmp <- do.call("rbind.data.frame", tmp)
tmp$R <- as.factor(tmp$R)
tmp$A <- as.factor(tmp$A)
tmp$A <- factor(tmp$A, levels = rev(levels(tmp$A)))
ggplot(tmp, aes(R, group=A, fill=A))+
  geom_bar()

#getting abundances for all 4 treatments: sub/det, sub/rand, add/det, add/rand
getabund4 <- function(dist, rand=F){
  #load function
  getabund <- function(r, n){ 
    #r=richness, n=total individuals, mu=mean of log distribution
    
    #get counts of each species
    n1 = n+max(r) #buffering the number of species so there aren't errors due to rounding
    
    #get distribution of each species. comes from a log distribution (mu=6, sd=.4) and 6 species are pulled. replicated 1000 times and the mean values taken.
    N <- 6
    sims <- 1000
    m2 <- matrix(NA, nrow = sims, ncol = N)
    for(i in 1:sims){
      y2 <- rlnorm(1:N, meanlog = 1, sdlog = .5) %>% 
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
  design <- tibble(r=dist, w=floor(25/r), l=floor(50/r), n=w*l, sp=c( 1,2,4,6))
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

#treatments: sub/det, sub/rand, add/det, add/rand
#interplant spacings (cm)
#default distances
D.add <- rev(c(1.75, 1.92, 2.5, 3.1))
#D.add <- rev(c(1.75, 1.9, 2.28, 3)) #old distances with the wrong log-dist
D.sub <- 1.75

A_D <- data.frame(getabund4(D.add, F), dens="additive", rand="deterministic")
A_S <- data.frame(getabund4(D.add, T), dens="additive", rand="stochastic")
S_D <- data.frame(getabund4(D.sub, F), dens="substitutive", rand="deterministic")
S_S <- data.frame(getabund4(D.sub, T), dens="substitutive", rand="stochastic")

plotdesign <- rbind(A_D, A_S, S_D, S_S)

ggplot(plotdesign, aes(R, group=A, fill=A))+
  geom_bar()+
  labs(fill="species", x="richness", y="# individuals")+
  scale_fill_manual(values=pal) +
  theme(strip.background = element_rect(fill="white",color="black"))+
  facet_wrap(~rand+dens)

#exporting the dimensions of the experiment to csv
design <- tibble(d=rep(c("add","sub"), each=4), r=c(1.75, 1.92, 2.5, 3.1, 1.75, 1.75, 1.75, 1.75), w=floor(25/r), l=floor(50/r), n=w*l, sp=c( 6, 4, 2, 1, 6, 4, 2, 1))
design <- design %>% 
  group_by(d) %>% 
  arrange(sp, .by_group=T)
write.csv(design, "simulation/outputs/design.csv", row.names = F)
################################
#plotting

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

