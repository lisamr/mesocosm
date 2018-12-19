#run the whole script `all_sim_functions2.Rmd` to load data and functions
B1 <- read.csv("simulation/outputs/rates111111.csv")
B2 <- read.csv("simulation/outputs/rates1.7.5.4.2.1.csv")
B3 <- read.csv("simulation/outputs/rates1.3.2.100.csv")

#simulation function that packages together simulation functions and responses
s2 <- function(Rich, A, beta, design){
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
  RV1 <- response(sim, HP)
  foi <- FOI(sim, HP)
  RV2 <- response.sum(RV1, foi)
  
  #report
  return(list("format.comp"=h, "HPraster"=HP, "simulate"=sim, "resp1"=RV1, "resp2"=RV2))
}

#load data. include density and random treatment
Random=F
density="add"
data <- prepdata(density, Random)
#data$A  %>% group_by(R) %>% table()

S2 <- s2(4, data, B, design) #can ignore warnings. Has to do with binding results together
animate(S2$HPraster, S2$simulate)

############################################################
#An attempt to generalize
############################################################
#load packaged simulation function ex/s2

#simulation and results functions ("f.sim" below)
results2 <- function(f.sim=s2,B, nreps=10, density="sub", Random=F){
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
    data <- prepdata(density, Random)
    res1 <- lapply(1:4, function(R) s2(R, data, B, design))
    ident <- data.frame(rep=i, dens=density, rand=Random, rich=R)
    bres[[i]] <- bind.response(res1)
    bres2[[i]] <- bind.respsum(res1)
  }
  bres <- plyr::ldply(bres, data.frame)
  bres2 <- plyr::ldply(bres2, data.frame)
  
  return(list("responses"= bres, "response.summary"= bres2))
}
fb <- function(Beta){
  SD <-  results2(s2,B=Beta, nreps=10, "sub", F) #substitutive, deterministic
  SR <-  results2(s2, B=Beta, nreps=10, "sub", T) #substitutive, random
  AD <-  results2(s2, B=Beta, nreps=10, "add", F) #additive, deterministic
  AR <-  results2(s2, B=Beta, nreps=10, "add", T) #additive, random
  rbind(SD$responses, SR$responses, AD$responses, AR$responses)
}
results2(s2, B3, 2, "sub", F)

#results with different "beta" values 
#resB1 <- fb(B1) #high competency across species
resB2 <- fb(B2) #medium decline in competency across species
resB3 <- fb(B3) #big decine in competency across species

########################################################################plot results
#######################################################################
#resB2 <- read.csv("simulation/outputs/resultsB2.csv")
#resB3 <- read.csv("simulation/outputs/resultsB3.csv")

head(resB2)
#exposure over time
ggplot(resB3 %>% filter(species==1), aes(time, n.E, group=interaction(rep, rich), color=rich))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens+rand)

#adding in community competency (sum(competency x n.I)). competency will be defined as the relative magnitude of transmissions across species. 

resB2.1 <- resB2 %>% 
  group_by(rep, dens, rand, rich) %>% 
  mutate(rel.n=n/n[species == "tot"]) %>% 
  mutate("comp"=Rename(c(1, .7, .5, .4, .2, .1), c(1:6), species)) %>%
  mutate("comp_Abund"=comp*n, "comp_relAbund"=comp*rel.n, "comp_nI"=comp*n.I) 

resB3.1 <- resB3 %>% 
group_by(rep, dens, rand, rich) %>% 
  mutate(rel.n=n/n[species == "tot"]) %>% 
  mutate("comp"=Rename(c(1, .3, .2, .1, 0, 0), c(1:6), species)) %>%
  mutate("comp_Abund"=comp*n, "comp_relAbund"=comp*rel.n, "comp_nI"=comp*n.I) 

head(resB2.1)
#variables to look at: pI, comp_Abund, comp_relAbund, comp_nI
resplot <- resB2.1 %>% 
  filter(time==max(time)) %>% 
  group_by(rep, dens, rand, rich) %>% 
  summarise(variable=sum(comp_nI, na.rm = T))
ggplot(resplot, aes(rich, variable, group=rep))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens+rand)
ggplot(resB2.1 %>% filter(species=="tot", time==max(time)), aes(rich, pI, group=interaction(rep)))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens+rand)

resplot <- resB3.1 %>% 
  filter(time==max(time)) %>% 
  group_by(rep, dens, rand, rich) %>% 
  summarise(pI=pI[species=="tot"],
            comp_Abund=sum(comp_Abund, na.rm = T),
            comp_relAbund=sum(comp_relAbund, na.rm = T),
            comp_nI=sum(comp_nI, na.rm = T)
            ) 
ggplot(resplot, aes(rich, pI, group=rep, color=rand, shape=dens))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens+rand)+
  labs(x = "Richness", y="Total proportion infected")
ggplot(resplot, aes(rich, comp_Abund, group=rep, color=rand, shape=dens))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens+rand)+
  labs(x = "Richness", y="Sum(competency.i x abundance.i)")
ggplot(resplot, aes(comp_Abund, pI, group=rep, color=rand, shape=dens))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens+rand) +
  labs(x = "Sum(competency.i x abundance.i)", y="Total proportion infected")
ggplot(resplot, aes(comp_Abund, pI, group=rep, color=rand, shape=dens))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens) +
  labs(x = "Sum(competency.i x abundance.i)", y="Total proportion infected")

#how does exposure change with abundance?
#the following code doesn't work. I need total # exposed, but the code double counts the exposed plants. 
resplot2 <- resB3.1 %>% 
  group_by(rep, species, dens, rand, rich) %>% 
  summarise(abund = n[time == 1],
            n.Ex = sum(n.E),
            p.Ex = n.Ex/abund)
ggplot(filter(resplot2, species == "tot"), aes(abund, p.Ex, group=rep, color=rich, shape=dens))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens + rand) +
  labs(x = "# plants", y="proportion exposed")
ggplot(filter(resplot2, species == "tot"), aes(abund, n.Ex, group=rep, color=rich, shape=dens))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens + rand) +
  labs(x = "# plants", y="# exposed")
ggplot(filter(resplot2, species == 1), aes(abund, n.Ex, group=rep, color=rich, shape=dens))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens + rand) +
  labs(x = "# plants", y="# exposed")

ggplot(resB3.1 %>% filter(species=="tot", time==max(time)), aes(rich, pI, group=interaction(rep)))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens+rand)

#write.csv(resB2.1, "simulation/outputs/resultsB2.csv")
#write.csv(resB3.1, "simulation/outputs/resultsB3.csv")

###
#Make competency reflect Bii, which changes with distance. Values from monospecific simulations. Bii=average exposure/infection rate. See "betaii.R" for how those values were generated.
