#run the whole script `all_sim_functions2.Rmd` to load data and functions

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
Random=T
density="add"
data <- prepdata(density, Random)
#data$A  %>% group_by(R) %>% table()

S2 <- s2(1, data, B, design)

############################################################
#An attempt to generalize
############################################################
#load packaged simulation function ex/s2

#simulation and results functions ("f.sim" below)
results2 <- function(f.sim=s2, nreps=10, density="sub", Random=F){
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
SD <-  results2(s2, nreps=10, "sub", F) #substitutive, deterministic
SR <-  results2(s2, nreps=10, "sub", T) #substitutive, random
AD <-  results2(s2, nreps=10, "add", F) #additive, deterministic
AR <-  results2(s2, nreps=10, "add", T) #additive, random

#View(SD$response.summary %>% filter(rich==1))

########################################################################
########################################################################
#Plot it!
#summary plots
#deterministic
plotpI <- function(data, main=""){
  results=data$response.summary
  results1 <- results %>% 
    filter(species=="tot")
  ggplot(results1, aes(rich, pI, group=rep)) +
    geom_point()+
    geom_line()+
    ggtitle(main)
}
resultslist <- list(SD, SR, AD, AR)
p1 <- plotpI(resultslist[[1]], "substitutive, deterministic")
p2 <- plotpI(resultslist[[2]], "substitutive, random")
p3 <- plotpI(resultslist[[3]], "additive, deterministic")
p4 <- plotpI(resultslist[[4]], "additive, random")

multiplot(p1, p2, p3, p4, cols = 2)
#results over time
results=SD$responses
results=SR$responses
results=AD$responses
results=AR$responses

ggplot(filter(results, species=="tot"), aes(time, pI, group=rep)) +
  geom_point()+
  geom_line() +
  facet_wrap(~rich)
ggplot(filter(results, species %in% c(1, 2)), aes(time, pI, group=rep, color=species)) +
  geom_point()+
  geom_line() +
  facet_wrap(~rich)
