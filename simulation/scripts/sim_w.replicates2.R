#run the whole script `all_sim_functions2.Rmd` to load data and functions
B1 <- read.csv("simulation/outputs/rates111111.csv")
B2 <- read.csv("simulation/outputs/rates1.7.5.4.2.1.csv")
B3 <- read.csv("simulation/outputs/rates1.3.2.100.csv")

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

#load data. include density and random treatment
Random=F
density="add"
data <- prepdata(density, Random)
#data$A  %>% group_by(R) %>% table()

S2 <- s2(1, data, B, design, 40) #can ignore warnings. Has to do with binding results together
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
fb <- function(Beta, n=10){
  SD <-  results2(s2,B=Beta, nreps=n, "sub", F) #substitutive, deterministic
  SR <-  results2(s2, B=Beta, nreps=n, "sub", T) #substitutive, random
  AD <-  results2(s2, B=Beta, nreps=n, "add", F) #additive, deterministic
  AR <-  results2(s2, B=Beta, nreps=n, "add", T) #additive, random
  response <- rbind(SD$responses, SR$responses, AD$responses, AR$responses)
  response.summary <- rbind(SD$response.summary, SR$response.summary, AD$response.summary, AR$response.summary)
  
  return(list(responses=response, response.summary=response.summary))
}
results2(s2, B3, 2, "sub", F)

#results with different "beta" values 
#resB1 <- fb(B1) #high competency across species
resB2 <- fb(B2, 10) #medium decline in competency across species
resB3 <- fb(B3, 10) #big decine in competency across species

########################################################################plot results
#######################################################################
#resB2 <- read.csv("simulation/outputs/resultsB2.csv")
#resB3 <- read.csv("simulation/outputs/resultsB3.csv")

#adding in community competency (sum(competency x n.I)). competency will be defined as the relative magnitude of transmissions across species. 

#resB2.1 <- resB2[[1]] %>% group_by(rep, dens, rand, rich) %>% mutate(rel.n=n/n[species == "tot"]) %>% mutate("comp"=Rename(c(1, .7, .5, .4, .2, .1), c(1:6), species)) %>% mutate("comp_Abund"=comp*n, "comp_relAbund"=comp*rel.n, "comp_nI"=comp*n.I) 
ggplot <- function(...) ggplot2::ggplot(...) + scale_fill_manual(values=pal) + scale_color_manual(values=pal) 
ggplot <- function(...) ggplot2::ggplot(...) 


resB3.1 <- resB3[[1]] %>% 
group_by(rep, dens, rand, rich) %>% 
  mutate(rel.n=n/n[species == "tot"]) %>% 
  mutate("comp"=Rename(c(1, .3, .2, .1, 0, 0), c(1:6), species)) %>%
  mutate("comp_Abund"=comp*n, "comp_relAbund"=comp*rel.n, "comp_nI"=comp*n.I) 

#plot (just B3 for now)
resplot <- resB3.1 %>% 
  filter(time==max(time)) %>% 
  group_by(rep, dens, rand, rich) %>% 
  summarise(pI=pI[species=="tot"],
            comp_Abund=sum(comp_Abund, na.rm = T),
            comp_Abund.N=sum(comp_Abund, na.rm = T)/n[species=="tot"],
            comp_nI=sum(comp_nI, na.rm = T)            ) 
ggplot(resplot, aes(rich, pI, color=dens, shape=rand))+
  geom_point()+
  geom_line(aes(group=rep))+
  geom_smooth(method = "lm", formula = y ~poly(x,3), color="black", lwd=.5)+
  facet_wrap(~rand+dens)+
  labs(color="density", shape="stochastic?", y="proportion infected", x="richness")
ggplot(resplot, aes(rich, comp_Abund, color=dens, shape=rand))+
  geom_point()+
  geom_line(aes(group=rep))+
  #geom_smooth(method = "lm", formula = y ~poly(x,2), color="black")+
  facet_wrap(~rand+dens)+
  labs(color="density", shape="stochastic?", y="sum(competency x abundance)", x="richness")
ggplot(resplot, aes(rich, comp_Abund.N, color=dens, shape=rand))+
  geom_point()+
  geom_line(aes(group=rep))+
  #geom_smooth(method = "lm", formula = y ~poly(x,2), color="black")+
  facet_wrap(~rand+dens)+
  labs(color="density", shape="stochastic?", y="sum(competency x abundance)", x="richness")
ggplot(filter(resplot), aes(comp_Abund, pI, color=dens))+
  geom_point(aes(shape=rand))+
  #geom_line(aes(group=rep))+
  #geom_smooth(aes(group=rand),color="black")+
  geom_smooth(color="black")+
  facet_wrap(~dens)+
  labs(color="density", shape="stochastic?", x="sum(competency x abundance)", y="proportion infected")

#how does exposure change with abundance?
resB3.2 <- resB3$response.summary
resB3.2$pE <- resB3.2$n.E / resB3.2$n
head(resB3.2)
ggplot(filter(resB3.2, species == "tot"), aes(n, pE, group=rep, color=rich, shape=dens))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens + rand) 
ggplot(filter(resB3.2, species == "tot"), aes(n, pI, group=rep, color=rich, shape=dens))+
  geom_point()+
  facet_wrap(~rich)
ggplot(filter(resB3.2, species == "tot"), aes(rich, pE, group=rep, color=rich, shape=dens))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens + rand) 
ggplot(filter(resB3.2, species == "tot"), aes(pE, pI, group=rep, color=rich, shape=dens))+
  geom_point()
  #geom_line()
  #facet_wrap(~dens ) 

#richness vs. PI and richness vs FOI
ggplot(filter(resB3.2, species == "tot"), aes(rich, pI, group=rep, color=rich, shape=dens))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens + rand) 
ggplot(filter(resB3.2, species == "tot"), aes(rich, FOI, group=rep, color=rich, shape=dens))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens + rand, scales = "free") 



#get exposure weighted by competency vs comm. comp.
#vector of competency
resplot3 <-  resB3.2 %>% 
  mutate(pI=NULL) %>% 
  left_join(resplot,
            by = c("rep", "dens", "rand", "rich") ) 
mag <- c(1, .3, .2, .1, 0, 0, NA)
resplot3 <- resplot3 %>% 
  mutate(comp.i = Rename(values_to = mag, index_from = as.character(c(1, 2, 3, 4, 5, 6, "tot")), from_column = resplot3$species),
         n.Eweighted = n.E * comp.i) %>% 
  group_by(rep, dens, rand, rich) %>% 
  mutate(n.Eweighted = ifelse(species == "tot", sum(n.Eweighted, na.rm = T), n.Eweighted ),
         p.Eweighted = n.Eweighted/n[species == "tot"])

#exposure drives infection across all treatments
ggplot(filter(resplot3, species == "tot"), aes(rich, comp_Abund, color=rich, shape=dens, group=rep))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens+rand)
ggplot(filter(resplot3, species == "tot"), aes(p.Eweighted, pI, color=rich, shape=dens))+
  geom_point()
ggplot(filter(resplot3, species == "tot"), aes(pE, pI, color=rich, shape=dens))+
  geom_point()

#richness vs percent infected of species1?
resplot4 <- resplot3 %>%  
  group_by(rep, dens, rand, rich) %>%  filter(any(species==1)) %>% 
  summarise(n.1 = n[species==1],
            nI.1 = n.I[species==1],
            nE.1 = n.E[species==1],
            pE.1 = pE[species==1],
            p.Eweighted = p.Eweighted[species==1],
            pI.1 = pI[species==1],
            N = n[species=="tot"]) 

#seems like it's maybe a function of  exposure and density
ggplot(resplot4, aes(rich, pI.1, color=pE.1, group=rep, shape=dens))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens+rand)+
  ylim(0,1)
ggplot(filter(resplot4, rand==F), aes(rich, pI.1, color=pE.1, group=rep, shape=dens))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens+rand)+
  ylim(0,1)
ggplot(filter(resplot4, rand==F), aes(rich, pE.1,group=rep, shape=dens))+
  geom_point()+
  geom_line()+
  facet_wrap(~dens+rand)+
  ylim(0,1)
ggplot(filter(resplot4, rand==F), aes(pE.1, pI.1, color=rich))+
  geom_point()
ggplot(filter(resplot4, rand==F), aes(p.Eweighted, pI.1, color=rich))+
  geom_point()
ggplot(resplot4, aes(N, pI.1, color=pE.1))+
  geom_point()


#as avg community competency goes up, proportion exposed goes up
ggplot(filter(resplot3, species == "tot"), aes(pcomp_Abund, p.Eweighted, color=rich, shape=dens))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~dens ) 
ggplot(filter(resplot3, species == "tot"), aes(pcomp_Abund, pE, color=rich, shape=dens))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~dens ) 
names(resplot3)
#which explains why pI decreses with higher commuity competency
ggplot(filter(resplot3, species == "tot"), aes(rich, pI, color=rich, shape=dens))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~dens + rand ) 


write.csv(resB3$response.summary, "simulation/outputs/B3.respsum.csv", row.names = F)
write.csv(resB3$responses, "simulation/outputs/B3.resp.csv", row.names = F)
#write.csv(resB2, "simulation/outputs/resultsB2.csv")
#write.csv(resplot3, "simulation/outputs/resplot3B3.csv", row.names = F)

###
#Make competency reflect Bii, which changes with distance. Values from monospecific simulations. Bii=average exposure/infection rate. See "betaii.R" for how those values were generated.
