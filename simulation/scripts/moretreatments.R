#adding 3 more treatments: replicating substitutive design at each of the additive densities

design2 <- read.csv("simulation/outputs/design2.csv")
design2$r[design2$r==2.28] <- 2.29 #for some reason R sucks and hates 2.28. changing it slightly.

design.128 <- design2 %>% filter(d=="add"| d=="sub" & n==128)
design.210 <- design2 %>% filter(d=="add"| d=="sub" & n==210)
design.338 <- design2 %>% filter(d=="add"| d=="sub" & n==338)
design.392 <- design2 %>% filter(d=="add"| d=="sub" & n==392)

prepdata2 <- function(dens, rand, design){
  #dens="add" or "sub", rand=F or T
  D <- design %>% filter(d==dens) %>% pull(r)#(planting distances)
  A <- getabund4(D, rand, 1) 
  list(A=A, D=D, dens=dens, rand=rand)
}

#load data. include density and random treatment
Random=F
density="sub"
data <- prepdata2(density, Random, design.392)

S2 <- s2(2, data, B, design.392) #can ignore warnings. Has to do with binding results together
animate(S2$HPraster, S2$simulate)


#simulation and results functions ("f.sim" below)
results3 <- function(f.sim=s2,B, nreps=10, density="sub", Random=F, design){
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
    data <- prepdata2(density, Random, design)
    res1 <- lapply(1:4, function(R) s2(R, data, B, design))
    ident <- data.frame(rep=i, dens=density, rand=Random, rich=R)
    bres[[i]] <- bind.response(res1)
    bres2[[i]] <- bind.respsum(res1)
  }
  bres <- plyr::ldply(bres, data.frame)
  bres2 <- plyr::ldply(bres2, data.frame)
  
  return(list("responses"= bres, "response.summary"= bres2))
}
fb2 <- function(Beta, n=10, design){
  SD <-  results3(s2,B=Beta, nreps=n, "sub", F, design) #substitutive, deterministic
  SR <-  results3(s2, B=Beta, nreps=n, "sub", T, design) #substitutive, random
  response <- rbind(SD$responses, SR$responses)
 
  return(response)
}
results3(s2, B, 2, "sub", F, design.392)

#bind all the treatments together. resnormal is what has been orignally run and the res with numbers are those with additional substitutive treatments at different densities.
resnormal <- fb(B, 10)[[1]] %>% mutate(densN=NA)
res128 <- fb2(B, 10, design.128) %>% mutate(densN="sub128")
res210 <- fb2(B, 10, design.210) %>% mutate(densN="sub210")
res338 <- fb2(B, 10, design.338) %>% mutate(densN="sub338")
res392 <- fb2(B, 10, design.392) %>% mutate(densN="sub392")
l <- list(resnormal, res128, res210, res338, res392)
resbind <- bind_rows(l)
head(resbind)

results <- resbind %>% 
  mutate("comp"=Rename(c(1, .3, .2, .1, 0, 0), c(1:6), species)) %>%
  mutate("comp_Abund"=comp*n) %>% 
  filter(time==max(time)) %>%
  group_by(rep, dens, rand,rich, densN) %>% 
  summarise(N = n[species=="tot"],
            pI=pI[species=="tot"],
            n.I=n.I[species=="tot"],
            comp_Abund=sum(comp_Abund, na.rm = T),
            comp_Abund.N=sum(comp_Abund, na.rm = T)/N) 

ggplot(filter(results, is.na(densN)), aes(rich, pI, color=dens, shape=rand))+
  geom_point()+
  geom_line(aes(group=rep))+
  #geom_smooth(method = "lm", formula = y ~poly(x,3), color="black", lwd=.5)+
  facet_wrap(~rand+dens)+
  labs(color="density", shape="stochastic?", y="proportion infected", x="richness")
ggplot(filter(results, is.na(densN)), aes(rich, n.I, color=dens, shape=rand))+
  geom_point()+
  geom_line(aes(group=rep))+
  #geom_smooth(method = "lm", formula = y ~poly(x,3), color="black", lwd=.5)+
  facet_wrap(~rand+dens)+
  labs(color="density", shape="stochastic?", y="# infected", x="richness")
ggplot(filter(results, is.na(densN)), aes(rich, comp_Abund, color=dens, shape=rand))+
  geom_point()+
  geom_line(aes(group=rep))+
  #geom_smooth(method = "lm", formula = y ~poly(x,2), color="black")+
  facet_wrap(~rand+dens)+
  labs(color="density", shape="stochastic?", y="sum(competency x abundance)", x="richness")
ggplot(filter(results, is.na(densN)), aes(rich, comp_Abund.N, color=dens, shape=rand))+
  geom_point()+
  geom_line(aes(group=rep))+
  #geom_smooth(method = "lm", formula = y ~poly(x,2), color="black")+
  facet_wrap(~rand+dens)+
  labs(color="density", shape="stochastic?", y="sum(competency x abundance)/N", x="richness")
ggplot(filter(results, is.na(densN)), aes(comp_Abund.N, pI, color=as.factor(rich)))+
  geom_point(aes(shape=rand))+
  #geom_line(aes(group=rep))+
  #geom_smooth(aes(group=rand),color="black")+
  geom_smooth(color="black")+
  facet_wrap(~dens)+
  labs(color="density", shape="stochastic?", x="sum(competency x abundance)/N", y="proportion infected")

#relative community comp vs proportion infected, grouped by density.
head(results)
ggplot(filter(results, !N==300), aes(comp_Abund.N, pI, group=N, color=as.factor(N)))+
  geom_point()+
  geom_smooth(se = F, lwd=.5)

#PI generally increases with density, to a point. It increases, then plateues or goes back down at higher densities. how do you show that...
#bin comp_Abund.N into groups
results$comp_Abund.Ncut <-  cut_number(results$comp_Abund.N, 4)

ggplot(filter(results, !N==300), aes(N, pI, group=comp_Abund.Ncut, color=(comp_Abund.Ncut)))+
  geom_point()+
  geom_smooth(se=F, method='lm', formula = y~poly(x,2))
ggplot(filter(results, !N==300), aes(N, pI, color=comp_Abund.N))+
  geom_point()

#try fitting that plot to a parameter surface to see how density and community competency interact to affect incidence
require(plot3D)

# linear fit
fit <- lm(pI ~ N*comp_Abund.N, data = results)
fit1 <- lm(pI ~ N+comp_Abund.N, data = results)
fit2 <- lm(pI ~ poly(comp_Abund.N, 3) + N, data = results)
fit3 <- lm(pI ~ poly(comp_Abund.N, 3) + poly(N,2), data = results)
fit4 <- lm(pI ~ poly(comp_Abund.N, 3) * poly(N,2), data = results)
fit5 <- lm(pI ~ poly(comp_Abund.N, 3) * poly(N,3), data = results)
AIC(fit, fit1, fit2, fit3, fit4, fit5)

#fit 5 is the best
plot(fit5)
hist(resid(fit5))
plot(fitted.values(fit5) ~ results$pI) + abline(0,1) #overestimates pI at higher values by a little.

# predict on x-y grid, for surface
N.pred <- seq(min(results$N), max(results$N), length.out = 30 )
comp.pred <- seq(min(results$comp_Abund.N), max(results$comp_Abund.N), length.out = 30 )
xy <- expand.grid(N = N.pred, comp_Abund.N = comp.pred)

pI.pred <- matrix (nrow = 30, ncol = 30, 
                    data = predict(fit5, newdata = data.frame(xy), interval = "prediction"))

# predicted z-values, fitted points for droplines to surface
fitpoints <- predict(fit5) 

scatter3D(z = results$pI, x = results$N, y = results$comp_Abund.N,
          theta = 20, phi = 20, ticktype = "detailed",
          pch=16, cex=.8, col = viridis(12, direction=1),
          xlab = "N", ylab = "comp_Abund.N", zlab = "pI", clab = "pI", 
          surf = list(x = N.pred, y = comp.pred, z = pI.pred, 
                      facets = NA, alpha=.7, fit = fitpoints),
          colkey = list(length = 0.8, width = 0.4),            
          main = "Proportion infected")


# TRY: 
# require(plot3Drgl)
# plotrlg()
  
  
  
  
  