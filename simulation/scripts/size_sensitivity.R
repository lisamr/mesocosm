#sensitivity analysis
#need to see how number of individuals within each mesocosm affects statistical power. need to think about this for tractibility and because infection distances might vary in my system. easiest way to do this is to make the lattices smaller rather than fussing with infection distances. 
library(dplyr)
rm(list=ls())
#load functions from 'all_sim_functions2.Rmd'
#load this function too
plot_pipeline <- function(sim, plot1){
  
  #analyze data
  dat <- sim$responses
  
  dat <- dat %>% 
    group_by(rep, dens, rand, rich) %>% 
    mutate(rel.n=n/n[species == "tot"]) %>% 
    mutate("comp"=Rename(c(1, .3, .2, .1, 0, 0), c(1:6), species)) %>%
    mutate("comp_Abund"=comp*n, "comp_relAbund"=comp*rel.n) 
  
  #data at time final
  dat30 <- dat %>% 
    filter(time==30, species=="tot") %>%  #look at total infections at time final
    mutate(trt = interaction(dens, rand),
           rep2 = interaction(trt, rep), #rep ID in each trt group
           rep3 = case_when(rand==F~"det", T~as.character(rep2)),
           comp=NULL, comp_Abund=NULL, comp_relAbund=NULL)
  
  #get summary values for community competency and merge
  sumstats <- dat %>% 
    filter(time==max(time)) %>% 
    group_by(rep, dens, rand, rich) %>% 
    summarise(comp_Abund=sum(comp_Abund, na.rm = T),
              comp_Abund.N=sum(comp_Abund, na.rm = T)/n[species=="tot"]) 
  dat30 <- left_join(dat30, sumstats, by = c("rep", "dens", "rand", "rich"))
  dat30 <- ungroup(dat30)
  
  ##################################################
  #visualize results
  ##################################################
  #1. disease vs richness
  plotA=ggplot(dat30, aes(rich, n.I, group=interaction(rand, dens)))+
    geom_point()+
    geom_line(aes(group=rep))+
    facet_wrap(~rand+dens)+
    background_grid(major = "xy", minor = "none") +
    labs(x="richness", y="total infected")
  plotB=ggplot(dat30, aes(rich, pI, group=interaction(rand, dens)))+
    geom_point()+
    geom_line(aes(group=rep))+
    facet_wrap(~rand+dens)+
    background_grid(major = "xy", minor = "none") +
    labs(x="richness", y="proportion infected")
  p1 <- plot_grid(plotA, plotB, labels = c("A", "B"))
  
  #2. community competency vs disease
  p2 <- ggplot(dat30, aes(comp_Abund, n.I))+
    geom_point()+
    background_grid(major = "xy", minor = "none") +
    labs(x="community competency (p)", y="total infected")+
    geom_smooth(method='lm', formula = y ~ poly(x, 3))
  
  if(plot1==T){p1}else{p2} #choose which one to plot
}

#design communities with different sizes. go from a square 1:1 to exisiting 2:1. 
fdesign <- function(L) {
  d <- tibble(
  d=rep(c("add","sub"), each=4), 
  r=c(1.75, 1.92, 2.39, 3, 1.75, 1.75, 1.75, 1.75), 
  w=floor(25/r), 
  l=floor(L/r), 
  n=w*l, 
  sp=c( 6, 4, 2, 1, 6, 4, 2, 1),
  L=L) %>% 
    group_by(d) %>% 
    arrange(sp, .by_group=T)
  d$r[d$r==2.39] <- 2.4 #for some reason R sucks and hates 2.39.
  return(as.data.frame(d))
}

#examine replicates over a range of differently sized containers
Seqs <- seq(10, 50, by=5)
designs <- lapply(Seqs, fdesign)

#test out 1 new design
i <- 1
des=designs[[i]]
data <- prepdata2(dens = 'sub', rand = F, L=unique(des$L), des)
#run simulation. works :)
S2 <- s2(Rich=4, data, B, des, t=30)
animate(S2$HPraster, S2$simulate)

#simulate 30 time steps, 10 reps per treatment, 4 treatments.
sims <- fb3(Beta = B, t = 30, L =des$L[1], n =2, design = des)
plot_pipeline(sims, F)

#loop over lots
plotlist1 <- list(NULL)
plotlist2 <- list(NULL)
for(i in 1:length(designs)){
  sims <- fb3(Beta = B, t = 30, L =designs[[i]]$L[1], n =10, design = designs[[i]])
  plotlist1[[i]] <- plot_pipeline(sims, T)
  plotlist2[[i]] <- plot_pipeline(sims, F)
}

#save plots
pdf("simulation/images/sensitivity_panelplots.pdf",onefile = TRUE, width = 7, height = 5)
sapply(1:length(designs), function(x) print(plotlist1[[x]]))
dev.off()
pdf("simulation/images/sensitivity_comcomp.pdf",onefile = TRUE, width = 7, height = 5)
sapply(1:length(designs), function(x) print(plotlist2[[x]]))
dev.off()