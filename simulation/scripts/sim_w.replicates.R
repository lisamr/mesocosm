#generalizing script so I can get response variables out for lots of replicates
#*Run Rmd script with all the functions*

#load data
#competency values
CV <- read.csv("simulation/outputs/competencyvalues.csv")
#species abundances
A <- getabund2(n=100, mu=1)

#run simulation
s1 <- function(a, x){
  #a=abundance dataframe; x=col in the dataframe
  t <- 15
  h <- format.comp(CV, a[,x], t)
  HP <- HPraster(h, 10, 10, 10)
  sim <- simulate(h, HP, t)
  #animate(HP, sim)

#get response variables
RV1 <- response(sim, HP)
RV2 <- response.sum(RV1)

#report
return(list("format.comp"=h, "HPraster"=HP, "simulate"=sim, "response"=RV1, "respons.sum"=RV2))
}

#get results for each of the 4 richness levels
res1 <- lapply(1:4, function(x) s1(A, x))

#assemble for binding results into 2 tables, instead of lists
R <- apply(A, 2, function(x) length(unique(x)))
ident <- data.frame(rep=1, dens="sub", assemb="det", rich=R)

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

bind.response(res1) %>% head()
bind.respsum(res1) %>% head()

#Try to do it again with more than 1 replicate
R <- apply(A, 2, function(x) length(unique(x)))

#1. response
bres <- list(NULL)
n.reps=4
for(i in 1:n.reps){
  res1 <- lapply(1:4, function(x) s1(A, x))
  ident <- data.frame(rep=i, dens="sub", assemb="det", rich=R)
  bres[[i]] <- bind.response(res1)
}
bres <- plyr::ldply(bres, data.frame)
head(bres)
summary(bres)

#2. response summary
bres.sum <- list(NULL)
for(i in 1:n.reps){
  res1 <- lapply(1:4, function(x) s1(A, x))
  ident <- data.frame(rep=i, dens="sub", assemb="det", rich=R)
  bres.sum[[i]] <- bind.respsum(res1)
}
bres.sum <- plyr::ldply(bres.sum, data.frame)

############################################################
#An attempt to generalize
############################################################
#competency values
CV <- read.csv("simulation/outputs/competencyvalues.csv")
A <- getabund2(n=100, mu=1)

#simulation function ("f.sim" below)
s1 <- function(a, x){
  #a=abundance dataframe; x=col in the dataframe
  t <- 15
  h <- format.comp(CV, a[,x], t)
  HP <- HPraster(h, 10, 10, 10)
  sim <- simulate(h, HP, t)
  #animate(HP, sim)
  
  #get response variables
  RV1 <- response(sim, HP)
  RV2 <- response.sum(RV1)
  
  #report
  return(list("format.comp"=h, "HPraster"=HP, "simulate"=sim, "response"=RV1, "respons.sum"=RV2))
}

results <- function(A, f.sim){
  nreps=10
  R <- apply(A, 2, function(x) length(unique(x)))
  
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
  for(i in 1:nreps){
    res1 <- lapply(1:4, function(x) f.sim(A, x))
    ident <- data.frame(rep=i, dens="sub", assemb="det", rich=R)
    bres[[i]] <- bind.response(res1)
  }
  bres <- plyr::ldply(bres, data.frame)
  #2. response summary
  bres2 <- list(NULL)
  for(i in 1:nreps){
    res1 <- lapply(1:4, function(x) f.sim(A, x))
    ident <- data.frame(rep=i, dens="sub", assemb="det", rich=R)
    bres2[[i]] <- bind.respsum(res1)
  }
  bres2 <- plyr::ldply(bres2, data.frame)
  
  return(list("responses"= bres, "response.summary"= bres2))
}

#get results
R1 <- results(A, s1)
Rt1 <- R1$responses
Rsum1 <- R1$response.summary

#summary plots
head(Rsum1)
Rsum2 <- Rsum1 %>% 
  filter(species=="tot")
head(Rsum2)

ggplot(Rsum2, aes(rich, pI, group=rep)) +
  geom_point()+
  geom_line()

ggplot(Rsum2, aes(rich, max.dI, group=rep)) +
  geom_point()+
  geom_line()

#responses over time
head(Rt1)
ggplot(filter(Rt1, species=="tot"), aes(time, pI, group=rep)) +
  geom_point()+
  geom_line() +
  facet_wrap(~rich)

ggplot(filter(Rt1, species=="tot"), aes(time, dI, group=rep)) +
  geom_point()+
  geom_line() +
  facet_wrap(~rich)

ggplot(filter(Rt1, species=="tot"), aes(time, n.E, group=rep)) +
  geom_point()+
  geom_line() +
  facet_wrap(~rich)
