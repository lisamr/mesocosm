#generalizing script so I can get response variables out for lots of replicates
#*Run Rmd script with all the functions*

#load data
#competency values
CV <- read.csv("simulation/outputs/competencyvalues.csv")
#species abundances
A <- getabund2(n=100, mu=1)

#run simulation
s1 <- function(a, x){
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
A <- getabund2(n=100, mu=1)

#simulation function ("f.sim" below)
s1 <- function(a, x){
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

#functions to bind results
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

results <- function(A, n.reps, f.sim){

  #richness levels for each community
  R <- apply(A, 2, function(x) length(unique(x)))
  
  #1. response
  bres <- list(NULL)
  ident <- data.frame(dens="sub", assemb="det", rich=R, rep=i)
  for(i in 1:nreps){
    res1 <- lapply(1:4, function(x) f.sim(A, x))
    bres[[i]] <- bind.response(res1)
  }
  bres <- plyr::ldply(bres, data.frame)
  
  #2. response summary
  bres.sum <- list(NULL)
  for(i in 1:nreps){
    res1 <- lapply(1:4, function(x) f.sim(A, x))
    ident <- data.frame(rep=i, dens="sub", assemb="det", rich=R)
    bres.sum[[i]] <- bind.respsum(res1)
  }
  bres.sum <- plyr::ldply(bres.sum, data.frame)
  
  list("responses"=bres, "response.sum"=bres.sum)
}

results1 <- results(A = A, nreps=5, f.sim = s1)
results1$responses %>% summary()
