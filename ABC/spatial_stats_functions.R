require(spatstat)
require(sp)
library(tidyverse)

rm(list=ls())

source('ABC/ABC_functions.R')


#generate simulation data using the script 'simple_simulation.R' to test of summary statistics

#setup community grid----
r <- 1.7 #interplanting distance in cm
pinoc <- .1 #percent inoculated
s <- 1 #number of spp
spp <- c("sp1") #species name
tfinal <- 20 #time steps

#simulate communities in a 9.5in square tray in hexagonal grid
grid <- sample_community(spp, pinoc, r) 
agents <- as.data.frame(grid) %>% mutate(ID = row_number())
N <- nrow(agents)

#distance matrix
d <- as.matrix(dist(cbind(agents$x, agents$y), method = 'euclidean', diag = TRUE, upper = TRUE))

head(agents)

#simulate data
D <- f_sim1(.75, 1, .5) #"True" data
Dstar <- f_sim1(.75, .5, .5) #simulated data


#quick comparison
#compare "data" with simulation
tmp <- data.frame(I = colSums(D=="I"),
                  Istar = colSums(Dstar=="I"),
                  S = colSums(D=="S"),
                  Sstar = colSums(Dstar=="S"),
                  C = colSums(D=="C"),
                  Cstar = colSums(Dstar=="C"),
                  t = 1:tfinal)

par(mfrow = c(1,3))
plot(tmp$t, tmp$I, type = 'o', col = 'blue', main = 'I', ylim=c(0,170))
lines(tmp$t, tmp$Istar, type = 'o', col = 'red')
plot(tmp$t, tmp$S, type = 'o', col = 'blue', main = 'S->I', ylim=c(0,170))
lines(tmp$t, tmp$Sstar, type = 'o', col = 'red')
plot(tmp$t, tmp$C, type = 'o', col = 'blue', pch=16, main = 'C->I')
lines(tmp$t, tmp$Cstar, type = 'o', col = 'red', pch=16)
par(mfrow=c(1,1))









#O-ring-----

#Convert to a point pattern
npp <- ppp(x=agents$x, y = agents$y, c(min(agents$x), max(agents$x)),
           c(min(agents$y), max(agents$y)))
marks(npp) <- as.factor(D[,4]) #states at an arbitrary time point, used for example. 


#O-ring statistic: proportional to the pair correlation function, which tells you the density of points across a range of distances from Infecteds to infecteds. Not sure how method changes the results, but 'spar' is a smoothing parameter. 
K12 <- Kcross(npp, 'I', 'I')
g12 <- pcf(K12,  method="d", spar=0.7)
lambda2 <- summary(npp)$marks['I',"intensity"]
Oring <- eval.fv(lambda2*g12)

par(mfrow=c(1,2))
plot(npp)
plot(Oring)
par(mfrow=c(1,1))

#peaks are at the first 6 unique distances. 
plot(Oring, add=T)
R <- sort(unique(round(c(d), 8)))
abline(v=R[1:6], lwd=2, lty=2, col='blue')


#explore effects of methods with the pcf function. Methods will interact with smoothing parameter, so beware. 
g12a <- pcf(K12,  method="a")
g12b <- pcf(K12,  method="b")
g12c <- pcf(K12,  method="c")
g12d <- pcf(K12,  method="d")
par(mfrow=c(2,2))
plot(g12a)
plot(g12b)
plot(g12c)
plot(g12d)

#spar is a smoothing parameter.
g12.9 <- pcf(K12,  method="a", spar=.8)
g12.7 <- pcf(K12,  method="a", spar=.7)
g12.3 <- pcf(K12,  method="a", spar=.3)
g120 <- pcf(K12,  method="a", spar=0)
par(mfrow=c(2,2))
plot(g12.9)
plot(g12.7)
plot(g12.3)
plot(g120)
par(mfrow=c(1,1))

#check out how pcf changes over time. location of peaks stay the same, but the amplitude increases. 

f_Oring2 <- function(ppp, t, Spar=.7){
  marks(npp) <- as.factor(D[,t]) #update time
  K12 <- Kcross(npp, 'I', 'I') #calculate o-ring stat 
  g12 <- pcf(K12,  method="a", spar=Spar)
  lambda2 <- summary(npp)$marks['I',"intensity"]
  Oring <- eval.fv(lambda2*g12)
  df <- data.frame(Oring = Oring$pcf, r = Oring$r, t = t)
  
  return(df)
}

tmpO <- lapply(2:20, function(t) f_Oring2(npp, t)) %>% 
  bind_rows()


#how does the O-ring statistic change over time and distance?
ggplot(tmpO, aes(r, Oring, group = t, color = t)) +
  geom_line() +
  scale_color_viridis_c()

#whats the best distance? could find the peaks, which will change with planting density.
peaks <- function(x){
  which(diff(sign(diff(x)))==-2)+1
}

#find peaks
peakstable <- tmpO %>% 
  group_by(t) %>% 
  summarise(r1 = r[peaks(Oring)[1]],
            r2 = r[peaks(Oring)[2]],
            r3 = r[peaks(Oring)[3]],
            r4 = r[peaks(Oring)[4]],
            r5 = r[peaks(Oring)[5]])
x <- tmpO %>% filter(t==10) %>% pull(Oring)

pk_r <- as.vector(round(peakstable[nrow(peakstable),-1], 3))

tmpO %>% 
  mutate(r = round(r, 3)) %>% 
  filter(r %in% pk_r) %>% 
  filter(!duplicated(.)) %>% 
  ggplot(., aes(t, Oring, color = as.factor(r), group = r)) +
  geom_point()+
  geom_line() +
  scale_color_viridis_d()
#it looks a lot like disease progress curve (sum(I_t)), so not sure if this is any better.


#Ripley's K----

#cumulative number of events across a range of distances. 
fK2 <- function(Marks, correction=c("border", "isotropic", "Ripley", "translate")){
  marks(npp) <- as.factor(Marks)
  Kest(subset(npp, marks == 'I'), correction) 
}

par(mfrow=c(2,2))
fK2(D[,5]) %>% plot
fK2(D[,10]) %>% plot
fK2(D[,15]) %>% plot
fK2(D[,20]) %>% plot
par(mfrow=c(1,1))

#compare K over time
fKtime <- function(t, C){
  tmp <- fK2(D[,t], C)
  d <- data.frame(r=tmp[[1]], k = tmp[[3]], t=t)
  return(d)
}

#how does k change over time and distance?
tmp <- lapply(2:20, function(t) fKtime(t, 'Ripley')) %>% 
  bind_rows()
ggplot(tmp, aes(r, k, group = t, color = t)) +
  geom_line() +
  scale_color_viridis_c()

#whats the best distance?
tmp %>% 
  mutate(r = round(r, 1)) %>% 
  filter(r %in% c(2.5, 3.2, 3.6, 4.5)) %>% 
  filter(!duplicated(.)) %>% 
  ggplot(., aes(t, k, color = r, group = r)) +
  geom_point()+
  geom_line() +
  scale_color_viridis_c()















#candidate summary statistics----



# 1. number of infecteds

I <- apply(D, 2, nI)
Istar <- apply(Dstar, 2, nI)

SSE(I, Istar)
chisq(I, Istar)

sqrt((apply(D, 2, nI) - apply(Dstar, 2, nI))^2) %>% plot


# 2. nearest infected neighbor from Minter 2019
#not going to be good if most infections are from nearest neighbor!

uniq.dist <- sort(unique(round(c(d), 8)))
uniq.dist

nobs <- function(D){
  n.obs <- matrix(NA, nrow = length(uniq.dist), ncol = ncol(D))
  for(i in 1:ncol(D)){
    inf <- which(D[,i]=="I")
    if(length(inf)>0){
      NinfN <- apply(d[inf, inf], 1, function(x) round(min(x[x!=0]), 8))
    }else NinfN <- 0
    n.obs[,i] <- sapply(1:length(uniq.dist), function(a) sum(NinfN==uniq.dist[a]))
  }
  return(n.obs)
}

sqrt(colSums((nobs(D) - nobs(Dstar))^2)) %>% plot

NN <- nobs(D)
NNstar <- nobs(Dstar)







# 3. o-ring statistic


#create point pattern object for dataset
npp <- ppp(x=agents$x, y = agents$y, c(min(agents$x), max(agents$x)),c(min(agents$y), max(agents$y)))


#arbitrarily choose r=2
O <- sapply(1:20, function(t) f_Oring(npp, D, t, 2))
Ostar <- sapply(1:20, function(t) f_Oring(npp, Dstar, t, 2))


#doesn't really seem to differ from sum(I_t), but the SSE is different!
plot(O, type = 'o', col='blue', main = 'Oring')
points(Ostar, type = 'o', col='red')




# 4. ripley's k (doesn't currently work when no infecteds. need to fix, like I did with O-ring)

#arbitrarily choose r=2.5
rK <- sapply(1:20, function(t) fK(npp, D, t, R = 2.5, C = 'isotropic'))
rKstar <- sapply(1:20, function(t) fK(npp, Dstar, t, R = 2.5, C = 'isotropic'))
plot(rK, type = 'o', col='blue', main = 'Oring', ylim=c(min(c(rKstar, rK)), max(c(rKstar, rK))))
points(rKstar, type = 'o', col='red')



#compare all 4 summary statistics visually
par(mfrow=c(2,2))
plot(I, type = 'o', col='blue', main = '# infected' )
points(Istar, type = 'o', col='red')

plot(NN, type = 'o', col='blue', main = 'nearest infected neighbor', sub = 'doesnt work, matrix')
points(NNstar, type = 'o', col='red')

plot(O, type = 'o', col='blue', main = 'Oring', ylim=c(min(c(Ostar, O)), max(c(Ostar, O))))
points(Ostar, type = 'o', col='red')

plot(rK, type = 'o', col='blue', main = 'K function', ylim=c(min(c(rKstar, rK)), max(c(rKstar, rK))))
points(rKstar, type = 'o', col='red')



#plot difference function over time.
dI <- data.frame(d = sqrt((I - Istar)^2), t=1:length(I), metric = 'I')
dNN <-data.frame(d = sqrt(colSums((NN - NNstar)^2)), t=1:ncol(NN), metric = 'NN') 
dO <-data.frame(d = sqrt((O - Ostar)^2), t=1:length(O), metric = 'O-ring') 
drK <-data.frame(d = sqrt((rK - rKstar)^2), t=1:length(rK), metric = 'ripleys K') 
df <- bind_rows(dI, dNN, dO, drK)
head(df)
ggplot(df, aes(t, d, group = metric)) +
  geom_point() +
  geom_line(aes(color = metric)) +
  facet_wrap(~metric, scales = 'free_y')



#more summary stats----

# 5. matching states at each location?
#(sum_j ((x_j - y_j))) / length(y) #%matching?
sum(as.numeric(D == Dstar)) / length(D)
colSums(matrix(as.numeric(D == Dstar), nrow = dim(D)[1], ncol = dim(D)[2])) %>% plot(type='l')


#sum of statistics for S and C?

D <- f_sim1(.75, 1, .5) #"True" data
Dstar <- f_sim1(.75, .5, .5) #simulated data


#quick comparison
#compare "data" with simulation
tmp <- data.frame(I = colSums(D=="I"),
                  Istar = colSums(Dstar=="I"),
                  S = colSums(D=="S"),
                  Sstar = colSums(Dstar=="S"),
                  C = colSums(D=="C"),
                  Cstar = colSums(Dstar=="C"),
                  t = 1:tfinal)

par(mfrow = c(1,3))
plot(tmp$t, tmp$I, type = 'o', col = 'blue', main = 'I')
lines(tmp$t, tmp$Istar, type = 'o', col = 'red')
plot(tmp$t, tmp$S, type = 'o', col = 'blue', main = 'S->I')
lines(tmp$t, tmp$Sstar, type = 'o', col = 'red')
plot(tmp$t, tmp$C, type = 'o', col = 'blue', pch=16, main = 'C->I')
lines(tmp$t, tmp$Cstar, type = 'o', col = 'red', pch=16)
par(mfrow=c(1,1))

# 1. n.C + n.S
sqrt(sum((colSums(D=='C') - colSums(Dstar=='C'))^2)) + sqrt(sum((colSums(D=='S') - colSums(Dstar=='S'))^2))

par(mfrow=c(1,2))
plot((colSums(D=='C') - colSums(Dstar=='C'))^2, type='l')
plot((colSums(D=='S') - colSums(Dstar=='S'))^2, type='l')

#make C and S on a more similar scale. 
C <- colSums(D=='C')
Cstar <- colSums(Dstar=='C')
xC <- ((C - Cstar)^2)
S <- colSums(D=='S')
Sstar <- colSums(Dstar=='S')
xS <- (abs(S - Sstar))

plot(xC, type='l')
plot(xS, type='l')

sqrt(sum(xC)) +  sqrt(sum(xS))


# 2. o-ring for C and S?
palv <- viridis::viridis(20)

f_Oring(npp, D, 5) %>% plot
f_Oring(npp, D, 5, 'C') %>% plot
f_Oring(npp, D, 2, 'C') %>% plot(type='l')
f_Oring(npp, D, 4, 'C') %>% points(type='l')

par(mfrow=c(1,3))
test <- lapply(1:20, function(t) f_Oring(npp, D, t, 'C')) 
test[[2]]%>% plot(type='l', ylim = c(-.01, 1), main = 'C')
for(i in 2:19) test[[i]]%>% lines(col = palv[i])

test <- lapply(1:20, function(t) f_Oring(npp, D, t, 'S')) 
test[[2]]%>% plot(type='l', ylim = c(-.01, 1), main = 'S')
for(i in 2:19) test[[i]]%>% lines(col = palv[i])

test <- lapply(1:20, function(t) f_Oring(npp, D, t)) 
test[[2]] %>% plot(type='l', ylim = c(-.01, 1), main = 'I')
for(i in 2:19) test[[i]]%>% lines(col = palv[i])
par(mfrow=c(1,1))

#does radius matter when you're finding the difference?
test <- lapply(1:20, function(t) f_Oring(npp, D, t, 'C')) 
teststar <- lapply(1:20, function(t) f_Oring(npp, Dstar, t, 'C')) 

#difference across peaks. use time 2 to determine pks
pks <- peaks(test[[2]])
plot(test[[2]], type='l')
points(pks, test[[2]][pks], col='red')


#difference across time at 5 peaks
test2 <- sapply(1:20, function(x) test[[x]][pks])
test2[is.na(test2)] <- 0
test2star <- sapply(1:20, function(x) teststar[[x]][pks])
test2star[is.na(test2star)] <- 0

test2diff <- abs(test2 - test2star)

test2diff[1,] %>% plot(type='l', xlab='time', ylab='oring', col=palz[1])
test2diff[2,] %>% lines(type='l', col=palz[2])
test2diff[3,] %>% lines(type='l', col=palz[3])
test2diff[4,] %>% lines(type='l', col=palz[4])
test2diff[5,] %>% lines(type='l', col=palz[5])
dists <- round(sort(unique(round(c(d), 8)))[2:6], 2)
palz <- RColorBrewer::brewer.pal(length(dists), 'BrBG')
legend('topright', legend = dists, fill = palz)


#how would you use o-ring for C and S?
peak1 <- sort(unique(round(c(d), 8)))[2] #first peak
C.O <- sapply(1:20, function(t) f_Oring(npp, D, t, 'C', R = peak1)) 
Cstar.O <- sapply(1:20, function(t) f_Oring(npp, Dstar, t, 'C', R = peak1)) 
SSE(C.O, Cstar.O)
S.O <- sapply(1:20, function(t) f_Oring(npp, D, t, 'S', R = peak1)) 
Sstar.O <- sapply(1:20, function(t) f_Oring(npp, Dstar, t, 'S', R = peak1)) 
SSE(S.O, Sstar.O)

SSE(C.O, Cstar.O) + SSE(S.O, Sstar.O) #total difference?


     