library(spatstat)

head(agents)
#Convert to a point pattern
npp <- ppp(x=agents$x, y = agents$y, c(min(agents$x), max(agents$x)),
    c(min(agents$y), max(agents$y)))
marks(npp) <- as.factor(D[,4]) #states at an arbitrary time point, used for example. 



#how does #Infecteds change over time? are these spatial summaries any better than counting infecteds?
I <- colSums(D=='I')
plot(1:20, I, type = 'o', pch=16)

#O-ring-----

#O-ring statistic: proportional to the pair correlation function, which tells you the density of points across a range of distances from Infecteds to infecteds. Not sure how method changes the results, but 'spar' is a smoothing parameter. 
K12 <- Kcross(npp, 'I', 'I')
g12 <- pcf(K12,  method="d", spar=0.7)
lambda2 <- summary(npp)$marks['I',"intensity"]
Oring <- eval.fv(lambda2*g12)

par(mfrow=c(1,2))
plot(npp)
plot(Oring)
par(mfrow=c(1,1))


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

f_Oring <- function(ppp, t, Spar=.7){
  marks(npp) <- as.factor(D[,t]) #update time
  K12 <- Kcross(npp, 'I', 'I') #calculate o-ring stat 
  g12 <- pcf(K12,  method="a", spar=Spar)
  lambda2 <- summary(npp)$marks['I',"intensity"]
  Oring <- eval.fv(lambda2*g12)
  df <- data.frame(Oring = Oring$pcf, r = Oring$r, t = t)
  
  return(df)
}

tmpO <- lapply(2:20, function(t) f_Oring(npp, t)) %>% 
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
plot(x)
points(pks, x[pks], col='red', pch = 16 )

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
fK <- function(Marks, correction=c("border", "isotropic", "Ripley", "translate")){
  marks(npp) <- as.factor(Marks)
  Kest(subset(npp, marks == 'I'), correction) 
}

par(mfrow=c(2,2))
fK(D[,5]) %>% plot
fK(D[,10]) %>% plot
fK(D[,15]) %>% plot
fK(D[,20]) %>% plot
par(mfrow=c(1,1))

#compare K over time
fKtime <- function(t, C){
  tmp <- fK(D[,t], C)
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
  
