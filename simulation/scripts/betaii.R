#trying to get g-matrix values from communities. G-matrix elements are betaij * average duration of infection. I don't have betaij, so i'll estimate betaii instead from the simulations. 

#betaii can equal 1) perc.infect at time-final OR 2) avg.exposure/inf.rate

#run simulations with monospecific communities to get Bii
dat <- prepdata("add", F)
data1 <- dat

#get results out for monospecific runs at different distances
#B2: rel magnitude=c(1, .7, .5, .4, .2, .1 )
res=tibble(species=0, n=0, n.I=0, pI=0, max.dI=0, tmax.dI=0, B=0, dist=0)
for(d in 1:length(dat$D)){ #for every distance
  for(s in 1:6){ #for every species
    data1$A$A <- s
    nreps=10
    for(i in 1:nreps){#for every rep
      mono1 <- s2(d, data1, B2, design)
      a <- mono1$resp2[1,]
      a$B <- "B2"
      a$dist <-dat$D[d]
      res <- rbind(res, a)  #get results
    }
  }
}
#B3: rel magnitude=c(1, .3, .2, .1, 0, 0 )
for(d in 1:length(dat$D)){ #for every distance
  for(s in 1:6){ #for every species
    data1$A$A <- s
    nreps=10
    for(i in 1:nreps){#for every rep
      mono1 <- s2(d, data1, B3, design)
      a <- mono1$resp2[1,]
      a$B <- "B3"
      a$dist <- dat$D[d]
      res <- rbind(res, a)  #get results
    }
  }
}

res4 <- res
res4 <- res4 %>% filter(!B==0)
head(res4)

res4 <- res4 %>%
  mutate(new.I=n.I-round(.1*n), pI2=new.I/(n-round(.1*n)))

res4 %>%
  group_by(B, species, dist) %>% 
  summarise(n=length(species), R0i=mean(pI2), se=sd(pI2)/sqrt(n) )

ggplot(filter(res4), aes(species, pI2, group=interaction(dist, species), fill=as.factor(dist)))+
  geom_boxplot()+
  facet_wrap(~B)

#Maybe B_ii is actually average exposure/overall infection rate
dat <- prepdata("add", F)
data1 <- dat
df <- tibble(avg.E=NA, inf.rate=NA, Bii=NA, species=NA, beta=NA, dist=NA)
for(d in 1:length(dat$D)){ #for every distance
  for(s in 1:6){ #for every species
    data1$A$A <- s
    nreps=10
    for(i in 1:nreps){#for every rep
      mono1 <- s2(d, data1, B2, design)
      
      avE <- mono1$resp1 %>% 
        filter(species==s) %>% 
        summarise(avg.E=mean(n.E)) 
      
      df1 <- mono1$resp1 %>% 
        filter(species==s, time %in% c(1, max(time))) %>% 
        summarise(avg.E=avE$avg.E, inf.rate=diff(n.I)/diff(time), Bii=inf.rate/avg.E) %>% 
        as.data.frame() %>% 
        mutate(species=s, beta="B2", dist=dat$D[d])
      
      df <- rbind(df, df1)  #get results
    }
  }
}
for(d in 1:length(dat$D)){ #for every distance
  for(s in 1:6){ #for every species
    data1$A$A <- s
    nreps=10
    for(i in 1:nreps){#for every rep
      mono1 <- s2(d, data1, B3, design)
      
      avE <- mono1$resp1 %>% 
        filter(species==s) %>% 
        summarise(avg.E=mean(n.E)) 
      
      df1 <- mono1$resp1 %>% 
        filter(species==s, time %in% c(1, max(time))) %>% 
        summarise(avg.E=avE$avg.E, inf.rate=diff(n.I)/diff(time), Bii=inf.rate/avg.E) %>% 
        as.data.frame() %>% 
        mutate(species=s, beta="B3", dist=dat$D[d])
      
      df <- rbind(df, df1)  #get results
    }
  }
}
df <- df %>% filter(!is.na(beta))
#rbind/add in distance 2 cm (substittuve design)
dat <- prepdata("sub", F)
data1 <- dat
for(d in 1){ #for just the first distance
  for(s in 1:6){ #for every species
    data1$A$A <- s
    nreps=10
    for(i in 1:nreps){#for every rep
      mono1 <- s2(d, data1, B2, design)
      
      avE <- mono1$resp1 %>% 
        filter(species==s) %>% 
        summarise(avg.E=mean(n.E)) 
      
      df1 <- mono1$resp1 %>% 
        filter(species==s, time %in% c(1, max(time))) %>% 
        summarise(avg.E=avE$avg.E, inf.rate=diff(n.I)/diff(time), Bii=inf.rate/avg.E) %>% 
        as.data.frame() %>% 
        mutate(species=s, beta="B2", dist=dat$D[d])
      
      df <- rbind(df, df1)  #get results
    }
  }
}
for(d in 1){ #for just the first distance
  for(s in 1:6){ #for every species
    data1$A$A <- s
    nreps=10
    for(i in 1:nreps){#for every rep
      mono1 <- s2(d, data1, B3, design)
      
      avE <- mono1$resp1 %>% 
        filter(species==s) %>% 
        summarise(avg.E=mean(n.E)) 
      
      df1 <- mono1$resp1 %>% 
        filter(species==s, time %in% c(1, max(time))) %>% 
        summarise(avg.E=avE$avg.E, inf.rate=diff(n.I)/diff(time), Bii=inf.rate/avg.E) %>% 
        as.data.frame() %>% 
        mutate(species=s, beta="B3", dist=dat$D[d])
      
      df <- rbind(df, df1)  #get results
    }
  }
}


ggplot(filter(df), aes(species, Bii, group=interaction(dist, species), fill=as.factor(dist)))+
  geom_boxplot()+
  facet_wrap(~beta)

R0i <- df %>%
  group_by(beta, species, dist) %>% 
  summarise(n=length(species), R0i=mean(Bii), se=sd(Bii)/sqrt(n) )

write.csv(R0i, "simulation/outputs/R0i.csv", row.names = F)
