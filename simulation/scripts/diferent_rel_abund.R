#need to see how relative community competency, richness, and percent infected vary. I don't have plots that are both high richness and high rel. C.C., so I'll need to add some plots witht that combo in order to show that CC controls infections, not something inherent about richness.  

#alter the relative abundances of the design. the function `prepdata2` is the one that assigns relative abundances.

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
i <- 5
des=designs[[i]]
data <- prepdata2(dens = 'sub', rand = F, L=unique(des$L), des)
#run simulation. works :)
S2 <- s2(Rich=4, data, B, des, t=30)
animate(S2$HPraster, S2$simulate)

#check out the innerworkings of prepdata2
prepdata2#getabund4 must define the relative abundances
getabund4#might want to manually adjust the rel abundances rather than relying on prepdata2

getabund4(des$r[des$d=='sub'], 50, 'sub')
