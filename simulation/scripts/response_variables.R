#libraries and display
library(raster)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
#display.brewer.all()
pal <-rev(brewer.pal(6, "YlOrRd"))  

#simulation evaluation
#I need reports on response variables for each simulation set

#Response variables
#max dI/dt; time to max dI/dt; %infected; all variables reported for each species and all together. 
#G matrix (communitity R0) = dominant eiganvalue of the N x N matrix. a=Bij*Pij/di. components are transmission rates and average duration of infection period. B=transmission rates, P=#infecting species, d=death rate, which I assume is 1 for this system. 

#load data
comp2 <- read.csv("simulation/outputs/competencyvalues.csv")
abund <- read.csv("simulation/outputs/hostabund.csv")

#simulate
set.seed(2)
h1 <- format.comp(plants = abund$R2)
HP <- HPraster(h1, n = 10)
vv <- simulate(h1, HP, steps=20, f)
animate(HP, vv, pause=0.2, col=pal)

########################################################################
#percent infected
########################################################################
PI <- function(data, HP){
#data=output from simulate; HP=output from HPraster
  
#get data
dat <- data[[1]] #state of each cell by time: N x time
spp <- as.vector(values(HP$host)) #vector of species present in correct

#make tables of infecteds
I <- (dat * spp)
tmp1 <- apply(I, 2, table) %>% t() %>% as.data.frame()
tmp1 <- tmp1[-1] #include only infecteds
tmp1$tot <- apply(tmp1, 1, sum)

#lookup table for species abundances
tmp2 <- table(spp) %>%  as.data.frame()
tmp4 <- cbind.data.frame("tot", sum(tmp2$Freq))
names(tmp4) <- names(tmp2)
tmp2 <- rbind(tmp2, tmp4)

#generate percent infecteds
tmp3 <- matrix(NA, ncol=ncol(tmp1), nrow = nrow(tmp1))
for(i in 1:ncol(tmp1)){
  sp <- names(tmp1)[i]
  sp
  freq <- tmp2[match(sp, tmp2$spp),2]
  freq
  tmp3[,i] <- tmp1[,i]/freq
}

#join with infected table
tmp3 <- as.data.frame(tmp3)
names(tmp3) <-  paste("PI", names(tmp1), sep = ".")
finaltab <- cbind.data.frame(time=as.numeric(row.names(tmp1)), tmp1, tmp3)

return(finaltab)
}

h1 <- format.comp(plants = abund$R4)
HP <- HPraster(h1, n = 10)
vv <- simulate(h1, HP, steps=20, f)
animate(HP, vv,pause = .1, col = pal)
pi1 <- PI(data = vv, HP)
head(pi1)
#melt PI data plot with ggplot

