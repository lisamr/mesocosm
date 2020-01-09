#trying to find a better way to enter data. it would be nice to enter data in the form of a matrix, which conserves the spatial arrangements of the plant cells, rather than a list which requires recalculating which cell you're keeping track of. 
#create 3 matrices per tray: 1. spatial arrangement of plants with their ID, 2. germination times, 3. state changes from S to I.  

library(tidyverse)
library(readxl)
rm(list=ls())

#read in data----
pos <- read_xlsx('GH_data/dataentry_matrix.xlsx', sheet=1, na = "NA")
germ <- read_xlsx('GH_data/dataentry_matrix.xlsx', sheet=2, na = "NA")
inf <- read_xlsx('GH_data/dataentry_matrix.xlsx', sheet=3, na = "NA")

#goal: turn the data in matrix form into dataframe form

#seperate each tray into a list----
sep_list <- function(dat){
  rowstarts <- grep("tray", dat$col1)
  ntrays <- length(rowstarts)-1
  dlist <- lapply(1:ntrays, function(i) {
    dat <- dat[(rowstarts[i]+1):(rowstarts[i+1]-1), ]
    dat[,1] <- as.numeric(dat$col1)
    dat
  })
  names(dlist) <- dat$col1[rowstarts][1:ntrays]
  return(dlist)
}
#run list function for each data sheet
posL <- sep_list(pos)
germL <- sep_list(germ)
infL <- sep_list(inf)

#turn matrix into df----
dfL <- lapply(1:length(posL), function(x) {
  d <- data.frame(
    tray = x,
    pos = as.vector(as.matrix(posL[[x]])),
    germ = as.vector(as.matrix(germL[[x]])),
    inf = as.vector(as.matrix(infL[[x]]))
  )
  filter(d, !is.na(pos))
  
})
head(dfL)

#turn list back into dataframe
df <- bind_rows(dfL)

#visualize changes over time----
#germination
df_germ <- df %>% 
  group_by(tray) %>% 
  count(germ) %>% 
  mutate(cumsum = cumsum(n),
         prop = cumsum/cumsum[which(cumsum==max(cumsum))]) 
#plot
df_germ %>% 
  filter(!is.na(germ)) %>% 
ggplot(., aes(germ, prop, color = as.factor(tray))) +
  geom_point() +
  geom_line() 

#infections
df_inf <- df %>% 
  group_by(tray) %>% 
  count(inf) %>% 
  mutate(cumsum = cumsum(n),
         prop = cumsum/cumsum[which(cumsum==max(cumsum))])
#plot
df_inf %>% 
  filter(!is.na(inf)) %>% 
ggplot(., aes(inf, prop, color=as.factor(tray))) +
  geom_point() +
  geom_line()


