#analyzing inoculation experiment to understand how distance, soil, isolate*species interactions affects secondary transmission
#plants were inoculated 7 days after sowing. inoculum has 2 days to be colonized in PDA plates. 

library(tidyverse)
library(readxl)

#create ID spreadsheet (only need to do once)----
#create spreadsheet to record data
#columns:isolate	species	soil	distance	tray	replicate	state	day	date	notes
library(lubridate)
isolate <- c('Rz3', 'Rz4')
spp <- c('Aru', 'Cel')
soil <- c('coir_lite', 'UC_mod', 'sand')
distance <- c(1,2,4,5)
rep <- 1:10
day <- 0:14 #time since first plant inoculated
date <- seq.Date(ymd(20191213), by='day', length.out = length(day))
names(date) <- day#make a named vector to recode day with date

#make table detailing which trays each treatment is in
expand_grid(isolate, spp, soil, distance, rep) %>% 
  mutate(tray=NA,
         position=NA) %>% 
  write.csv(., 'GH_output/tinkering_stage/dist_secondary_trans_position.csv', row.names = F)

rm(list=ls())
theme_set(theme_bw())

#read in data----
#read in that table (had to look at photos to fig out positions)
ID <- readxl::read_xlsx('GH_data/dist_secondary_trans_position.xlsx')
print(ID, width = Inf)

#merge the position ID table with the data on germination and infection status
pos <- read_xlsx('GH_data/dist_secondary_trans_matrix.xlsx', sheet=1, na = "NA")
germ_d <- read_xlsx('GH_data/dist_secondary_trans_matrix.xlsx', sheet=2, na = "NA")
germ_r <- read_xlsx('GH_data/dist_secondary_trans_matrix.xlsx', sheet=3, na = "NA")
inf_d <- read_xlsx('GH_data/dist_secondary_trans_matrix.xlsx', sheet=4, na = "NA")
inf_r <- read_xlsx('GH_data/dist_secondary_trans_matrix.xlsx', sheet=5, na = "NA")

#seperate each tray into a list----
sep_list <- function(dat){
  #add a line at end of df to designate the end
  dat <- add_row(dat)
  dat$col1[nrow(dat)] <- "last_tray"
  #now define how to seperate the trays and put into list
  rowstarts <- grep("tray", dat$col1)
  ntrays <- length(rowstarts)-1
  dlist <- lapply(1:ntrays, function(i) {
    dat <- dat[(rowstarts[i]+1):(rowstarts[i+1]-1), ]
    dat[,1] <- as.integer(dat$col1)
    dat
  })
  names(dlist) <- dat$col1[rowstarts][1:ntrays]
  return(dlist)
}
#run list function for each data sheet
posL <- sep_list(pos)
germ_dL <- sep_list(germ_d)
germ_rL <- sep_list(germ_r)
inf_dL <- sep_list(inf_d)
inf_rL <- sep_list(inf_r)

#turn matrix into df----
dfL <- lapply(1:length(posL), function(x) {
  d <- data.frame(
    tray = x,
    position = as.vector(as.matrix(posL[[x]])),
    germ_d = as.vector(as.matrix(germ_dL[[x]])),
    germ_r = as.vector(as.matrix(germ_rL[[x]])),
    inf_d = as.vector(as.matrix(inf_dL[[x]])),
    inf_r = as.vector(as.matrix(inf_rL[[x]]))
  )
  filter(d, !is.na(position))
  
})
head(dfL)

#turn list back into dataframe
df <- bind_rows(dfL)

#merge experiment data with ID attributes----
df2 <- right_join(ID, df, by=c("tray", 'position')) 

#visualize changes over time----
#germination
df_germ <- df2 %>% 
  group_by(spp) %>% 
  count(germ_d) %>% 
  mutate(cumsum = cumsum(n), prop = cumsum/cumsum[which(cumsum==max(cumsum))]) 
#plot
df_germ %>% 
  filter(!is.na(germ_d)) %>% 
  ggplot(., aes(germ_d, prop, fill = as.factor(spp))) +
  geom_col(position = 'dodge2') 

#infections
df_inf <- df2 %>% 
  group_by(spp, soil) %>% 
  count(inf_d) %>% 
  mutate(cumsum = cumsum(n), prop = cumsum/cumsum[which(cumsum==max(cumsum))]) 
#plot
df_inf %>% 
  filter(!is.na(inf_d)) %>% 
  ggplot(., aes(inf_d, prop, color = soil)) +
  geom_point() +
  geom_line() +
  facet_grid(~spp)

#analyze the data----

