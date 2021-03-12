#organizing which trays are done on each day

#not going to use the same organization as last time since it wasn't exactly stratified. Best way is to split them up numerically into thirds. 

rm(list = ls())
library(tidyverse)

#load data----------------------------------------
mapdf <- read_csv('post_quarantine_trials/output/big_experiment_design_02082021.csv')

#can ignore warning. column "rep" is mostly numbers, but the control trays are called 'control'. I did that so it would print with a good description. 
problems(mapdf) %>% 
  distinct(col, expected, actual, file)



#check out data-----------------------------------
head(mapdf); dim(mapdf)
mapdf$trayID %>% unique %>% sort #should be 178 trays, non-consecutive, in 1s (regular), 200s (extras), 300s (controls)

#count how many seeds of each species
nrow(mapdf) #OMG, 39592 seeds in total!  
mapdf %>% 
  count(spID) %>% arrange(-n)
#1 radish    11914
#2 arugula    8119
#3 basil      5369
#4 red_rom    5172
#5 green_rom  4718
#6 butter     4300



#assign which days you'll plant each tray---------
ntrays <- length(unique(mapdf$trayID))
master <- mapdf %>% 
  distinct(trayID, rand, dens, richness, rep, SD) %>% 
  arrange(trayID) %>% 
  mutate(day = rep(1:3, times = ceiling(ntrays/3))[1:ntrays])

#do controls all on the last day
master$day[master$trayID>300] <- 3

master %>% 
  count(day)

#create datasheets for monitoring-----------------
library(lubridate)
dates <- c('02/08/21', '02/09/21', '02/10/21')
dates <- mdy(dates)
df_trays <- data.frame(tray = master$trayID, 
           day = master$day,
           date_sowed = dates[master$day]) 

#add days and notes columns
df_monitor <- df_trays
daynames <- paste0('day', 0:6*3)
df_monitor[,daynames] <- NA
df_monitor$notes <- NA
df_monitor <- df_monitor %>% arrange(day)
head(df_monitor)


#create datasheet for watering--------------------

df_trays

#add col for each day where weight and sec will be recorded
df_water <- df_trays
daynames <- paste0('day', 3:8*2 + 1)
df_water[,daynames] <- NA
df_water <- df_water %>% arrange(day)
head(df_water)


#Export-------------------------------------------
write_csv(df_monitor, 'post_quarantine_trials/output/monitor_trays_Feb2021.csv', na = '')
write_csv(df_water, 'post_quarantine_trials/output/water_trays_Feb2021.csv', na = '')

