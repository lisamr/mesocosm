#create plant level and tray level datasets


#LOAD SOURCE CODE---------------------------------

rm(list=ls())
source('IBM/scripts/IBM_functions.R')


#IMPORT DATA AND SPATIAL DATAFRAME----------------
spdf_list <- readRDS('experiment_2021/output/big_experiment_spdf_list_02082021.RDS')
plantdf <- read_csv('experiment_2021/data/plantdata_Feb2021.csv') 
traydf <- readxl::read_xlsx('experiment_2021/data/traydata_Feb2021.xlsx') %>% 
  mutate(plate = as.integer(plate))
temps <- readxl::read_xlsx('experiment_2021/data/temps.xlsx', sheet = 1)
dome_locations <- readxl::read_xlsx('experiment_2021/data/temps.xlsx', sheet = 2)


#assign treatments to trays-----------------------

trayID_consec <- plantdf %>% 
  distinct(trayID) %>% 
  mutate(trayID_consec = row_number())#consecutive ID number

#trays excluding controls and extras
tray_trt <- plantdf %>% 
  left_join(trayID_consec) %>% 
  group_by(trayID) %>% 
  mutate(nplants = n()) %>% 
  ungroup() %>% 
  distinct(trayID,trayID_consec, rand, dens, richness, rep, SD, nplants) %>%
  filter(SD == .5 | is.na(SD), !is.na(rep))

#1. sub/det, 2. add/det, 3. sub/stoch, 4. add/stoch. all should have (at least) 40 trays per treatment. I omitted some trays because they were repeated across treatment...need to manually reassign.
add_det <- tray_trt %>% 
  filter(rand == 'det', dens == 'add') %>% 
  select(-rand, -dens)
sub_det <- tray_trt %>% 
  filter((rand == 'det'& dens == 'sub') | 
           (rand == 'det' & dens == 'add' & richness == 4))%>%
  select(-rand, -dens)
add_stoch <- tray_trt %>% 
  filter((trayID %in% 1:10 ) | #low dens radish. only need 2
           (rand == 'stoch' & dens == 'add') ) %>% 
  select(-rand, -dens)
sub_stoch <- tray_trt %>% 
  filter((rand == 'stoch' & dens == 'sub') |
           (trayID %in% 21:30 ) )%>% #med dens radish. only need 1
  select(-rand, -dens)

#single species trays with 238 plants
single_species <- tray_trt %>% 
  filter(trayID %in% c(21:40) | 
           trayID > 200 & trayID < 300) %>% 
  select(-rand, -dens)


treatments <- list(add_det = add_det, 
                   sub_det = sub_det, 
                   add_stoch=add_stoch, 
                   sub_stoch=sub_stoch, 
                   single_species=single_species)
names(treatments)






#cleaned up plant data----------------------------

#get last infection state
plantdf2 <- plantdf %>% 
  left_join(trayID_consec) %>% 
  mutate(day_infected = case_when(
    is.na(day_infected) ~ '99', 
    day_infected %in% c('x', "X") ~ "NA",
    T ~ day_infected), 
    day_infected = as.integer(day_infected), 
    state_final = ifelse(day_infected <= 18, 'I', 'S')) %>% 
  select(-NOTES, -rand, -dens, NOTES) %>% 
  rename(notes = NOTES)

#check for mistakes and fix in the excel spreadsheet
plantdf2$day_infected %>% unique #should only by 6, 9, 12, 15, 18, 99, NA
acceptable <- c(6, 9, 12, 15, 18, 99, NA)
plantdf2 %>% 
  filter(!day_infected %in% acceptable) %>% 
  print(width = Inf)


print(plantdf2, width = Inf)





#Create states matrix-----------------------------

#good for plotting and using the IBM functions
#put into the form of the IBM simulation. 
#generate a matrix of state changes. That can be plotted. Also bind state matrix to the original attribute data. 
tfinal <-18
spp <- unique(plantdf2$spID) #names of species

#fill in the state matrix
state_mat <- matrix("S", ncol=tfinal, nrow=nrow(plantdf2))
for(i in 2:(tfinal)){
  state_mat[,i] <- state_mat[,i-1]#current time same as last time
  I <- which(plantdf2$state_final=='I' & plantdf2$day_infected==i)#find infecteds
  state_mat[I,i] <- "I"#update state if infected
}
#replace first column with challenged status and rows as NA
state_mat[,1] <- plantdf2$state0_v2 #some mistakes were made in planting, so use corrected 'v2' 
state_mat[is.na(plantdf2$state_final),] <- NA

#split state matrix by tray
rows <- plantdf2 %>%
  mutate(row=row_number()) %>%  
  group_by(trayID_consec) %>% 
  summarise(firstrow=first(row), lastrow=last(row))

state_mat_list <- list(NULL)
for(i in 1:nrow(rows)){
  whichtray <- rows %>% filter(trayID_consec == rows$trayID_consec[i])
  state_mat_list[[i]] <- state_mat[c(whichtray$firstrow:whichtray$lastrow),]
}

#plot a tray
which_tray <- function(x) which(trayID_consec$trayID == x)
plotS_I(state_mat_list[[which_tray(144)]])
plot_spread_map(spdf_list[[which_tray(144)]], state_mat_list[[which_tray(144)]], animate = F, alpha_intensity = .1)
#anim <- plot_spread_map(spdf_list[[which_tray(144)]], state_mat_list[[which_tray(144)]], animate = T, filename = 'experiment_2021/figures/gif_tray67.gif')



#Collate tray level data--------------------------
temps2 <- temps %>% 
  left_join(dome_locations) %>% 
  mutate(x1 = ifelse(side == 'R', x + .3, x)) %>% 
  select(trayID:side, x, x1, y, temp_PM, temp_AM) 
traydf2 <- traydf %>% 
  rename(trayID = tray) %>% 
  left_join(trayID_consec) %>% 
  left_join(temps2) %>% 
  arrange(trayID) %>% 
  select(trayID, trayID_consec, everything())

#summarized info for each tray
tray_species1 <- plantdf2 %>% 
  group_by(trayID) %>% 
  summarise(n_allP = n(), #planted
            n_allC = sum(!is.na(day_infected)), #counted (germinated or not removed or eaten by bugs)
            I_all = sum(state_final == 'I', na.rm = T)
            )

tray_species2 <- plantdf2 %>% 
  group_by(trayID, spID) %>% 
  summarise(n = sum(!is.na(day_infected)), 
            I = sum(state_final == 'I', na.rm = T)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = spID, values_from = c(n, I)) %>% 
  replace(., is.na(.), 0)

richnesslevels <- plantdf %>% 
  distinct(trayID, richness, SD)

tray_df3 <- left_join(traydf2, richnesslevels) %>% 
  left_join(., left_join(tray_species1, tray_species2)) %>% 
  select(-notes, notes)   #puts notes at the end

#Export data--------------------------------------

saveRDS(state_mat_list, 'experiment_2021/output/state_mat_list.RDS') 
saveRDS(treatments, 'experiment_2021/output/treatments_list.RDS')
write_csv(plantdf2, 'experiment_2021/output/plant_level_data.csv')
write_csv(tray_df3, 'experiment_2021/output/tray_level_data.csv')




