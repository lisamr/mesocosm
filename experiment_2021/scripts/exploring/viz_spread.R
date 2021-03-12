#exploratory first look at the data


#LOAD SOURCE CODE---------------------------------

rm(list=ls())
source('IBM/scripts/IBM_functions.R')
library(cowplot)
theme_set(theme_classic() + theme(plot.title = element_text(hjust = .5))) 

#IMPORT DATA AND SPATIAL DATAFRAME----------------
spdf_list <- readRDS('experiment_2021/output/big_experiment_spdf_list_02082021.RDS')
state_mat_list <- readRDS('experiment_2021/output/state_mat_list.RDS')
treatments_list <- readRDS('experiment_2021/output/treatments_list.RDS')
plant_level_data <- read_csv('experiment_2021/output/plant_level_data.csv')
tray_level_data <- read_csv('experiment_2021/output/tray_level_data.csv')



#plot a tray--------------------------------------
which_tray <- function(x) which(tray_level_data$trayID == x)
tr <- 41 
tfinal <- 18
spp <- unique(plant_level_data$spID)
plotS_I(state_mat_list[[which_tray(tr)]])
plot_spread_map(spdf_list[[which_tray(tr)]], state_mat_list[[which_tray(tr)]], animate = F)
#anim <- plot_spread_map(spdf_list[[64]], state_mat_list[[64]], animate = T, filename = 'experiment_2021/figures/gif_tray67.gif')


#diversity-disease for 4 treatments---------------
print(plant_level_data, width = Inf)
print(tray_level_data, width = Inf)


names(treatments_list)
#[1] "add_det"        "sub_det"        "add_stoch"     
#[4] "sub_stoch"      "single_species"

#function to plot incidence for the 4 trts
plot_4trts <- function(trt, I, n, loess = NULL){
  p <- tray_level_data %>% 
    mutate(incidence = !!I / !!n) %>% 
    filter(trayID %in% treatments_list[[{{trt}}]]$trayID) %>% 
    ggplot(., aes(x = richness, y = incidence)) +
    geom_jitter(height = 0, width = .1, alpha = .5, size = 2, color = 'slateblue') +
    labs(title = names(treatments_list)[trt],
         y = '%infected') + 
    scale_y_continuous(limits = c(0, 1))
  if(!is.null(loess)) p <- p +  geom_smooth(color = grey(.5), lwd = .5)
  return(p)
}

#Community-wide disease incidence
plot_all <- lapply(1:4, function(i) plot_4trts(i, quo(I_all), quo(n_allC), T))
plot_grid(plotlist = plot_all)

#disease incidence of radish
plot_radish <- lapply(1:4, function(i) plot_4trts(i, quo(I_radish), quo(n_radish), T))
plot_grid(plotlist = plot_radish)

#disease incidence of arugula
plot_arugula <- lapply(1:4, function(i) plot_4trts(i, quo(I_arugula), quo(n_arugula), T))
plot_grid(plotlist = plot_arugula)



#Plot infections vs. radish density---------------

plot_vs_radishdens <- function(I, n, Title, df = tray_level_data){
  df %>% 
    filter(trayID < 300) %>% 
    ggplot(., 
           aes(n_radish, !!I/!!n, 
               Infected = !!I, Suscept = !!n - !!I)) +
    geom_point(alpha = .5) +
    geom_smooth(
      method="glm",
      method.args=list(family="binomial"),
      formula = cbind(Infected, Suscept) ~ x
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(y = '% infected', title = Title)
}

p2 <- list(
  plot_vs_radishdens(quo(I_radish), quo(n_radish), 'radish'),
  plot_vs_radishdens(quo(I_arugula), quo(n_arugula), 'arugula'),
  plot_vs_radishdens(quo(I_basil), quo(n_basil), 'basil'),
  plot_vs_radishdens(quo(I_red_rom), quo(n_red_rom), 'red rom'),
  plot_vs_radishdens(quo(I_green_rom), quo(n_green_rom), 'green rom'),
  plot_vs_radishdens(quo(I_butter), quo(n_butter), 'butter')
)

plot_grid(plotlist = p2)


#plot invidence vs overall densiy----------------
plot_vs_overalldens <- function(I, n, Title, df = tray_level_data){
  df %>% 
    filter(trayID < 300) %>% 
    ggplot(., 
           aes(n_allC, !!I/!!n, 
               Infected = !!I, Suscept = !!n - !!I)) +
    geom_point(alpha = .5) +
    geom_smooth(
      method="glm",
      method.args=list(family="binomial"),
      formula = cbind(Infected, Suscept) ~ x
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(y = '% infected', title = Title)
}


p3 <- list(
  plot_vs_overalldens(quo(I_radish), quo(n_radish), 'radish'),
  plot_vs_overalldens(quo(I_arugula), quo(n_arugula), 'arugula'),
  plot_vs_overalldens(quo(I_basil), quo(n_basil), 'basil'),
  plot_vs_overalldens(quo(I_red_rom), quo(n_red_rom), 'red rom'),
  plot_vs_overalldens(quo(I_green_rom), quo(n_green_rom), 'green rom'),
  plot_vs_overalldens(quo(I_butter), quo(n_butter), 'butter')
)

plot_grid(plotlist = p3)


#Viz single species trays-------------------------

singlesp_trays <- plant_level_data %>% 
  filter(trayID %in% treatments_list[[5]]$trayID) %>% 
  distinct(trayID, trayID_consec, spID) %>% 
  mutate(spID = factor(spID, levels = c('radish', 'arugula', 'basil', 'red_rom', 'green_rom', 'butter')))

I_N <- function(tray_consec){
  I <- colSums(state_mat_list[[tray_consec]] == 'I', na.rm = T) 
  N <- colSums(!is.na(state_mat_list[[tray_consec]])) 
  tmp <- cbind(trayID_consec = tray_consec, time = 1:length(I), I, N)
  return(tmp[c(3,6,9,12,15,18),])
}

tmp <- lapply(singlesp_trays$trayID_consec, I_N)
tmp <- as.data.frame(do.call(rbind, tmp))
singlesp_trays2 <- singlesp_trays %>% 
  left_join(tmp)
ggplot(singlesp_trays2, aes(time, I/N, group = trayID, color = spID)) +
  geom_point(alpha = .9) +
  geom_line(alpha = .9) +
  facet_wrap(~spID) +
  scale_color_brewer(palette = 'RdYlBu')



#Visualize disease progress curves of all trays-----

#get attribute and infection data from each tray
tmp <- df2 %>% group_by(trayID, spID, rep) %>% 
  summarise(N = sum(!is.na(state0))) %>% ungroup()

tmp2 <- state_matdf %>% 
  group_by(trayID) %>% 
  summarise_at(vars(V1:V18), .funs = function(x) sum(x %in% 'I')) %>% 
  pivot_longer(V1:V18, names_prefix = 'V', names_to = 'Time', values_to = 'I') %>% 
  mutate(Time = as.integer(Time)) %>% 
  filter(Time %in% seq(3, tfinal, by = 3))

#bind together
df3 <- left_join(tmp, tmp2) %>% 
  mutate(spID = as.factor(spID),
         spID = fct_relevel(spID, spp), 
         propI = I/N)

#viz
ggplot(df3, aes(Time, I, group = trayID, color = spID)) +
  geom_point(alpha = .8, size = 1) +
  geom_line() + 
  facet_wrap(~spID) +
  scale_color_brewer(palette =  'RdYlBu')













