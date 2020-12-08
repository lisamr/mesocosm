#host competency trial in Nov 2020
#1. PLANT AT FULL CAPACITY. SMALLER TRAYS LEADS TO CONSISTENTLY LOWER COMPETENCY VALUES.
#2. VARIATION AFTER APPROXIMATELY 4 AND 6 SEEMS TO LEVEL OFF for low and high, respectively. 5 for high might be okay too, but just depends on time and seed quantities.

rm(list=ls())
source('IBM/scripts/IBM_functions.R')
library(cowplot)
library(lme4)
library(tidyverse)
theme_set(theme_classic())


#DESIGN PARAMETERS----
#distances wanted
Dist <- c(1.7)

#spp used (high, med, low competency species)
spp <- c('radish', 'mustard', 'pac choy', 'arugula', 'basil', 'green romaine', 'butter lettuce', 'red romaine', 'clover')

#design
#6 trays for high/med competency species, 4 for low
design <- rbind(expand.grid(species = spp[1:5], rep = 1:6),
                expand.grid(species = spp[-(1:5)], rep = 1:4)) %>% 
  arrange(species) %>% 
  mutate(distance = Dist,
         trayID = 1:nrow(.)) 


#MAKE TRAYS----
#interplanting distance is just 1.7cm

get_grid <- function(x){
  design_x <- design[x,]
  which_spp <- which(spp == design_x$species)
  grid <- sample_community(which_spp, .1, Dist)
  grid$trayID <- design_x$trayID
  grid$rep <- design_x$rep
  return(grid)
}

#simulate trays 
set.seed(2020)
#normal dimensions
grid_list <- list(NULL)
for(i in 1:nrow(design)){
  grid_list[[i]] <- get_grid(i)
}

#print maps, 2 per page.
plot_list <- lapply(grid_list, function(x){
  plot_maps(x) + scale_fill_manual(values = 'grey80')
})
#pdf('post_quarantine_trials/figures/maps_competency_Nov2020.pdf')
#plot_list
#dev.off()

#output of dataframe
d <- lapply(grid_list, function(x) x@data) %>% bind_rows() 
write_csv(d, 'post_quarantine_trials/data/competency_Nov2020_empty.csv')

#ID labels start in left corner and move horizontally, ending in upper right corner.
d %>% filter(trayID==1) %>% 
  ggplot(., aes(x, y)) +
  geom_point(size=.1) +
  geom_text(aes(label = ID))



#LOAD RECORDED DATA----
df <- as.data.frame(read_csv('post_quarantine_trials/data/competency_Nov2020_recorded.csv'))
head(df)
spp <- c('radish', 'arugula', 'basil', 'green romaine',  'red romaine', 'butter lettuce')

#day_infected == 99 (still susceptible), == 98 (didnt germinate)
df2 <- df %>% 
  filter(spID %in% spp,
         trayID != 29) %>%
  mutate(day_infected = ifelse(is.na(day_infected), '99', day_infected),
         day_infected = as.numeric(recode(day_infected, 'x' = '98', '*' = '98')),
         state0 = ifelse(day_infected==98, NA, state0)) 
head(df2)

df2 %>% group_by(spID) %>% 
  summarise(ntrays = length(unique(trayID)))


#Visualize with IBM functions-----

#generate a matrix of state changes. 
tfinal <-18

state_mat <- matrix(df2$state0, ncol=tfinal, nrow=nrow(df2))
for(i in 1:nrow(state_mat)){
  if(df2$day_infected[i]<=tfinal){
    state_mat[i,df2$day_infected[i]:tfinal] <- 'I'
  }
}
state_matdf <- bind_cols(trayID = df2$trayID, as.data.frame(state_mat))


#choose a tray
tray <- 21
state_mat_ind <- state_matdf %>% filter(trayID == tray) %>% dplyr::select(-trayID) %>%   as.matrix 
colnames(state_mat_ind) <- NULL

#plot it
plotS_I(state_mat_ind)
plot_spread_map(grid_list[[tray]], state_mat_ind, animate = F, alpha_intensity = .2) 



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





#estimate disease-----
#1. use a binomial model to calculate incidence at tfinal
#2. estimate speed of spread with a logistic model


#BINOMIAL MODEL----

df_tfinal <- df3 %>% 
  filter(Time == max(Time)) %>% 
  mutate(S = N - I)

m1 <- glm(cbind(I, S) ~ -1 + spID, data = df_tfinal, family = binomial())
summary(m1)

m1coefs <- as.data.frame(plogis(cbind(coef(m1), confint(m1)) ) )
names(m1coefs) <- c('mean', 'lower', 'upper')
rownames(m1coefs) <- NULL
m1coefs <- m1coefs %>% 
  mutate(spID = as.factor(spp),
         spID = fct_relevel(spID, spp))

p1 <- ggplot(m1coefs, aes(spID, mean )) +
  geom_pointrange(aes(ymin = lower, ymax = upper, color = (spID))) +
  scale_color_brewer(palette =  'RdYlBu') +
  #scale_color_gradientn(colours = palwes) +
  #scale_x_discrete(limits = rev) +
  scale_y_continuous(limits = c(0, 1), n.breaks = 6) +
  labs(y = 'Proportion infected', x = 'Species') 
p1 

#LOGISTIC FUNCTION----
#dN/dt = mumax*N(1-N/K)

#try to fit a custom defined function with nlme. Will be useful in the bigger experiment when you want to include covariates for richness*treatment. 
library(nlme)

#the formula for the models
F1 <- formula(propI ~ K*y0 / (y0 + (K - y0)*exp(-r*Time)) | spID)
m2.1 <- nlsList(F1, data=df3, start=list(K=1, y0=.001,r=0.5))
summary(m2.1)

#fix K so you can directly compare r
Kfixed = 1
F2 <- formula(propI ~ Kfixed*y0 / (y0 + (Kfixed - y0)*exp(-r*Time)) | spID) 
m2.2 <- nlsList(F2, data=df3, start=list(y0=.001,r=0.2))
summary(m2.2)


#see predictions. they aren't that different so I hope I can get away with fixing K=1.
df3$predicted1 <- predict(m2.1) #solid line
df3$predicted2 <- predict(m2.2) #dashed line
p2 <- ggplot(df3, aes(Time, propI, group = trayID, color = spID)) +
  geom_point() + 
  geom_line(aes(y = predicted1)) +
  geom_line(aes(y = predicted2), lty = 2) +
  facet_wrap(~spID) +
  scale_color_brewer(palette =  'RdYlBu') 


#get summary and plot
x <- confint(m2.2)
m2coefs <- as.data.frame(cbind(coef(m2.2)[,2], t(sapply(x, function(x) x['r',]))))
names(m2coefs) <- c('mean', 'lower', 'upper')
rownames(m2coefs) <- NULL
m2coefs <- m2coefs %>% 
  mutate(spID = as.factor(spp),
         spID = fct_relevel(spID, spp))

p3 <- ggplot(m2coefs, aes(spID, mean )) +
  geom_pointrange(aes(ymin = lower, ymax = upper, color = (spID))) +
  scale_color_brewer(palette =  'RdYlBu') +
  #scale_color_gradientn(colours = palwes) +
  scale_y_continuous(limits = c(0, .5)) +
  labs(y = 'Intrinsic growth rate', x = 'Species') 
p3 






#plot everything together------
library(cowplot)
grid1 <- plot_grid(
  p1 + theme(
    legend.position = 'none',  
    axis.title.x = element_blank(), 
    axis.text.x = element_blank()),
  p3 + theme(
    legend.position = 'none', 
    axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)), 
  nrow = 2, rel_heights = c(.7, 1))
grid1
finalplot <- plot_grid(p2+ theme(legend.position = 'none'),
          grid1, rel_widths = c(1, .6), scale = .9)

finalplot #solid line = all logistic curve pars, dashed = fixed K at 1. 
ggsave('post_quarantine_trials/figures/competency_trials.pdf' , finalplot, width = 6.5, height = 4)




