
comp <- read_table('post_quarantine_trials/data/prop_inf_compNov2020.txt')
comp <- data.frame(comp = comp$mean,
                   spID = c('radish', 'arugula', 'basil', 'green_rom', 'red_rom', 'butter'))
trays <- read_csv('GH_data/real_experiment/water_weight_Jan2020.csv')
df <- readxl::read_xlsx('GH_data/real_experiment/big_experiment_03302020_recorded.xlsx')


#clean up individual level data
df$spID %>% unique
df2 <- df %>% 
  mutate(spID = recode(spID, 
                       'pac_choy' = 'basil', 
                       'romaine' = 'green_rom',
                       'basil' = 'red_rom',
                       'clover' = 'butter')) %>% 
  left_join(comp, by = "spID")


#generate tray level data
design <- df2 %>% 
  distinct(trayID, rand, dens, richness, rep, SD)
trayinfo <- df2 %>% 
  #filter(!infected %in% c("NA", "X")) %>% 
  group_by(trayID) %>% 
  summarise(CC = sum(comp), nplants = n()) %>% 
  left_join(trays) %>% 
  left_join(design)
df2 <- df2 %>% left_join(trayinfo)
#ggplot(trayinfo, aes(CC, weight)) +geom_point() +geom_smooth(method = 'lm') #+facet_wrap(~day_planted)


#calculate perc inf for each species in each tray
percI_sp <- df2 %>% 
  filter(!infected %in% c("NA", "X")) %>% 
  group_by(trayID, richness, spID) %>% 
  summarise(I = sum(infected == "I"), 
            S = sum(infected == "S"),
            percI = I / (S + I)) %>% 
  left_join(trayinfo)

percI_sp%>% filter(percI > 0) %>% 
  ggplot(., aes(CC, percI, succ = I, fail = S)) +
  geom_point(alpha = .5) +
  facet_wrap(~spID) +
  geom_smooth(
    method="glm",
    method.args=list(family="binomial"),
    formula = cbind(succ, fail) ~ x
  )


#check out traylevel percent infected
percI_tray <- df2 %>% 
  filter(!infected %in% c("NA", "X")) %>% 
  group_by(trayID, richness) %>% 
  summarise(I = sum(infected == "I"), 
            S = sum(infected == "S"),
            percI = I / (S + I)) %>% 
  ungroup() %>% 
  filter(percI > 0) %>% 
  left_join(trayinfo) 
  
ggplot(percI_tray, aes(CC, percI, color = richness)) +
  geom_point(alpha = .7)+
  scale_color_viridis_c()
ggplot(percI_tray, aes(richness, percI)) +
  geom_jitter(height = 0, width = .2, alpha = .5)+
  facet_grid(~rand)



#calculate FOI[i]
f_FOI <- function(tray){
  tmp <- df2 %>% filter(trayID == tray)
  dmat <- as.matrix(dist(cbind(tmp$x, tmp$y), upper = T, diag = T))
  dmat2 <- exp(-1*(dmat^2))
  diag(dmat2) <- 0
  FOI = colSums(dmat2*tmp$comp) 
  return(FOI)
}
FOI <- sapply(unique(df2$trayID), f_FOI)
df2$FOI <- unlist(FOI)
pal <- RColorBrewer::brewer.pal(4, 'RdYlBu')
ggplot(df2 %>% filter(trayID == 90), aes(x, y)) +
  geom_point(aes(color = spID, size= FOI)) + 
  scale_color_manual(values = pal[c(2,3,4,1)]) +
  coord_equal()

ggplot(df2 %>% filter(trayID == 107), aes(infected=="I", FOI)) +
  geom_boxplot()

df2 %>% 
  filter(nplants == 238) %>% 
  ggplot(., aes(CC, FOI)) +
  geom_point(alpha = .3)
