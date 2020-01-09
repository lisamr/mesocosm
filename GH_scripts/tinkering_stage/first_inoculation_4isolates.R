#monitoring emergence and susceptibility to the four Rhizoc isolates
#plants inoculated after 9 days since sowing, inoculum in PDA for 3 days.

library(dplyr)
library(readxl)
library(ggplot2)
library(tidyr)
library(plotly)

rm(list=ls())
setwd('/Users/lisarosenthal/Box/mesocosm expt/mesocosm.git/')

#set ggplot theme
theme_set(theme_bw()+theme(panel.grid.minor = element_blank()))

#supply columns you want
species <- c("romaine", "endive", "sunflower", "PC_F1", 'PC_micro', "radish", "arugula", "celosia", "beets_micro", "beets_F1", "cilantro", "scallion", "basil", "corn", "popcorn", "wheat", "pea", "clover", "mustard")
family <- c('asteraceae', 'asteraceae', 'asteraceae', 'brassica', 'brassica', 'brassica', 'brassica', 'brassica', 'amaranth', 'amaranth', 'umbellifer', 'allium', 'lamiaceae', 'poaceae', 'poaceae', 'poaceae', 'fabaceae', 'fabaceae', 'brassica')
obs_ease <- c(3,3,2.5, 3, 3, 1, 2, 3, 2, 2, 3, 2, 3, 1, 1, 2, 2, 3, 2)
isolate <- c("Rz1_143N", 'Rz2_RS293', 'Rz3_B1429', 'Rz4_Tat1')

#make dataframe
df <- expand.grid(species=species, isolate=isolate) %>% 
  arrange(species) %>% 
  mutate(family=rep(family, each=4),
         obs_ease=rep(obs_ease, each=4)) %>% 
  select(family, everything()) %>%  #put family at beginning
  #columns to be entered in excel...
  mutate(nplanted=NA,
         emerged_day6=NA,
         emerged_day9=NA,
         infected_day3=NA,
         infected_day4=NA,
         infected_day5=NA,
         infected_day6=NA,
         notes=NA) 
head(df)  
#write.csv(df, 'GH_output/tinkering_stage/firstinoculations_4isolates.csv', row.names = F)


##############################
#After entering in data into excel, need to know how species-isolates differ in susceptibility
df1 <- read_xlsx('GH_data/firstinoculations_4isolates.xlsx')

#first examine emergence. I really don't want any species that suck at germinating.
sum_emergence <- df1 %>% select(family, species, isolate, nplanted, emerged_day6, emerged_day9) %>% 
  group_by(family, species) %>% 
  summarise(total=sum(nplanted),
            day6=sum(emerged_day6),
            day9=sum(emerged_day9)) 

#convert wide to tall for plotting
#notes for gather/spread
#key is the column names being converted into a column of variables
#value is the name of the values
#subtract the columns you don't want melted
sum_emergence <- sum_emergence %>% 
  gather('day', 'n_emerged', -c(family, species, total)) %>% 
  mutate(perc_emerged = n_emerged/total, 
         day=as.integer(gsub('day', '', day)))

#plot emergence
pemerge <- ggplot(sum_emergence, aes(day, perc_emerged, color=species, group=species)) +
  geom_point()+
  geom_line()+
  geom_hline(yintercept=.9, color='firebrick3', lty=2)
plotly::ggplotly(pemerge)

#hard to tell from the plot which ones have emergence rates >.9
#clover also should be in there. had a mishap in the gh.
sum_emergence %>% filter(day==9, perc_emerged>=.9)
sum_emergence %>% filter(day==9) %>% arrange(-perc_emerged)
good_germ <- sum_emergence %>% filter(day==9, perc_emerged>=.9) %>% pull(species)
df1$goodgerm <- ifelse(df1$species %in% good_germ, 1, 0)

#secondly, investigate infection rates
sum_inf <- df1 %>% select(family, species, isolate, emerged_day9, infected_day3, infected_day4, infected_day5, infected_day6) %>% 
  gather('day', 'n_inf', -c(family, species, isolate, emerged_day9)) %>% 
  rename(total=emerged_day9) %>% 
  mutate(day=as.integer(gsub('infected_day', '', day)),
         n_inf=as.numeric(n_inf),
         p_inf=n_inf/total)
head(sum_inf)

#visualize infections 
sum_inf %>% filter(!is.na(p_inf)) %>% 
  ggplot(., aes(day, p_inf, group=isolate)) +
  geom_point() +
  geom_line(aes(color=isolate)) +
  facet_wrap(~species) +
  labs(title='first inoc 4 isolates',
       subtitle = 'plants inoc day9')

sum_inf %>% filter(!is.na(p_inf)) %>% 
  ggplot(., aes(day, p_inf, group=species, color=species)) +
  geom_point() +
  geom_line(aes()) +
  facet_wrap(~isolate) +
  scale_color_viridis_d()

#viz ranked susceptibility curves
test <- sample(1:100, 10)
rank(-test)
rank_suscept <- sum_inf %>% 
  filter(!is.na(p_inf)) %>% 
  filter(day==max(day)) %>% 
  group_by(isolate) %>% 
  mutate(rank = rank(-p_inf))
p1 <- ggplot(rank_suscept, aes(rank, p_inf, fill=species)) +
  geom_col(position=position_dodge2(), width=1) +
  facet_wrap(~isolate) +
  scale_fill_viridis_d()
ggplotly(p1)
#ggsave('GH_plots/rank_suscept.png', p1, device = 'png')
