#analyzing inoculation experiment to understand how distance, soil, isolate*species interactions affects secondary transmission
#plants were inoculated 7 days after sowing. inoculum was colonized for 3 days on PDA plates. 
library(forcats)
library(tidyverse)
library(readxl)

rm(list=ls())
theme_set(theme_bw())

#read in data----
#read in that table (had to look at photos to fig out positions)
ID <- read_xlsx('GH_data/dist_primary_trans.xlsx', sheet = 1)
print(ID, width = Inf)

#merge the position ID table with the data on germination and infection status
pos <- read_xlsx('GH_data/dist_primary_trans.xlsx', sheet=2, na = "NA")
germ <- read_xlsx('GH_data/dist_primary_trans.xlsx', sheet=3, na = "NA")
inf <- read_xlsx('GH_data/dist_primary_trans.xlsx', sheet=4, na = "NA")

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
germL <- sep_list(germ)
infL <- sep_list(inf)


#turn matrix into df----
dfL <- lapply(1:length(posL), function(x) {
  d <- data.frame(
    tray = x,
    position = as.vector(as.matrix(posL[[x]])),
    germ = as.vector(as.matrix(germL[[x]])),
    inf = as.vector(as.matrix(infL[[x]]))
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
df_germ2 <- df2 %>% 
  group_by(species) %>% 
  count(germ) %>% 
  mutate(cumsum = cumsum(n), prop = cumsum/cumsum[which(cumsum==max(cumsum))]) 
df_germ2 %>% 
  filter(!is.na(germ)) %>% 
  ggplot(., aes(germ, prop, color = as.factor(species))) +
  geom_point() +
  geom_line() 

#infections
df_inf <- df2 %>% 
  filter(germ<=7) %>% 
  group_by(species) %>% 
  count(inf) %>% 
  mutate(cumsum = cumsum(n), prop = cumsum/cumsum[which(cumsum==max(cumsum))]) 
#plot
df_inf %>% 
  filter(!is.na(inf)) %>% 
  ggplot(., aes(inf, prop, color=species)) +
  geom_point() +
  geom_line() 

#analyze the data----
#does soil affect infection rate?
library(brms)
library(ggridges)
library(tidybayes)

#prep data. just look at final time infecteds.
df_inf2 <- df2 %>% 
  mutate(
    inoculated_TF = ifelse(germ==4, T, F),
    inf_TF = ifelse(inf>0, T, F)) %>% 
  group_by(isolate, species, soil) %>% 
  summarise(ninf = sum(inf_TF, na.rm = T),
            ninoc = sum(inoculated_TF, na.rm = T))

#check out priors
f1 <- bf(ninf | trials(ninoc) ~ species*soil, family = binomial)
f2 <- bf(ninf | trials(ninoc) ~ species + soil, family = binomial)
get_prior(f1, df_inf2)

m1 <- brm(data = df_inf2, family = binomial, formula = f1,
          prior = c(prior(normal(0, 1), class = b), prior(normal(0, 1), class = Intercept)), iter = 4000, warmup = 1000, chains = 4, cores = 4)
m2 <- brm(data = df_inf2, family = binomial, formula = f2,
          prior = c(prior(normal(0, 1), class = b), prior(normal(0, 1), class = Intercept)), iter = 4000, warmup = 1000, chains = 4, cores = 4)
#check out model
plot(m1)
pp_check(m1) #model fit looks good

#coef plot
m1;m2 #no differences between species or soil

#predictions
newd <- expand.grid(species=c('arugula', 'PC_F1'), ninoc=1, soil=c("coir_lite", "sand", "UC", "UC_mod"))
fitted <- add_fitted_draws(newd, m2)
#plot
ggplot(fitted, aes(x=species, y=.value, color=soil)) +
  stat_pointinterval(.width = .9, position = position_dodge2(width = .3)) +
  scale_y_continuous(limits=c(0,1))
