#making sense of 1010 trial planted 2/17/20. see page 75 in lab notebook.

library(tidyverse)

#enter data
tray <- c(2,3,4,5,1,8,9,6,10,7)
richness <- c(1,2,4,6,6,6,6,1,1,1)
n.ind <- c(224, 224, 224, 224, 323, 224, 224, 99, 99, 224)
dist <- c(1.75, 1.75, 1.75, 1.75, 1.4, 1.75, 1.75, 2.7, 2.7, 1.75)
dens <- c('sub', 'sub', 'sub', 'sub', 'add', 'sub', 'sub', 'add', 'add', 'sub')
n.I <- c(82, 60, 54, 35, 77, 77, 45, 19, 15, 41)
species <- list(
  c('radish'),
  c('radish', 'arugula'),
  c('radish', 'arugula', 'PC', 'romaine'),
  c('radish', 'arugula', 'PC', 'romaine', 'clover', 'basil'),
  c('radish', 'arugula', 'PC', 'romaine', 'clover', 'basil'),
  c('radish', 'mustard', 'arugula', 'PC', 'romaine', 'clover'),
  c('mustard', 'arugula', 'PC', 'romaine', 'clover', 'basil'),
  c('radish'),
  c('mustard'),
  c('mustard') )

#put into dataframe
df <- data.frame(tray, richness, dist, dens, n.ind, n.I)
df$species <- species
df <- df %>% mutate(percI = n.I/n.ind) %>% arrange(tray)

#PLOT----
#PART 1: the following should help me figure out best composition order.
#substitutive, radish first. Looks great! Going with this order. 
df %>% filter(tray %in% c(2,3,4,5)) %>% 
ggplot(., aes(richness, n.I)) +
  geom_point() + geom_line()

#substitutive, mustard first. Does the oppositite of what I want. :(
df %>% filter(tray %in% c(7,9)) %>% 
  ggplot(., aes(richness, n.I)) +
  geom_point() + geom_line()
#substitutive, radish then mustard. Not different enough.
df %>% filter(tray %in% c(2,8)) %>% 
  ggplot(., aes(richness, n.I)) +
  geom_point() + geom_line()

#PART 2: the following should help me figure out how steep the density gradient needs to be. For additive assembly.
#radish first: 1sp, 2.7cm -> 4spp., 1.75cm -> 6spp., 1.4cm
df %>% filter(tray %in% c(6, 4, 1)) %>% 
  ggplot(., aes(richness, n.I)) +
  geom_point() + geom_line()
df %>% filter(tray %in% c(6, 4, 1)) %>% 
  ggplot(., aes(richness, percI)) +
  geom_point() + geom_line()

#mustard first: 1sp, 2.7cm -> 6spp., 1.75cm
df %>% filter(tray %in% c(10, 9)) %>% 
  ggplot(., aes(richness, n.I)) +
  geom_point() + geom_line()
df %>% filter(tray %in% c(10, 9)) %>% 
  ggplot(., aes(richness, percI)) +
  geom_point() + geom_line()

