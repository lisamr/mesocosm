library(readxl)
library(tidyverse)
theme_set(theme_classic())

sh1 = read_xlsx('post_quarantine_trials/data/temps.xlsx', sheet = 1)
sh2 = read_xlsx('post_quarantine_trials/data/temps.xlsx', sheet = 2)

head(sh1)
head(sh2)

df <- left_join(sh1, sh2)

df %>% 
  mutate(x = x + .3 * (as.numeric(as.factor(side)) -1)) %>%
  ggplot(., aes(x, y, color = temp )) +
  geom_point(size = 4, shape = 'square') +
  scale_color_viridis_c()

hist(df$temp)
