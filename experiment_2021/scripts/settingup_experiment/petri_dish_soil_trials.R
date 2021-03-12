#petri dish soil trials: testing the effect of soil on distance decay patterns in rhizoctonia

library(tidyverse)
df <- expand_grid(
  soil = c('sand', 'coirlite'),
  distance = c(2,3,4),
  day = 2*1:4, #day sampled, new target seed added
  rep = 1:10,
  present = NA
) %>% 
  pivot_wider(values_from = present, names_from = day, names_prefix = 'day') %>% 
  mutate(ID = 1:nrow(.)) %>% 
  select(ID, everything())

write_csv(df, 'post_quarantine_trials/output/petri_dish_soil_trials.csv', na = '')
