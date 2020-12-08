library(tidyverse)
library(lubridate)
df <- read_csv('post_quarantine_trials/data/ibuttons/ibuttons_1026.csv')
theme_set(theme_bw())

str(df)
head(df)

df2 <- df %>% 
  rename(Time = `Date/Time`) %>% 
  mutate(Time = mdy_hm(Time), 
         Unit = NULL) 

p1 <- ggplot(df2, aes(Time, Value, group = button, color = button)) +
  geom_line()

plotly::ggplotly(p1)
