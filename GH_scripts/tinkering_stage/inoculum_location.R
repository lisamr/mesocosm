#does location of poppy seed matter for inoculum load? I plated poppy seeds from 3 different zones of a PDA Rhizoc plate and plated them on water agar. There were 5 reps. After ~21 hours I counted # hyphae under compound scope.  

#difference not much, but probably prudent to collect from outer 1/3 of the plate.

rm(list=ls())
library(brms)
library(tidybayes)

#import data
zone <- 1:3
plate <- 1:5
hyphae <- c(15, 29, 7, 63, 18, 23, 16, 29, 43, 26, 30, 44, 16, 10, 47)
df <- cbind(expand.grid(zone=zone, plate=plate), hyphae)
df <- df %>% mutate_at(.vars = vars(zone, plate), .funs = as.factor)
head(df)

#viz data
ggplot(df, aes(plate, hyphae, color=zone)) +
  geom_point()

#model variation
f1 <- bf(hyphae ~  zone + (1|plate))
get_prior(f1, df)
rnorm(10000, 4, 1) %>% exp %>% density %>% plot(xlim=c(0,100)) #intercept
rnorm(10000, 0, 2.5) %>% exp %>% density %>% plot(xlim=c(0,100))#zone effect
p1 <- c(set_prior('normal(4,1)', class='Intercept'),
      set_prior('normal(0,2.5)', class='b'),
      set_prior('exponential(1)', class='sd'))
m1 <- brm(f1, data=df, family = poisson, prior = p1)
summary(m1)

#check out fit. looks good.
pp_check(m1)

#how do the zones vary? not by much. zone3 is a little higher, but not statistically speaking.
parnames(m1)
coefs <- m1 %>% 
  spread_draws(b_Intercept, b_zone2, b_zone3) %>%
  mutate(zone1 = b_Intercept,
         zone2 = b_Intercept + b_zone2,
         zone3 = b_Intercept + b_zone3,
         b_Intercept = NULL,
         b_zone2 = NULL,
         b_zone3 = NULL) %>% 
  gather_variables() 

ggplot(coefs, aes(x = .variable, y=.value)) +
  stat_pointinterval()

#how did the plates vary? any patterns? not really.
m1 %>% 
  gather_draws(r_plate[plate,]) %>% 
  ggplot(aes(x = plate, y=.value)) +
  stat_pointinterval()

