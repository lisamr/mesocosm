order <- seq(1:10)
test <- function(x, col) {100 + (-5)*col[x] + round(rnorm(1, sd = 3))}
test2 <- function(x, col) {100 + col[x] + round(rnorm(1, sd = 3))}
species <- letters[1:10]
abund <- sapply(1:10, function(x) test(x, order))
comp <- sapply(1:10, function(x) test2(x, abund))
df <- data.frame(species, order, abund, comp)
head(df)

ggplot(df, aes(order, abund))+
  geom_point()+
  geom_smooth(method = 'lm')
ggplot(df, aes(abund, comp))+
  geom_point()+
  geom_smooth(method = 'lm')
ggplot(df, aes(order, comp))+
  geom_point()+
  geom_smooth(method = 'lm')

#scramble relationship between abund and comp
df$abund.r <- sample(abund)
sort(order(df$abund.r ))
ggplot(df, aes(order.r, abund.r))+
  geom_point()+
  geom_smooth(method = 'lm')
ggplot(df, aes(abund.r, comp))+
  geom_point()+
  geom_smooth(method = 'lm')
ggplot(df, aes(order.r, comp))+
  geom_point()+
  geom_smooth(method = 'lm')
