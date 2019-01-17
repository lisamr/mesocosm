#whats the relationship between competency vs. order and order vs. abundance?
library(ggplot2)
library(dplyr)

comp <- c(1, .3, .2, .1, 0, 0)
or <- 1:6
df <- data.frame(or, comp, abund=dlnorm(1:6, 1))
plot(dlnorm(1:10, 0), type='o')
#comp vs order: competency decreases with order
ggplot(df, aes(or, comp))+
  geom_point()+
  geom_line()+
  labs(y="competency", x="assembly order")
#order vs abundance: later species are rarer
#pool comes from a lognormal distribution
ggplot(df, aes(or, abund))+
  geom_point()+
  geom_line()+
  labs(y="proportion=dlnorm(1:6, 1)", x="assembly order")
#comp vs abund: more common species are more competent
ggplot(df, aes(abund, comp))+
  geom_point()+
  geom_line()

n=10
c <- replicate(n, sample(df$comp))  
c2 <- as.vector(c)
or2 <- rep(or, n)
abund2 <- rep(df$abund, n)
rep <- rep(1:n, each=6)
df2 <- data.frame(rep, or2, c2, abund2)
head(df2)
df1 <- cbind(rep=n+1, df)
names(df1) <- names(df2)
df2 <- rbind(df1, df2)
ggplot(df2, aes(or2, c2, group=rep, color=rep==11))+
  geom_point()+
  geom_line()+
  scale_color_manual(values=c("red", "blue")) +
  theme_classic()

