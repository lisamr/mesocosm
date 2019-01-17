#figuring out how to randomly sample assembly orders with more coverage
sp <- 1:6

sample(sp, 6, F)

#sample so the first number is equally represented across all species. And the second number is not represented twice. 

#get all the combinations of the first 2 numbers
test <- (t(combn(6, 2)))
test <- rbind.data.frame(test, cbind(test[,2], test[,1]))

#shuffle the rows and take the first 2 of each group
test2 <- test[sample(nrow(test)),]
test2 %>% 
  group_by(V1) %>% 
  filter((unique(V1)))
  View
unique(test2$V1)
