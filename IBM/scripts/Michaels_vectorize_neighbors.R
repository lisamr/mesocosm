#code from michael to vectorize which neighbors are infected
rm(list=ls())

#generate df with IDs and states
set.seed(1)
d <- data.frame(ID=1:10, state=sample(c(0,1), 10, replace = T))

#generate a neighbor matrix which could prob be done however you want if your neigh IDs dont change. you could also skip this step and go straight to having a df where each row is a pair of neighbors, with each col representing ID in neighbor pair
neighbors <- matrix(nrow = 10, ncol = 10, data=sample(c(F, T), size=100, replace = T))

#this is what were after. each row represents a pair of neighbors. you can subset in either direction using either rows or cols
neighbors <- as.data.frame(which(neighbors, arr.ind = T))
names(neighbors) <- c('neigh_1', 'neigh_2')

#count how many neighbors each ID has. ex/use neigh_1 col and get counts for how many total neighs each ID has
d$count_neighs <- as.data.frame(table(neighbors$neigh_1))[,2]
d

#now into the neighbor col, we can get the state for each neigh2 entry
neighbors$neigh_2_state <- d$state[neighbors$neigh_2]
neighbors

#finally we can say take the neighbors df, calculate sum of neigh_2_state values for each neigh_1 entry. so then we will have a count of how many infected neighs each neigh1 entry has, which will be sorted by neigh_1 ID by default. we turn that into a numeric vector and add that to the d df. 
d$count_inf_neighs <- as.numeric(with(neighbors, by(neigh_2_state, neigh_1, sum)))
d

#now we have a df that has every individual ID, their state, number of neighbors and number of infected neighbors. you can repeat the proces by calculating number of infected neighbors for each time step, and leave the rest of the neighs df unchanged, since the neighbor pairs dont change thru time. 