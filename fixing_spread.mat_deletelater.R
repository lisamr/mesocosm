#trash
MakeDispMatrix <- function(parms.df, locations, parms.obj) {
  
  distance.matrix <- as.matrix(dist(locations[,c("x","y")]))
  
  # Generate dispersal matrices for each class
  #creates an array of matrices for each class
  spread.matrices <- plyr::daply(treeparms.df,"classes", function(x) {
    do.call(x$kernel.fn, 
            unname(c(list(distance.matrix), 
                     as.list(na.omit(x[matchcols(x, "kernel.par[0-9]+")]))
            ))
    )
  })
  return(spread.matrices)
}
#end trash

#got this to work. 
test <- treeparms.df %>% filter(ageclass%in%c(1,2,3))
fmatspread <- function(x){
  kernel.pars <- x[matchcols(x, 'kernel.par[0-9]+')] %>% 
    na.omit %>% as.list #get kernel parameters into a list
  dmat <- list(distance.matrix)
  arguments <- c(dmat, kernel.pars) %>% unname()
  #do.call('adjacent.dispersal', arguments)
  do.call(as.character(x$kernel.fn), arguments)
}
test2 <- plyr::daply(test, c('species', 'ageclass'), fmatspread)
str(test2)
