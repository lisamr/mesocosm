library(gdata)
#' A dispersal function that only has 3 values, one for local, one for adjacent
#' cells, and zero for others.  Requires distances between cells to = 1
#' @export
adjacent.dispersal <- function(distance, local, adjacent) {
  ifelse(distance < 0.5, local,
         ifelse(distance < 1.5, adjacent,0)
  )
}

#'Generate a lattice of equally spaced locations and coordinates
#'@export
MakeLattice <- function(nx, ny, dist=1) {
  locations <- cbind(location = 1:(nx*ny), 
                     x = rep(seq(from=0, by=dist, length.out=nx), each=ny),
                     y = rep(seq(from=0, by=dist, length.out=ny), times=nx)
  )
  return(locations)
}

DensityDependence.Rfun <- function(pop, space) { 
  1 - sum(pop*space)
}


require(compiler)
#' Produced multiplier for new recruits given population and space vectors
#'@import compiler
#'@export
DensityDependence <- cmpfun(DensityDependence.Rfun)

#'Runs the disease model.  Outputs a large matrix of population by species, ageclass, location
#'@import plyr 
#'@importFrom tidyr separate_
#'@export
SODModel <- function(parms.df, locations, time.steps, init, df.out=TRUE) {
  
  #Tests
  if(!all.equal(dim(init), c(nrow(locations), 2*nrow(parms.df)))) {
    stop("Dimensions of initial values and parameters do not match")
  }
  
  parms.df$class <- 1:nrow(parms.df)
  n.locations <- nrow(locations)
  
  #Convert parameter data frame into list of parameters
  parms.obj <- MakeParmsList(parms.df)
  #Create dispersal matrices
  spread.matrices <- MakeDispMatrix(parms.df, locations, parms.obj)  
  
  # Create transition matrix
  # TODO: Make this into a function
  
  #Unload list of parms into memory
  for(i in 1:length(parms.obj)) assign(names(parms.obj)[i], parms.obj[[i]])
  
  tran.mat <- matrix(0, nrow=n.classes*2, ncol=n.classes*2)
  diags <- row(tran.mat) - col(tran.mat)
  fec.mat <- tran.mat
  diag(tran.mat) <- 1 - trans.vec - mort.vec
  
  tran.mat[diags==1] <- diag(tran.mat)[-(nrow(tran.mat))] * rep(c(1,0),n.classes)[-(nrow(tran.mat))]
  tran.mat[diags==-1] <- diag(tran.mat)[-1] * rep(c(1,0),n.classes)[-(nrow(tran.mat))]
  
  tran.mat[diags==2] <- trans.vec[1:(n.classes*2 - 2)]
  tran.mat[diags==3] <- trans.vec[1:(n.classes*2 -3)] * rep(c(1,0),n.classes)[1:(n.classes*2 - 3)]
  tran.mat[diags==1] <- tran.mat[diags==1] + trans.vec[1:(n.classes*2-1)] * rep(c(0,1), n.classes)[1:(n.classes*2 - 1)]
  
  # Step update is (trans.mat * force.mat + fec.mat) * population
  
  # Create empty data matrix and populate with initial values
  
  pop <- array(NA, dim=c(length(time.steps), n.locations, 2*n.classes), 
               dimnames=list(Time=time.steps,
                             Location=1:n.locations, 
                             Class=paste(rep(names(classes),each=2),
                                         rep(c("S","I"),n.classes),sep=",")))
  
  
  
  pop[1,,] <- init
  spore.burden <- matrix(NA, nrow=n.classes, ncol=n.locations)
  #Rprof("out.prof")
  
  for(time in time.steps[-(length(time.steps))]) {
    
    # First act in simulation step.  Given population at each location, calculate spore burden at each location
    for(class in classes) {
      spore.burden[class,] <- pop[time,,class*2] %*% spread.matrices[class,,]
    }
    
    force.infection <- waifw %*% spore.burden
    real.recovery <- (1-force.infection) * recover
    
    for(location in 1:n.locations) {
      Force <- force.infection[,location]
      real.rec <- real.recovery[,location]
      force.matrix <- matrix(rbind(c(1-Force, Force), c(real.rec, 1-real.rec)), 2*n.classes, 2*n.classes, byrow=TRUE)
      
      E <- DensityDependence(pop[time,location,], space)
      for(i in 1:n.species) {
        classindex <- 1:(2*ageclasses[i]) + (sum(ageclasses[0:(i-1)])*2)
        fec.mat[classindex[1],classindex] <- (E*recruit.vec + resprout.vec*mort.vec)[classindex]  #Fecundities + Death X Resprout probabilties
      }
      trans.mat <- tran.mat*force.matrix + fec.mat
      pop[time + 1,location,] <- trans.mat %*% pop[time,location,]
    }
  }
  if(df.out) {
    pop <- reshape::melt(pop, value.name="Population")
    pop <- tidyr::separate_(pop, "Class", into=c("Species", "AgeClass", "Disease"), sep=",", remove=TRUE, convert=TRUE)
  }
  return(pop)
  
}


#'@import gdata
MakeParmsList <- function(treeparms.df) {

  classes <- 1:nrow(treeparms.df)
  names(classes) <- paste(treeparms.df$species, treeparms.df$ageclass, sep=",")
  
  ageclasses <- as.vector(table(treeparms.df$species))
  names(ageclasses) <- 1:length(unique(treeparms.df$species))
  
  parms.obj <- list(
    classes = classes,
    ageclasses = ageclasses,
    n.species = length(unique(treeparms.df$species)),
    n.classes = nrow(treeparms.df),
    waifw = as.matrix(treeparms.df[,matchcols(treeparms.df, 
                                              with="waifw[0-9]+"),],
                      dimnames = list(names(classes),names(classes))),
    recruit.vec = as.vector(rbind(treeparms.df$S.recruit, 
                                  treeparms.df$I.recruit)),
    mort.vec = as.vector(rbind(treeparms.df$S.mortality, 
                               treeparms.df$I.mortality)),
    trans.vec = as.vector(rbind(treeparms.df$S.transition, 
                                treeparms.df$I.transition)),
    resprout.vec = as.vector(rbind(treeparms.df$S.resprout, 
                                   treeparms.df$I.resprout)),
    recover = treeparms.df$recover,
    space = treeparms.df$space
  )
  return(parms.obj)
}




#'Generate a dispersal matrix from locations.
#'@import plyr
MakeDispMatrix <- function(parms.df, locations, parms.obj) {
  
  distance.matrix <- as.matrix(dist(locations[,c("x","y")]))
  
  # Generate dispersal matrices for each class
  spread.matrices <- plyr::daply(parms.df,"class", function(x) {
    do.call(x$kernel.fn, 
            unname(c(list(distance.matrix), 
                     as.list(na.omit(x[matchcols(x, "kernel.par[0-9]+")]))
            ))
    )
  })
  return(spread.matrices)
}

#install the data
library(RCurl)
d1=getURL("https://raw.githubusercontent.com/noamross/SODDr/master/inst/paper_tree_parms_eq.csv")
d2=getURL("https://raw.githubusercontent.com/noamross/SODDr/master/inst/paper_tree_parms_first.csv")
d3=getURL("https://raw.githubusercontent.com/noamross/SODDr/master/inst/paper_tree_parms_steady.csv")
d4=getURL("https://raw.githubusercontent.com/noamross/SODDr/master/inst/test_tree_parms.csv")

tree_parms_eq <- read.csv(text=d1)
tree_parms_first <- read.csv(text=d2)
tree_parms_steady <- read.csv(text=d3)
test_tree_parms <- read.csv(text=d4)

