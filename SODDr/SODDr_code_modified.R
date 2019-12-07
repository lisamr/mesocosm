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



#'@import gdata
#'Turns the treeparms.df dataframe into a list that I guess is more readable to the SODmodel.
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




#'Generate a dispersal matrix from locations.I don't think the parms.obj argument is ever used.
#'@import plyr
MakeDispMatrix <- function(parms.df, locations, parms.obj) {
  
  #generate matrix of pairwise distances between every cell
  distance.matrix <- as.matrix(dist(locations[,c("x","y")]))
  
  #make a function to generate dispersal among cells
  fkernel <- function(x){
    #get kernel pars as list
    k_pars <- x[matchcols(x, "kernel.par[0-9]+")]
    k_pars <- as.list(na.omit(k_pars))
    #make a list of lists
    lists <- unname(c(list(distance.matrix), k_pars))
    #apply the kernel function with the supplied data in list form
    result <- do.call(x$kernel.fn, lists)
    return(result)
  }

  #generate dispersal matrix for each species-age, called "class" and put into array
  parms.df$class <- 1:nrow(parms.df) 
  spread.matrices <- array(NA, dim=c(length(parms.df$class), dim(distance.matrix)))
  for(i in 1:length(parms.df$class)){
    spread.matrices[i,,] <-  fkernel(parms.df[parms.df$class==i,])
  }
  
  return(spread.matrices)
}

#'Runs the disease model.  Outputs a large matrix of population by species, ageclass, location
#'@import plyr 
#'@importFrom tidyr separate_
#'@export
SODModel <- function(parms.df, locations, time.steps, init, df.out=TRUE) {
  #parameters
  
  #Tests
  if(!all.equal(dim(init), c(nrow(locations), 2*nrow(parms.df)))) {
    stop("Dimensions of initial values and parameters do not match")
  }
  
  n.locations <- nrow(locations)
  
  #Convert parameter data frame into list of parameters
  parms.obj <- MakeParmsList(parms.df)
  #Create dispersal matrices
  spread.matrices <- MakeDispMatrix(parms.df, locations, parms.obj)  
  
  # Create transition matrix
  # TODO: Make this into a function
  
  #Unload list of parms into memory. parms are assigned as objects within the SODModel and can be directly called as a result.
  for(i in 1:length(parms.obj)) {
    assign(names(parms.obj)[i], parms.obj[[i]])
  }
  
  #transition matrix of class-diseases (square matrix where length(rows)=n.classes*2)
  #matrix detailing demographic changes (i.e. growth and death. recruitment is in another mat)
  tran.mat <- matrix(0, nrow=n.classes*2, ncol=n.classes*2)
  diags <- row(tran.mat) - col(tran.mat)
  fec.mat <- tran.mat
  
  #trans.vec and mort.vec are values for S and I of each class. mortality rates determined from survival analysis and uninfected plots. transitions are the growth rates between size classes. Tanoak is the only one with non-zero values.
  #calculate stasis for individuals staying in same disease class (S->S, I->I)
  diag(tran.mat) <- 1 - trans.vec - mort.vec 
  
  #calculate stasis for individuals going from S->I (infecteds omitted)
  tran.mat[diags==1] <- (diag(tran.mat)* rep(c(1,0),n.classes))[-length(diag(tran.mat))]
  
  #calculate stasis for individuals going from I->S (susceptibles omitted)
  tran.mat[diags==-1] <- diag(tran.mat)[-1] * rep(c(1,0),n.classes)[-(nrow(tran.mat))]
  
  #calculate growth for individuals going from S->S and I->I
  tran.mat[diags==2] <- trans.vec[1:(n.classes*2 - 2)]
  
  #calculate growth for individuals going from S->I
  tran.mat[diags==3] <- trans.vec[1:(n.classes*2 -3)] * rep(c(1,0),n.classes)[1:(n.classes*2 - 3)]
  
  #calculate growth for individuals going from I->S
  tran.mat[diags==1] <- tran.mat[diags==1] + trans.vec[1:(n.classes*2-1)] * rep(c(0,1), n.classes)[1:(n.classes*2 - 1)]
  
  # Step update is (trans.mat * force.mat + fec.mat) * population
  
  # Create empty data matrix and populate with initial values 
  #dims = [1:time.steps, 1:n.locations, 1:2*n.classes]. 
  #keeps track of every disease-class at each location and each time step.
  pop <- array(NA, dim=c(length(time.steps), n.locations, 2*n.classes), 
               dimnames=list(Time=time.steps,
                             Location=1:n.locations, 
                             Class=paste(rep(names(classes),each=2),
                                         rep(c("S","I"),n.classes),sep=",")))
  
  #initialize model. add pops for first time step.
  pop[1,,] <- init
  
  #matrix for calculating spore burden from each of the classes in all of the lcoations?
  spore.burden <- matrix(NA, nrow=n.classes, ncol=n.locations)
  #Rprof("out.prof")
  
  for(time in time.steps[-(length(time.steps))]) {
    
    # First act in simulation step.  Given population at each location, calculate spore burden at each location
    for(class in classes) {
      #infected population of class_i at that time %*% dispersal matrix for class_i
      #spore.burden is a matrix for every class and location with values of how much spore there is. gets updated and overwritten every time step.
      spore.burden[class,] <- pop[time,,class*2] %*% spread.matrices[class,,]
    }
    
    #FOI=betas x the amount of spores. ex/[6,6] %*% [6,400], get matrix of every class at every location detailing FOI. value will be higher if surrounded by lots of spores (high spore burden)
    force.infection <- waifw %*% spore.burden
    #High FOI reduces rate of recovery.
    real.recovery <- (1-force.infection) * recover
    
    for(location in 1:n.locations) { #at every location...
      Force <- force.infection[,location] #FOI at that location for all of the classes
      real.rec <- real.recovery[,location] #recover at that location for all of the classes
      force.matrix <- matrix(rbind(c(1-Force, Force), c(real.rec, 1-real.rec)), 2*n.classes, 2*n.classes, byrow=TRUE) #a square matrix with values for S->S (1-force), I->I (1-recovery), S->I (force), I->S (recovery)
      
      
      E <- DensityDependence(pop[time,location,], space)#coef. for DD, proportion of unoccupied space. assumed all species affected equally? Low E means strong competition and low recruitment.
      for(i in 1:n.species) { #for every species...
        classindex <- 1:(2*ageclasses[i]) + (sum(ageclasses[0:(i-1)])*2)#identify the class index. each species will have 2*n.ageclasses.
        fec.mat[classindex[1],classindex] <- (E*recruit.vec + resprout.vec*mort.vec)[classindex]  #D-D Fecundities + Death X Resprout probabilties. recruitment = DD-recruitment + resprouts from deads. values are only filled in for the susceptible disease class (all recruits start as S)
      }
      trans.mat <- tran.mat*force.matrix + fec.mat #update the transition matrix 
      pop[time + 1,location,] <- trans.mat %*% pop[time,location,] #at each location and time, multiply population by transition probabilities to get population.
    }
  }
  if(df.out) {
    pop <- reshape::melt(pop, value.name="Population")
    pop <- tidyr::separate_(pop, "Class", into=c("Species", "AgeClass", "Disease"), sep=",", remove=TRUE, convert=TRUE)
  }
  return(pop)
  
}


