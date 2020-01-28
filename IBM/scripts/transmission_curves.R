x <- seq(0,20, by=.1)
beta_curve <- function(x){
  y <- .2*exp(-3*(log(x/11))^2)
  plot(x, y, type='l', ylim=c(0,.25))
}


beta_curve(x)

xx <- seq(0,9, by=.1)
alpha_curve <- function(x){
  y <- .4*(1-.3)^x
  
  plot(x, y, type='l')
}

alpha_curve(xx)
