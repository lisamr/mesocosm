#visualize transmission curves alpha, beta, and distance decay

#transmission curves
#mean values approximated from Otten et al. (2003)
alpha_curve <- function(x) .4*(1-.3)^x
alpha_curve <- function(x) .4*exp(-.7*x)
beta_curve <- function(x, amp=.2, gamma=3, tq=11) {
  #amp=amplitude
  #gamma=variance of curve
  #tq=time of peak
  amp*exp(-gamma*(log(x/tq))^2)
  }


#distance decay in probability of secondary transmission from Kleczkowski et al. (1997). increase sigma to decay a bit faster.
dist_decay <- function(x, a=14.9, d=.36, tau=3.73, sigma=.00099){
  C <- a*d*exp(d*tau)
  y <- C*exp(-sigma*x^2)

  #standardizing to be between 0 and 1, relative to a distance of 0. 
  y0 <- C*exp(-sigma*0^2)
  y/y0
  
  #standardize value so that at 20 mm, value is 1. reason is that beta is derived from interplanting distances of 20mm and I want everything to be relative to that. 
  #y20 <- C*exp(-sigma*20^2)
  #y/y20
}

distance2 <- seq(1,50, by=1) #mm
plot(distance2, dist_decay(distance2, sigma = .002), type='l')
which(dist_decay(distance2)<.1)#past 53mm essentially zero

#plot it
time <- seq(0,25, by=1) #days
alpha_curve <- function(x) .4*(1-.3)^x
alpha_curve <- function(x) .4*exp(-.3*x)
plot(time, alpha_curve(time), type='l')

par(mfrow=c(1,2))
distance <- seq(0,50, by=1) #mm
plot(time, beta_curve(time, gamma = 2, tq=15), type='l') 
plot(distance, dist_decay(distance), type='l')
par(mfrow=c(1,1))



#didn't use the function below, but fun to see.
#combine the beta and distance decay curves to make a beta curve as a function of time and distance. 
beta_t_x <- function(time, distance){
  
  #rows=distances in mm, cols=time in days
  B <- beta_curve(time)
  D <- dist_decay(distance)
  mat1 <- matrix(rep(B, length(D)), nrow = length(D), byrow = T)
  mat2 <- mat1*D #matrix of beta(x,t)
  return(mat2)
}
betamatrix <- beta_t_x(time, distance)

#plot. different lines represent distance decay in .5 cm increments (0 to 5 cm away).
pal <- viridis::plasma(length(distance))
plot(NULL, xlim=c(0,max(time)), ylim=c(0,.4), xlab="time", ylab='beta(t,x)')
lapply(1:length(distance), function(x) lines(time, betamatrix[x,], lwd=2, col=pal[x]))

#try 3d plot
plot_ly(z=~t(betamatrix)) %>% 
  add_surface() %>% 
  layout(
    title = "beta(time, distance)",
    scene = list(
      yaxis = list(title = "time (sec)", autorange = "reversed"),
      xaxis = list(title = "distance (mm)", autorange = "reversed"),
      zaxis = list(title = "beta (rate of transmission)")
    ))




