#' Computes tau(t), the pointwise standard deviation of the standardized and differentiated sample functions
#' 
#' @param x Matrix of sample functions (nrow=p, ncol=n, p=number of discretization point, n=sample size).  
#' @return tau_t Pointwise standard deviation of the standardized and differentiated sample functions.
#' @example 
#' p         <- 200 
#' N         <- 50
#' rangeval  <- c(0,1)
#' grid      <- make.grid(p, rangevals=rangeval)
#' mu        <- meanf.poly(grid, params = c(0,0)) 
#' 
#' # Generate random functions using a stationary 
#' # covariance function (homogeneous roughness (HR))
#' cov.m = make.cov.m(cov.f = covf.st.matern, grid=grid, 
#' cov.f.params=c(2,2,2))
#' X_HR  <-  make.sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
#' 
#' # Generate random functions using non-stationary 
#' # covariance function (increasing roughness (IR))
#' cov.m = make.cov.m(cov.f = covf.st.matern.warp.power, grid=grid, 
#' cov.f.params=c(1, 1, 1, 3))
#' X_IR  <-  make.sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
#' 
#' # Estimate tau(t):
#' tau_HR <- tau_fun(X_HR)
#' tau_IR <- tau_fun(X_IR)
#' 
#' # Plot data and estimated tau() functions
#' par(mfrow=c(2,2))
#' matplot(x=grid, y=X_HR, type="l", main="Homogeneous Roughness", 
#' ylab="X(t)", xlab="")
#' matplot(x=grid, y=X_IR, type="l", main="Increasing Roughness",  
#' ylab="X(t)", xlab="")
#' plot(x=grid, y=tau_HR,  type="l", main="Homogeneous Roughness", 
#' ylab="tau(t)", xlab="", ylim=range(tau_HR, tau_IR))
#' plot(x=grid, y=tau_IR,  type="l", main="Increasing Roughness",  
#' ylab="tau(t)", xlab="", ylim=range(tau_HR, tau_IR))
#' par(mfrow=c(1,1))
#' @export tau_fun
tau_fun <- function(x){
  x_scl    <- t(apply(x, 1, scale)) # pointwise standardization (mean=0, sd=1)
  x_p      <- apply(x_scl, 2,
                    FUN=function(yy){     # differentiation
                      xx <- seq(0,1,len=length(yy))
                      fn <- stats::splinefun(x = xx, y = yy, method = "natural")
                      pracma::fderiv(f = fn, x = xx, n = 1, h = diff(xx)[1], method = "central")}) 
  tau_t    <- apply(x_p, 1, sd)           # pointwise sd
  return(tau_t)
}