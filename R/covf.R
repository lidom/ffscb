#' Matern Covariance Function
#'
#' @param x1 First argument of cov(x1,x2).
#' @param x2 Second argument of cov(x1,x2).
#' @param params Matern covariance function parameters: params=c(nu, sigma). 
#' @export
covf_st_matern <- function(x1, x2, params = c(1,1)){
  nu    <- params[1]
  sigma <- params[2]  
  l     <- 1
  d     <- sqrt(2*nu)*abs(x1-x2)/l
  if (d>0) {sigma^2 * 2^(1-nu) / gamma(nu) * d^nu * besselK(d,nu)} else {sigma^2}
}


#' Modified Matern Covariance Function (varying roughness parameter)
#'
#' @param x1 First argument of cov(x1,x2). Caution: It is assumed that 0<=x1<=1.
#' @param x2 Second argument of cov(x1,x2). Caution: It is assumed that 0<=x2<=1.
#' @param params Covariance function parameters: params=c(nu1, nu2, sigma). 
#' @export
covf_nonst_matern <- function(x1,x2,params=c(3/2, 1/2, 1)){
  nu    <- params[1] + sqrt(max(x1,x2)) * (params[2] - params[1])
  sigma <- params[3]  
  l     <- 1 
  d     <- sqrt(2*nu)*abs(x1-x2)/l
  if (d>0) {sigma^2 * 2^(1-nu) / gamma(nu) * d^nu * besselK(d,nu)} else {sigma^2}
}


# Matern Covariance Function with power-fct warped time-scale
#
# @param x1 First argument of cov(x1,x2).
# @param x2 Second argument of cov(x1,x2).
# @param params Matern covariance function parameters: params=c(nu, l, sigma, exponent). 
# @export
# covf_st_matern.warp.power <- function(x1,x2,params=c(1/2,1,1,2)){ #params <- c(nu,l,sigma,power)
#   covf_st_matern(warpf.power(x1,params[4]),warpf.power(x2,params[4]),params=params[1:3])
# }


# Matern Covariance Function with inverse power-fct warped time-scale
#
# @param x1 First argument of cov(x1,x2).
# @param x2 Second argument of cov(x1,x2).
# @param params Matern covariance function parameters: params=c(nu, l, sigma, exponent). 
# @export
# covf_st_matern.warp.power_inv <- function(x1,x2,params=c(1/2,1,1,2)){
#   if(any(c(x1,x2)<0) | any(c(x1,x2)>1)){stop("'warp.power_inv' is only for 0<= x1,x2 <=1")}
#   covf_st_matern(warpf.power_inv(x1,params[4]),warpf.power_inv(x2,params[4]),params=params[1:3])
# }


# Matern Covariance Function with inverse logit warped time-scale
#
# @param x1 First argument of cov(x1,x2).
# @param x2 Second argument of cov(x1,x2).
# @param params Matern covariance function parameters: params=c(nu, l, sigma). 
# @export
# covf_st_matern.warp.logit_inv <- function(x1,x2,params=c(1/2,1,1)){ #params <- c(nu,l,sigma,power)
#   covf_st_matern(warpf.logit_inv(x1),warpf.logit_inv(x2),params=params[1:3])
# }

# Matern Covariance Function with sigmoid warped time-scale
#
# @param x1 First argument of cov(x1,x2).
# @param x2 Second argument of cov(x1,x2).
# @param params Matern covariance function parameters: params=c(nu, l, sigma). 
# @param rangeval Range of evaluation points (default: rangeval=c(0,1)).
# @export
# covf_st_matern.warp.sigmoid <- function(x1,x2,params=c(1/2,1,1),rangeval=c(0,1)){ #params <- c(nu,l,sigma,power)
#   covf_st_matern(warpf.sigmoid(x1,rangeval),warpf.sigmoid(x2,rangeval),params=params)
# }


# warpf.power     <- function(t,params=c(10)){ t^params[1] }
# warpf.power_inv <- function(t,params=c(2)){-((-(t-1))^params[1]-1)}
# warpf.logit_inv <- function(t){ exp(t)/(1+exp(t)) }
# warpf.sigmoid   <- function(t, rangeval=c(0,1)){
#   t <- (t-rangeval[1])/(rangeval[2]-rangeval[1]) # to [ 0,1]
#   t <- (t-.5)*7                                 # to [-5,5]
#   exp(t)/(1+exp(t))
# }