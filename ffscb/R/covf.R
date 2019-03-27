#' Matern Covariance Function
#'
#' @param x1 First argument of cov(x1,x2).
#' @param x2 Second argument of cov(x1,x2).
#' @param params Matern covariance function parameters: params=c(nu, l, sigma). 
#' @export
covf.st.matern <- function(x1,x2,params=c(1/2,1,1)){
  nu <- params[1] ; l <- params[2] ; sigma <- params[3]  #params=c(nu,l,sigma)
  d  <- sqrt(2*nu)*abs(x1-x2)/l
  if (d>0) {sigma^2 * 2^(1-nu) / gamma(nu) * d^nu * besselK(d,nu)} else {sigma^2}
}

#' Matern Covariance Function with power-fct warped time-scale
#'
#' @param x1 First argument of cov(x1,x2).
#' @param x2 Second argument of cov(x1,x2).
#' @param params Matern covariance function parameters: params=c(nu, l, sigma, exponent). 
#' @export
covf.st.matern.warp.power <- function(x1,x2,params=c(1/2,1,1,2)){ #params <- c(nu,l,sigma,power)
  covf.st.matern(ffscb:::warpf.power(x1,params[4]),ffscb:::warpf.power(x2,params[4]),params=params[1:3])
}


#' Matern Covariance Function with inverse power-fct warped time-scale
#'
#' @param x1 First argument of cov(x1,x2).
#' @param x2 Second argument of cov(x1,x2).
#' @param params Matern covariance function parameters: params=c(nu, l, sigma, exponent). 
#' @export
covf.st.matern.warp.power_inv <- function(x1,x2,params=c(1/2,1,1,2)){
  if(any(c(x1,x2)<0) | any(c(x1,x2)>1)){stop("'warp.power_inv' is only for 0<= x1,x2 <=1")}
  covf.st.matern(ffscb:::warpf.power_inv(x1,params[4]),ffscb:::warpf.power_inv(x2,params[4]),params=params[1:3])
}


#' Matern Covariance Function with inverse logit warped time-scale
#'
#' @param x1 First argument of cov(x1,x2).
#' @param x2 Second argument of cov(x1,x2).
#' @param params Matern covariance function parameters: params=c(nu, l, sigma). 
#' @export
covf.st.matern.warp.logit_inv <- function(x1,x2,params=c(1/2,1,1)){ #params <- c(nu,l,sigma,power)
  covf.st.matern(ffscb:::warpf.logit_inv(x1),ffscb:::warpf.logit_inv(x2),params=params[1:3])
}

#' Matern Covariance Function with sigmoid warped time-scale
#'
#' @param x1 First argument of cov(x1,x2).
#' @param x2 Second argument of cov(x1,x2).
#' @param params Matern covariance function parameters: params=c(nu, l, sigma). 
#' @export
covf.st.matern.warp.sigmoid <- function(x1,x2,params=c(1/2,1,1),rangeval=c(0,1)){ #params <- c(nu,l,sigma,power)
  covf.st.matern(ffscb:::warpf.sigmoid(x1,rangeval),ffscb:::warpf.sigmoid(x2,rangeval),params=params)
}


warpf.power     <- function(t,params=c(10)){ t^params[1] }
warpf.power_inv <- function(t,params=c(2)){-((-(t-1))^params[1]-1)}
warpf.logit_inv <- function(t){ exp(t)/(1+exp(t)) }
warpf.sigmoid   <- function(t, rangeval=c(0,1)){
  t <- (t-rangeval[1])/(rangeval[2]-rangeval[1]) # to [ 0,1]
  t <- (t-.5)*7                                 # to [-5,5]
  exp(t)/(1+exp(t))
}