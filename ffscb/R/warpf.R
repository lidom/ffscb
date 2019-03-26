#' @export
warpf.power     <- function(t,params=c(10)){ t^params[1] }

#' @export
warpf.power_inv <- function(t,params=c(2)){-((-(t-1))^params[1]-1)}

#' @export
warpf.logit_inv <- function(t){ exp(t)/(1+exp(t)) }

#' @export
warpf.sigmoid <- function(t, rangeval=c(0,1)){
  t <- (t-rangeval[1])/(rangeval[2]-rangeval[1]) # to [ 0,1]
  t <- (t-.5)*7                                 # to [-5,5]
  exp(t)/(1+exp(t))
}