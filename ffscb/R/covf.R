#' @export
covf.st.matern <- function(x1,x2,params=c(1/2,1,1)){
  nu <- params[1] ; l <- params[2] ; sigma <- params[3]  #params=c(nu,l,sigma)
  d  <- sqrt(2*nu)*abs(x1-x2)/l
  if (d>0) {sigma^2 * 2^(1-nu) / gamma(nu) * d^nu * besselK(d,nu)} else {sigma^2}
}

#' @export
covf.st.matern.warp.power <- function(x1,x2,params=c(1/2,1,1,10)){ #params <- c(nu,l,sigma,power)
  covf.st.matern(warpf.power(x1,params[4]),warpf.power(x2,params[4]),params=params[1:3])
}


#' @export
covf.st.matern.warp.power_inv <- function(x1,x2,params=c(1/2,1,1,2)){
  if(any(c(x1,x2)<0) | any(c(x1,x2)>1)){stop("'warp.power_inv' is only for 0<= x1,x2 <=1")}
  covf.st.matern(warpf.power_inv(x1,params[4]),warpf.power_inv(x2,params[4]),params=params[1:3])
}


#' @export
covf.st.matern.warp.logit_inv <- function(x1,x2,params=c(1/2,1,1)){ #params <- c(nu,l,sigma,power)
  covf.st.matern(warpf.logit_inv(x1),warpf.logit_inv(x2),params=params[1:3])
}

#' @export
covf.st.matern.warp.sigmoid <- function(x1,x2,params=c(1/2,1,1),rangeval=c(0,1)){ #params <- c(nu,l,sigma,power)
  covf.st.matern(warpf.sigmoid(x1,rangeval),warpf.sigmoid(x2,rangeval),params=params)
}