#' Takes eigen decomposition of bifd covariance object.
#'
#' @param cov.fd Covariance operator as a bifd object. Use \link{Data2bifd} to convert covariance matrix into a bifd object.
#' @export
eigen.fd <- function(cov.fd){
  BtB       <- inprod(cov.fd$sbasis,cov.fd$tbasis) ## sbisis and tbasis should be the same.
  B         <- chol(BtB)
  coefs     <- cov.fd$coefs
  coefsU    <- B%*%coefs%*%t(B)
  e.coefsU  <- eigen(coefsU)
  vcoefs    <- solve(B) %*% e.coefsU$vectors
  harmonics <- fd(vcoefs,cov.fd$sbasis)
  rtnobj    <- list(values=e.coefsU$values, harmonics=harmonics)
  ##
  class(rtnobj) <- "eigen.fd"
  return(rtnobj)
}
