#' Makes confidence bands
#'
#' @param x Functional parameter estimate (for instance, the empirical mean function). It can be either a vector or \link{fd} object from \link{fda}.
#' @param cov.x Cov(x), in which x is the functional estimator (for instance, the covariance function of the empirical mean function). It can be either matrix or \link{bifd} object from \link{fda}. The eigen decomposition of Cov(X) can be used instead.
#' @param tau Pointwise standard deviation of the standardized and differentiated sample functions. Can be estimated by tau_fun().
#' @param t0 Parameter t0 of the fast and fair simultaneous confidence bands.
#' @param df Degrees of freedom parameter for the t-distribution based bands 'FFSCB.t', 'KR.t', and 'naive.t'. If x is the empirical mean function, set df=n-1, where n denotes the sample size.
#' @param type The band(s) to be constructed.
#' \itemize{
#'   \item FFSCB.z : Fast'n'Fair (adaptive) simultaneous confidence band based for a Gaussian functional parameter estimate.
#'   \item FFSCB.t : Fast'n'Fair (adaptive) simultaneous confidence band based for a t-distributed functional parameter estimate.
#'   \item KR.z : The constant simultaneous confidence band based on the classical Kac-Rice (KR) formula for Gaussian random functions.
#'   \item KR.t : The constant simultaneous confidence band based on the classical Kac-Rice (KR) formula for t-distributed random functions.
#'   \item BEc : The suggested modified Scheffe style band from hyper-ellipsoie Ec, which uses up to the very last dimension.
#'   \item Bs : Parametric bootstrap simultaneous confidence band, similar to the one appeard in Degras(2011) (for comparision purpose)
#'   \item naive.t : A collection of point-wise t-intervals. (for comparision purpose)
#' }
#' @param conf.level A vector of confidence levels for the bands to achieve.
#' @param grid.size This determines on how fine grid the bands will be constructed before converted as an `fd' object. This parameter is used only when 'x' is fd object and 'cov.x' is bifd object.
#' @param Bs.sim.size This determines bootstrap sample size for Bs
#' @param n_int Number of intervals for the piecewise linear confidence bounds.
#' @param tol tolerance 'tol' parameter used by stats::uniroot(). If tol=NULL, we use tol=.Machine$double.eps^0.32 which increases the default accuracy used by uniroot().
#' @return confidence_band Either a collection of vector valued bands or `fd' object whose objectname is changed to confidence_band.
#' @references 
#' \itemize{
#'    \item Liebl, D. and Reimherr, M. (2019). Fast and fair simultaneous confidence bands.
#'    \item  Choi, H. and Reimherr, M. (2018). A geometric approach to confidence regions and bands for functional parameters. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 80 239-260.
#' }
#' @examples
#' # Generate a sample
#' p          <- 200 
#' N          <- 80 
#' grid       <- make_grid(p, rangevals=c(0,1))
#' mu0        <- meanf_poly(grid,c(0,1))   
#' names(mu0) <- grid
#' mu         <- meanf_poly(grid,c(0,1.1)) 
#' names(mu)  <- grid
#' cov.m      <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2/2,1))
#' sample     <- make_sample(mu,cov.m,N)
#'
#' # Compute the estimate and its covariance
#' hat.mu     <- rowMeans(sample)
#' hat.cov    <- crossprod(t(sample - hat.mu)) / N
#' hat.cov.mu <- hat.cov / N
#' 
#' # Compute the tau-parameter (for the KR- and FFSCB-bands)
#' hat.tau    <- tau_fun(sample)
#'
#' # Make and plot confidence bands
#' b <- confidence_band(x=hat.mu, cov.x=hat.cov.mu, tau=hat.tau, df=N-1,
#'                      type=c("FFSCB.t", "Bs","BEc","naive.t"),
#'                      conf.level  = 0.95)
#' plot(b)
#' lines(x=grid, y=mu0, lty=2)
#' @export 
confidence_band <- function(x, 
                            cov.x, 
                            tau         = NULL, 
                            t0          = NULL, 
                            df          = NULL, 
                            type        = c("FFSCB.z", "FFSCB.t", "KR.z", "KR.t", "BEc", "Bs", "naive.t"), 
                            conf.level  = 0.95, 
                            grid.size   = 200,
                            Bs.sim.size = 10000, 
                            n_int       = 10,
                            tol         = NULL){
  ### Check the data type ###
  if (inherits(x,"fd") & (inherits(cov.x,"bifd") | inherits(cov.x,"pca.fd") | inherits(cov.x,"eigen.fd"))) datatype="fd" else if
  ((inherits(x,"numeric") | inherits(x,"matrix"))  & (inherits(cov.x,"matrix") | inherits(cov.x,"list") | inherits(cov.x,"eigen") )) datatype="vector" else stop ("The format of data is unknown")
  
  ### Evaluate x and cov.x ###
  # evaluate x and cov.x if datatype is "fd".
  # Since all functions for generating bands evaluate fd object inside, we do this here and just use vector/matrix version
  if (datatype=="fd") {
    if (!inherits(cov.x,"bifd")) {
      J    <- min(sum(cov.x$values > 0),dim(cov.x$harmonics$coefs)[2])
      coef <- cov.x$harmonics$coefs[,c(1:J)] %*% diag(cov.x$values[c(1:J)]) %*% t(cov.x$harmonics$coefs[,c(1:J)])
      cov.x  <- fda::bifd(coef,cov.x$harmonics$basis,cov.x$harmonics$basis)
    }
    evalgrid <- make_grid(p=grid.size, rangevals=x$basis$rangeval)
    cov.m    <- fda::eval.bifd(evalgrid, evalgrid, cov.x)
    x.v      <- fda::eval.fd(evalgrid, x)
  } else {
    if (inherits(cov.x,"list")) {
      J     <- sum(cov.x$values > 0)
      cov.m <- cov.x$vectors[,c(1:J)] %*% diag(cov.x$values[c(1:J)]) %*% t(cov.x$vectors[,c(1:J)])
    } else cov.m <- cov.x;
    x.v <- x
  }
  p <- dim(cov.m)[1]
  if (!isSymmetric(cov.m)) cov.m <- (cov.m + t(cov.m))/2  # force cov.m to be symmetric
  
  ## Take eigen decomposition if BEc is used.
  if ("BEc" %in% type) {
    cor.m       <- stats::cov2cor(cov.m)
    eigen.cov.m <- eigen(cov.m) ; eigen.cov.m$values[ eigen.cov.m$values < 0 ] <- 0 # trim negative eigenvalues.
    eigen.cor.m <- eigen(cor.m) ; eigen.cor.m$values[ eigen.cor.m$values < 0 ] <- 0 # trim negative eigenvalues.
  }
  
  if( sum(c("FFSCB.z", "FFSCB.t", "KR.z", "KR.t") %in% type) > 0 & is.null(tau)) {
    stop("The procedures FFSCB.z, FFSCB.t, KR.z, and KR.t need a tau parameter.")
  }
  if( sum(c("FFSCB.t", "KR.t", "naive.t") %in% type) > 0 & is.null(df)) {
    stop("The procedures FFSCB.t, KR.t, and naive.t need a df parameter.")
  }
  
  ## Parameter estimate in first column
  result           <- as.matrix(x.v, ncol=1) 
  colnames(result) <- c("x")
  
  if(is.null(tol)){
    tol         <- .Machine$double.eps^0.32 # increases the default accuracy for uniroot() (.Machine$double.eps^0.25) used by the FFSCB functions
  }
  
  ## Take loop for conf.level
  for (i in c(1:length(conf.level))){
    level <- conf.level[i]
    # Find number of fpc to use.
    if (!(level > 0 & level < 1)) stop("conf.level should have values between 0 and 1")
    
    # Make bands
    # If 'ALL' is included then run all the tests
    
    if ("Bs" %in% type) {
      tmp.colnames     <- c(colnames(result), paste0("Bs.u.",level), paste0("Bs.l.",level))
      Bs               <- make_band_Bs(cov=cov.m, conf.level=level, sim.size=Bs.sim.size) 
      result           <- cbind(result, x.v + Bs, x.v - Bs)
      colnames(result) <- tmp.colnames
    }
    
    if ("BEc" %in% type) {
      tmp.colnames     <- c(colnames(result), paste0("BEc.u.",level), paste0("BEc.l.",level))
      BEc              <- make_band_Ec(eigen=eigen.cov.m, conf.level=level) 
      result           <- cbind(result, x.v + BEc, x.v - BEc)
      colnames(result) <- tmp.colnames
    }
    
    if ("naive.t" %in% type) {
      tmp.colnames     <- c(colnames(result), paste0("naive.t.u.",level), paste0("naive.t.l.",level))
      naive.t          <- make_band_naive_t(cov=cov.m, conf.level=level, df=df)
      result           <- cbind(result, x.v + naive.t, x.v - naive.t)
      colnames(result) <- tmp.colnames
    }
    
    if ("KR.z" %in% type){
      tmp.colnames     <- c(colnames(result), paste0("KR.z.u.",level), paste0("KR.z.l.",level))
      KR.z             <- make_band_KR_z(tau=tau, diag.cov=diag(cov.m), conf.level=level)
      result           <- cbind(result, x.v + KR.z, x.v - KR.z)
      colnames(result) <- tmp.colnames
    }
    
    if ("KR.t" %in% type){
      tmp.colnames     <- c(colnames(result), paste0("KR.t.u.",level), paste0("KR.t.l.",level))
      KR.t             <- make_band_KR_t(tau=tau, diag.cov=diag(cov.m), df=df, conf.level=level)
      result           <- cbind(result, x.v + KR.t, x.v - KR.t)
      colnames(result) <- tmp.colnames
    }
    
    if ("FFSCB.z" %in% type){
      tmp.colnames     <- c(colnames(result), paste0("FFSCB.z.u.",level), paste0("FFSCB.z.l.",level))
      FFSCB.z          <- .make_band_FFSCB_z(tau=tau, t0=t0, diag.cov=diag(cov.m), conf.level=level, n_int=n_int, tol=tol)
      result           <- cbind(result, x.v + FFSCB.z$band, x.v - FFSCB.z$band)
      colnames(result) <- tmp.colnames
    }
    
    if ("FFSCB.t" %in% type){
      tmp.colnames     <- c(colnames(result), paste0("FFSCB.t.u.",level), paste0("FFSCB.t.l.",level))
      if(df <= 100){
        FFSCB.t          <- .make_band_FFSCB_t(tau=tau, t0=t0, diag.cov=diag(cov.m), df=df, conf.level=level, n_int=n_int, tol=tol)
      }else{
        FFSCB.t          <- .make_band_FFSCB_z(tau=tau, t0=t0, diag.cov=diag(cov.m),        conf.level=level, n_int=n_int, tol=tol)
      }
      result           <- cbind(result, x.v + FFSCB.t$band, x.v - FFSCB.t$band)
      colnames(result) <- tmp.colnames
    }
    
  }
  if (datatype=="fd") {
    result.fd <- fda::Data2fd(evalgrid, result, basisobj=x$basis)
    class(result.fd) <- "confidence_band"
    return(result.fd)
  } else {
    class(result) <- "confidence_band"
    return(result)
  }
}




#' Makes confidence bands for fragmentary functional data
#'
#' @param x Functional parameter estimate (for instance, the empirical mean function). It can be either a vector or \link{fd} object from \link{fda}.
#' @param diag.cov.x diag(Cov(x)), in which x is the functional estimator (for instance, the covariance function of the empirical mean function). It can be either matrix or \link{bifd} object from \link{fda}. The eigen decomposition of Cov(X) can be used instead.
#' @param tau Pointwise standard deviation of the standardized and differentiated sample functions. Can be estimated by tau_fun().
#' @param t0 Parameter t0 of the fast and fair simultaneous confidence bands.
#' @param df Degrees of freedom parameter for the t-distribution based bands 'FFSCB.t', 'KR.t', and 'naive.t'. If x is the empirical mean function, set df=n-1, where n denotes the sample size.
#' @param type The band(s) to be constructed.
#' \itemize{
#'   \item FFSCB.z : Fast'n'Fair (adaptive) simultaneous confidence band based for a Gaussian functional parameter estimate.
#'   \item FFSCB.t : Fast'n'Fair (adaptive) simultaneous confidence band based for a t-distributed functional parameter estimate.
#'   \item KR.z : The constant simultaneous confidence band based on the classical Kac-Rice (KR) formula for Gaussian random functions.
#'   \item KR.t : The constant simultaneous confidence band based on the classical Kac-Rice (KR) formula for t-distributed random functions.
#' }
#' @param conf.level A vector of confidence levels for the bands to achieve.
#' @param n_int Number of intervals for the piecewise linear confidence bounds.
#' @param tol tolerance 'tol' parameter used by stats::uniroot(). If tol=NULL, we use tol=.Machine$double.eps^0.32 which increases the default accuracy used by uniroot().
#' @return confidence_band_fragm
#' @references 
#' \itemize{
#'    \item Liebl, D. and Reimherr, M. (2019). Fast and fair simultaneous confidence bands.
#' }
confidence_band_fragm <- function(x, 
                                  diag.cov.x, 
                                  tau         = NULL, 
                                  t0          = NULL, 
                                  df          = NULL, 
                                  type        = c("FFSCB.z", "FFSCB.t", "KR.z", "KR.t", "naive.t"), 
                                  conf.level  = 0.95, 
                                  n_int       = 10,
                                  tol         = NULL){
  ##
  
  if( sum(c("FFSCB.z", "FFSCB.t", "KR.z", "KR.t") %in% type) > 0 & is.null(tau)) {
    stop("The procedures FFSCB.z, FFSCB.t, KR.z, and KR.t need a tau parameter.")
  }
  if( sum(c("FFSCB.t", "KR.t", "naive.t") %in% type) > 0 & is.null(df)) {
    stop("The procedures FFSCB.t, KR.t, and naive.t need a df parameter.")
  }
  
  ## Parameter estimate in first column
  result           <- as.matrix(x, ncol=1) 
  colnames(result) <- c("x")
  
  if(is.null(tol)){
    tol         <- .Machine$double.eps^0.32 # increases the default accuracy for uniroot() (.Machine$double.eps^0.25) used by the FFSCB functions
  }
  
  ## Take loop for conf.level
  for (i in c(1:length(conf.level))){
    level <- conf.level[i]
    # Find number of fpc to use.
    if (!(level > 0 & level < 1)) stop("conf.level should have values between 0 and 1")
    
    # Make bands
    # If 'ALL' is included then run all the tests
    
    if ("naive.t" %in% type) {
      tmp.colnames     <- c(colnames(result), paste0("naive.t.u.",level), paste0("naive.t.l.",level))
      naive.t          <- make_band_naive_t_fragm(diag.cov=diag.cov.x, conf.level=level, df=df)
      result           <- cbind(result, x + naive.t, x - naive.t)
      colnames(result) <- tmp.colnames
    }
    
    if ("KR.z" %in% type){
      tmp.colnames     <- c(colnames(result), paste0("KR.z.u.",level), paste0("KR.z.l.",level))
      KR.z             <- make_band_KR_z(tau=tau, diag.cov=diag.cov.x, conf.level=level)
      result           <- cbind(result, x + KR.z, x - KR.z)
      colnames(result) <- tmp.colnames
    }
    
    if ("KR.t" %in% type){
      tmp.colnames     <- c(colnames(result), paste0("KR.t.u.",level), paste0("KR.t.l.",level))
      KR.t             <- make_band_KR_t(tau=tau, diag.cov=diag.cov.x, df=df, conf.level=level)
      result           <- cbind(result, x + KR.t, x - KR.t)
      colnames(result) <- tmp.colnames
    }
    
    if ("FFSCB.z" %in% type){
      tmp.colnames     <- c(colnames(result), paste0("FFSCB.z.u.",level), paste0("FFSCB.z.l.",level))
      FFSCB.z          <- .make_band_FFSCB_z(tau=tau, t0=t0, diag.cov=diag.cov.x, conf.level=level, n_int=n_int, tol=tol)
      result           <- cbind(result, x + FFSCB.z$band, x - FFSCB.z$band)
      colnames(result) <- tmp.colnames
    }
    
    if ("FFSCB.t" %in% type){
      tmp.colnames     <- c(colnames(result), paste0("FFSCB.t.u.",level), paste0("FFSCB.t.l.",level))
      if(df <= 100){
        FFSCB.t          <- .make_band_FFSCB_t(tau=tau, t0=t0, diag.cov=diag.cov.x, df=df, conf.level=level, n_int=n_int, tol=tol)
      }else{
        FFSCB.t          <- .make_band_FFSCB_z(tau=tau, t0=t0, diag.cov=diag.cov.x,        conf.level=level, n_int=n_int, tol=tol)
      }
      result           <- cbind(result, x + FFSCB.t$band, x - FFSCB.t$band)
      colnames(result) <- tmp.colnames
    }
  }
  class(result) <- "confidence_band_fragm"
  return(result)
}
