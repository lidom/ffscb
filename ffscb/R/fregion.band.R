#' Makes confidence bands
#'
#' @param x Functional parameter estimate. It can be either a vector or \link{fd} object from \link{fda}.
#' @param cov N * Cov(X), in which X is the functional estimator. It can be either matrix or \link{bifd} object from \link{fda}. The eigen decomposition of Cov(X) can be used instead.
#' @param tau Pointwise standard deviation of the standardized and differentiated sample functions. Can be estimated by tau_fun().
#' @param t0 
#' @param N It should be '1' if 'cov' is the covariance operator for X itself, which is the default value.
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
#' @param grid.size This determines on how fine grid the bands will be constructed before converted as an `fd' object. This parameter is used only when 'x' is fd object and 'cov' is bifd object.
#' @param Bs.sim.size This determines bootstrap sample size for Bs
#' @param n_int Number of intervals for the piecewise linear confidence bounds.
#' @param tol Controls the tolerance value used by stats::optimize(). The default (tol=NULL) leads to the functions' default values.
#' @return fregion.band Either a collection of vector valued bands or `fd' object whose objectname is changed to fregion.band.
#' @name band
#' @references 
#' \itemize{
#'    \item Liebl, D. and Reimherr, M. (2019). Fast and fair simultaneous confidence bands.
#'    \item  Choi, H. and Reimherr, M. (2018). A geometric approach to confidence regions and bands for functional parameters. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 80 239-260.
#' }
#' @examples
#' # Generate a sample
#' p <- 200 ; N <- 80 ; rangeval = c(0,1)
#' grid  <- make.grid(p, rangevals=rangeval)
#' mu0   <- meanf.poly(grid,c(0,1)) ; names(mu0) = grid
#' mu    <- meanf.poly(grid,c(0,1.1)) ; names(mu) = grid
#' cov.m <- make.cov.m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2/2,1,1))
#' x     <- make.sample(mu,cov.m,N)
#'
#' # Find the estimate and covariance
#' hat.mu       <- rowMeans(x)
#' hat.cov.m    <- crossprod(t(x - hat.mu)) / (N-1)
#' hat.tau.v    <- tau_fun(x)
#'
#' # Make and visualize/compare confidence bands
#' b <- fregion.band(x=hat.mu,cov=hat.cov.m,tau=hat.tau.v,N=N,
#'                   type=c("FFSCB.t", "Bs","BEc","naive.t"),conf.level=c(0.95))
#' plot(b)
#' @export 
fregion.band <- function(x, 
                         cov, 
                         tau         = NULL, 
                         t0          = NULL, 
                         N           = 1, 
                         type        = c("FFSCB.z", "FFSCB.t", "KR.z", "KR.t", "BEc", "Bs", "naive.t"), 
                         conf.level  = c(0.95), 
                         grid.size   = 200,
                         Bs.sim.size = 10000, 
                         n_int       = 10, 
                         tol         = NULL){
  ### Check the data type ###
  if (inherits(x,"fd") & (inherits(cov,"bifd") | inherits(cov,"pca.fd") | inherits(cov,"eigen.fd"))) datatype="fd" else if
     ((inherits(x,"numeric") | inherits(x,"matrix"))  & (inherits(cov,"matrix") | inherits(cov,"list") | inherits(cov,"eigen") )) datatype="vector" else stop ("The format of data is unknown")

  ### Evaluate x and cov ###
  # evaluate x and cov if datatype is "fd".
  # Since all functions for generating bands evaluate fd object inside, we do this here and just use vector/matrix version
  if (datatype=="fd") {
    if (!inherits(cov,"bifd")) {
      J    <- min(sum(cov$values > 0),dim(cov$harmonics$coefs)[2])
      coef <- cov$harmonics$coefs[,c(1:J)] %*% diag(cov$values[c(1:J)]) %*% t(cov$harmonics$coefs[,c(1:J)])
      cov  <- bifd(coef,cov$harmonics$basis,cov$harmonics$basis)
    }
    evalgrid <- make.grid(p=grid.size, rangevals=x$basis$rangeval)
    cov.m    <- eval.bifd(evalgrid, evalgrid, cov)
    x.v      <- eval.fd(evalgrid, x)
  } else {
    if (inherits(cov,"list")) {
      J     <- sum(cov$values > 0)
      cov.m <- cov$vectors[,c(1:J)] %*% diag(cov$values[c(1:J)]) %*% t(cov$vectors[,c(1:J)])
    } else cov.m <- cov;
    x.v <- x
  }
  p <- dim(cov.m)[1]
  if (!isSymmetric(cov.m)) cov.m <- (cov.m + t(cov.m))/2  # force cov.m to be symmetric

  ## Take eigen decomposition if BEc is used.
  if ("BEc" %in% type) {
    cor.m       <- cov2cor(cov.m)
    eigen.cov.m <- eigen(cov.m) ; eigen.cov.m$values[ eigen.cov.m$values < 0 ] <- 0 # trim negative eigenvalues.
    eigen.cor.m <- eigen(cor.m) ; eigen.cor.m$values[ eigen.cor.m$values < 0 ] <- 0 # trim negative eigenvalues.
  }
  
  if( sum(c("FFSCB.z", "FFSCB.t", "KR.z", "KR.t") %in% type) > 0 & is.null(tau)) {
    stop("The procedures FFSCB.z, FFSCB.t, KR.z, KR.t need tau.")
  }
  
  ## Parameter estimate in first column
  result           <- as.matrix(x.v,ncol=1) 
  colnames(result) <- c("x")
  
  ## Take loop for conf.level
  for (i in c(1:length(conf.level))){
    level <- conf.level[i]
    # Find number of fpc to use.
    if (!(level > 0 & level < 1)) stop("conf.level should have values between 0 and 1")

    # Make bands
    # If 'ALL' is included then run all the tests

    if ("Bs" %in% type) {
      tmp.colnames     <- c(colnames(result), paste0("Bs.u.",level), paste0("Bs.l.",level))
      Bs               <- make.band.Bs(cov=cov.m,conf.level=level,sim.size=Bs.sim.size) / sqrt(N)
      result           <- cbind(result, x.v + Bs, x.v - Bs);
      colnames(result) <- tmp.colnames
    }

    if ("BEc" %in% type) {
      tmp.colnames     <- c(colnames(result), paste0("BEc.u.",level), paste0("BEc.l.",level))
      BEc              <- make.band.BEc(eigen=eigen.cov.m, conf.level=level) / sqrt(N)
      result           <- cbind(result, x.v + BEc, x.v - BEc);
      colnames(result) <- tmp.colnames
    }

    if ("naive.t" %in% type) {
      tmp.colnames <- c(colnames(result), paste0("naive.t.u.",level), paste0("naive.t.l.",level))
      naive.t      <- make.band.naive.t(cov.m, conf.level=level, df=N-1) / sqrt(N)
      result       <- cbind(result, x.v + naive.t, x.v - naive.t);
      colnames(result) <- tmp.colnames
    }

    if ("KR.z" %in% type){
      tmp.colnames     <- c(colnames(result), paste0("KR.z.u.",level), paste0("KR.z.l.",level))
      KR.z             <- make.band.KR.z(tau=tau, conf.level=level) * sqrt(diag(cov.m)) / sqrt(N)
      result           <- cbind(result, x.v + KR.z, x.v - KR.z);
      colnames(result) <- tmp.colnames
    }

    if ("KR.t" %in% type){
      tmp.colnames     <- c(colnames(result), paste0("KR.t.u.",level), paste0("KR.t.l.",level))
      KR.t             <- make.band.KR.t(tau=tau, conf.level=level)*sqrt(diag(cov.m)) / sqrt(N)
      result           <- cbind(result, x.v + KR.t, x.v - KR.t);
      colnames(result) <- tmp.colnames
    }

    if ("FFSCB.z" %in% type){
      tmp.colnames     <- c(colnames(result), paste0("FFSCB.z.u.",level), paste0("FFSCB.z.l.",level))
      FFSCB.z          <- make.band.FFSCB.z(tau=tau, t0=t0, conf.level=level, n_int=n_int)*sqrt(diag(cov.m)) / sqrt(N)
      result           <- cbind(result, x.v + FFSCB.z, x.v - FFSCB.z);
      colnames(result) <- tmp.colnames
    }

    if ("FFSCB.t" %in% type){
      tmp.colnames     <- c(colnames(result), paste0("FFSCB.t.u.",level), paste0("FFSCB.t.l.",level))
      FFSCB.t          <- make.band.FFSCB.t(tau=tau, t0=t0, conf.level=level, N = N, n_int=n_int)*sqrt(diag(cov.m)) / sqrt(N)
      result           <- cbind(result, x.v + FFSCB.t, x.v - FFSCB.t);
      colnames(result) <- tmp.colnames
    }

  }
  if (datatype=="fd") {
    result.fd <- Data2fd(evalgrid, result, basisobj=x$basis)
    class(result.fd) <- "fregion.band"
    return(result.fd)
  } else {
    class(result) <- "fregion.band"
    return(result)
  }
}
