#' Performs hypothesis test(s) based on functional confidence region(s).
#'
#' @param x Functional parameter estimate. It can be either a vector or \link{fd} object from \link{fda}.
#' @param x0 Functional parameter under the null hypothesis. Zero function is assumed if it's not given.
#' @param cov N * Cov(X), in which X is the functional estimator. It can be either matrix or \link{bifd} object from \link{fda}. The eigen decomposition of Cov(X) can be used instead.
#' @param tau Pointwise standard deviation of the standardized and differentiated sample functions. Can be estimated by tau_fun().
#' @param t0 Parameter t0 of the fast and fair simultaneous confidence bands.
#' @param N It should be '1' if 'cov' is the covariance operator for X itself, which is the default value.
#' @param type This specifies which regions to be used for the tests. It should be a collection of the following: "FFSCB.t", "Ec".
#' \itemize{
#'   \item FFSCB.t : Test based on the fast 'n' fair simultaneous confidence band of Liebl and Reimherr (2019).
#'   \item Ec : Test based on the ellipsoid region in which radius are proportional to the square-root of corresponding eigenvalues as suggested in Choi and Reimherr (2018).
#' }
#' @param conf.level A vector of confidence levels for the bands to achieve.
#' @param n_int Number of intervals for the piecewise linear confidence bounds.
#' @param tol Controls the tolerance value used by stats::optimize(). The default (tol=NULL) leads to the functions' default values.
#' @param pc.cut It takes a vector of number of fPC to use in each HT. For integer values, fPC up to those values will be used. If it's a value from 0 to 1, this specifies the proportion of (estimated) variance that should be explained by the fPCs. If it is 0, all the available fPCs will be used as long as the size of eigenvalues are greater than .Machine$double.eps.
#' @param prec This determines the accuracy of \link{imhof}. One may try to modify this if p-value achieved in Ellipsoid form other than Epc gives negative value. It should the the form of c(epsabs, epsrel, limit).
#' @name testing
#' @references 
#' \itemize{
#'    \item Liebl, D. and Reimherr, M. (2019). Fast and fair simultaneous confidence bands.
#'    \item Choi, H. and Reimherr, M. (2018). A geometric approach to confidence regions and bands for functional parameters. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 80 239-260.
#' }
#' @examples
#' # Generate a sample
#' p <- 200 ; N <- 80 ; rangeval = c(0,1)
#' grid  <- make_grid(p, rangevals=rangeval)
#' mu0   <- meanf_poly(grid,c(0,1)) ; names(mu0) = grid
#' mu    <- meanf_poly(grid,c(0,1.1)) ; names(mu) = grid
#' cov.m <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2/2,1,1))
#' x     <- make_sample(mu,cov.m,N)
#'
#' # Find the estimate and covariance
#' hat.mu       <- rowMeans(x)
#' hat.cov.m    <- crossprod(t(x - hat.mu)) / (N-1)
#' hat.tau.v    <- tau_fun(x)
#'
#' # Compare different methods for Hypothesis testings.
#' fregion.test(x=hat.mu,x0=mu0,cov=hat.cov.m,tau=hat.tau.v,N=N,type=c("Ec"),
#' pc.cut=c(1,3,4,5,0.99,0.999))
#' @export
fregion.test <- function(x, 
                         x0          = 0, 
                         cov, 
                         tau         = NULL, 
                         t0          = NULL, 
                         N           = 1, 
                         type        = c("FFSCB.t", "Ec"), 
                         conf.level  = c(0.95), 
                         n_int       = 10, 
                         tol         = NULL,
                         pc.cut      = c(0.99), 
                         prec        = NULL) {
  ### Check the data type ###
  if (inherits(x,"fd") & (inherits(cov,"bifd") | inherits(cov,"pca.fd") | inherits(cov,"eigen.fd"))) datatype="fd" else if
  ((inherits(x,"numeric") | inherits(x,"matrix"))  & (inherits(cov,"matrix") | inherits(cov,"list") | inherits(cov,"eigen") )) datatype="vector" else stop ("The format of data is unknown")
  
  if ( sum(c("all", "All", "ALL") %in% type) > 0 ) type <- c("Ec","FFSCB.t")
  
  if( sum(c("FFSCB.z", "FFSCB.t", "KR.z", "KR.t") %in% type) > 0 & is.null(tau)) {
    stop("The procedures FFSCB.z, FFSCB.t, KR.z, KR.t need tau.")
  }
  
  ### If covariance is given as it is (not eigen decomposition), take eigen decomposition.
  e.cov <- cov
  if (inherits(cov,"bifd"))   e.cov <- eigen.fd(cov)
  if (inherits(cov,"matrix")) e.cov <- eigen(cov)
  
  # Trim cov to have all non-negative eigenvalues
  e.cov$values[e.cov$values < .Machine$double.eps] <- 0
  
  # Center the function x
  x_adj_null <- x - x0  # this works both for "fd" and vector.
  
  ################ Take loop for fpc.cut ##############
  result_all        <- vector("list",length = length(type))
  names(result_all) <- type
  ##
  if ("Ec" %in% type) {
    result <- NULL
    for (i in c(1:length(pc.cut))){
      cut <- pc.cut[i]
      # 5. Find number of fpc to use.
      if (cut < 0) stop("fpc.cut should be some positive fraction from 0 to 1, or integer greater than or equal to 1, or just 0 (to use all available pcs)")
      if (cut == 0 ) cut <- sum(e.cov$values > .Machine$double.eps) else { # if fpc.cut is 0, use all PCs
        if (cut < 1) cut <- get.req.n.pc(cut,e.cov$values) else {
          if (cut != round(cut)) {cut=round(cut) ; print("fpc.cut[",i,"] was rounded to the closest integer")}
        }
      }  # value 'Inf' can survive
      #####################################################
      ##
      pval.Ec             <- ffscb::get.pval.Ec(x=x_adj_null,N=N,eigen=e.cov,fpc.cut=cut,prec=prec)
      ##
      names.pval          <- ls(pattern=utils::glob2rx("pval.*"))
      pvalues             <- sapply(names.pval,get,inherits=FALSE,envir=environment()) #,envir=1
      names(pvalues)      <- sub("pval.","",names(pvalues))
      pvalues             <- c(pc.cut[i],cut,pvalues)
      result              <- rbind(result,pvalues)
      rownames(result)[i] <- i
    }
    colnames(result)[c(1,2)] <- c("pc.cut","pc.used")
    result_all[['Ec']] <- result
  }
  if ("FFSCB.t" %in% type) {
    if(length(x0)==1){ rep(x0, length(tau)) }
    diag.cov <- diag(cov)
    result_all[['FFSCB.t']] <- ffscb::get.pvalue.FFSCB.t(x=x,x0=x0,tau=tau,t0=t0,diag.cov=diag.cov,N=N,n_int=n_int)
  }
  return(result_all)
}
