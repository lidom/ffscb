#' Performs hypothesis test(s) based on functional confidence region(s).
#'
#' @param x Functional parameter estimate. It can be either a vector or \link{fd} object from \link{fda}.
#' @param x0 Functional parameter under the null hypothesis. Zero function is assumed if it's not given.
#' @param cov N * Cov(X), in which X is the functional estimator. It can be either matrix or \link{bifd} object from \link{fda}. The eigen decomposition of Cov(X) can be used instead.
#' @param N It should be '1' if 'cov' is the covariance operator for X itself, which is the default value.
#' @param type This specifies which regions to be used for the tests. It should be a collection of the following: "Enorm", "Epc", "Ec", "Ec1", "Rz", "Rz1", "Rzs", or "Rz1s".
#' \itemize{
#'   \item Enorm : The Hilbert space norm based test, which uses a `ball' in the function space(H).
#'   \item Epc : fPCA based test, which is a finite-dimensional chi-square ellipse in H
#'   \item Ec : The suggested test based on the ellipsoid region in which radius are proportional to the square-root of corresponding eigenvalues.
#'   \item Ec1 : The second suggestion for ellipsoie region, in which radius are proportional to the square-root of tails sums of eigenvalues. Doesn't require any smoothness assumption.
#'   \item Rz : The suggested rectangular region based test.
#'   \item Rz1 : The second suggested rectangular region based test.
#'   \item Rzs : The small sample version of Rz.
#'   \item Rz1s : The small sample version of Rz1.
#' }
#' @param pc.cut It takes a vector of number of fPC to use in each HT. For integer values, fPC up to those values will be used. If it's a value from 0 to 1, this specifies the proportion of (estimated) variance that should be explained by the fPCs. If it is 0, all the available fPCs will be used as long as the size of eigenvalues are greater than .Machine$double.eps.
#' @param prec This determines the accuracy of \link{imhof}. One may try to modify this if p-value achieved in Ellipsoid form other than Epc gives negative value. It should the the form of c(epsabs, epsrel, limit).
#' @param hat.cov An optional estimated covariance operator, which will be used when 'cov' is given as the true covariance and small sample version(s) are used. Usually not needed.
#' @param df Degree of freedom to use in small sample versions.
#' @param conf.level A vector of confidence levels for the bands to achieve.
#' @param band_pen Controls smoothness of adaptive band.  Note adaptive band still has proper coverage, but a smoother band may be desireable by the user.
#' @examples
#' # 1. Vector/matrix version
#'
#' # Generate a sample
#' p = 200 ; N = 80 ; rangeval = c(0,1)
#' grid = make.grid(p, rangevals=rangeval)
#' mu0 = meanf.poly(grid,c(0,1)) ; names(mu0) = grid
#' mu = meanf.poly(grid,c(0,1.2)) ; names(mu) = grid
#' cov.m = make.cov.m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2/2,1,1))
#' x = make.sample(mu,cov.m,N)
#'
#' # Find the estimate and covariance
#' hat.mu = rowMeans(x)
#' hat.cov.m = crossprod(t(x - hat.mu)) / (N-1)
#' e.hat.cov.m = eigen(hat.cov.m)   # <- This is optional and can be provide into the functions instead of hat.cov.m below.
#'
#' # Compare different methods for Hypothesis testings.
#' (a1 <- fregion.test(hat.mu,mu0,hat.cov.m,N,type=c("ALL"),pc.cut=c(1,3,4,5,0.99,0.999)))
#'
#'
#' # 2. fd/bifd version
#'
#' # create basis, convert vector/matrix into fd/bifd objects.
#' require(fda)
#' nbasis <- round(p*0.9)
#' fd.basis <- create.bspline.basis(rangeval,nbasis)
#' mu0.fd <- Data2fd(names(mu0),mu0,fd.basis)
#' mu.fd <- Data2fd(names(mu),mu,fd.basis)
#' x.fd <- Data2fd(rownames(x),x,fd.basis)
#' hat.mu.fd <- mean.fd(x.fd)
#' hat.cov.fd <- var.fd(x.fd)
#' e.hat.cov.fd <- eigen.fd(hat.cov.fd)   # <- This is optional and can be provide into the functions instead of hat.cov.fd below.
#'
#' # Compare different methods for Hypothesis testings.
#' (a1.fd=fregion.test(hat.mu.fd,mu0.fd,hat.cov.fd,N,type=c("ALL"),pc.cut=c(1,3,4,5,0.99,0.999)))
#'
#' @export

fregion.test <- function(x, x0=0, cov, N=1, type=c("Ec"), pc.cut=c(0.99), prec=NULL, hat.cov=NULL, df=NULL, conf.level=c(0.95),band_pen=1){
  ### 1. Check the data type ###
  if (inherits(x,"fd") & (inherits(cov,"bifd") | inherits(cov,"pca.fd") | inherits(cov,"eigen.fd"))) datatype="fd" else if
  ((inherits(x,"numeric") | inherits(x,"matrix"))  & (inherits(cov,"matrix") | inherits(cov,"list") | inherits(cov,"eigen") )) datatype="vector" else stop ("The format of data is unknown")

  ### 2. If covariance is given as it is (not eigen decomposition), take eigen decomposition.
  e.cov <- cov
  if (inherits(cov,"bifd")) e.cov <- eigen.fd(cov)
  if (inherits(cov,"matrix")) e.cov <- eigen(cov)

  # Trim cov to have all non-negative eigenvalues
  e.cov$values[e.cov$values < .Machine$double.eps] <- 0

  # 4. Center the function x
  x <- x - x0  # this works both for "fd" and vector.

  ################ Take loop for fpc.cut ##############
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

    # 6. Run tests ; If 'ALL' is included then run all the tests
    if ( sum(c("all", "All", "ALL") %in% type) > 0 ) type <- c("Enorm","Epc","Ec","Ec1","Rz","Rz1","Rzs","Rz1s")
    if ("Enorm" %in% type)  pval.Enorm <- fregion::get.pval.Enorm(x=x,N=N,eigen=e.cov,fpc.cut=cut,prec=prec)
    if ("Epc" %in% type)    pval.Epc   <- fregion::get.pval.Epc(x=x,N=N,eigen=e.cov,fpc.cut=cut)
    if ("Ec" %in% type)     pval.Ec    <- fregion::get.pval.Ec(x=x,N=N,eigen=e.cov,fpc.cut=cut,prec=prec)
    if ("Ec1" %in% type)    pval.Ec1   <- fregion::get.pval.Ec1(x=x,N=N,eigen=e.cov,fpc.cut=cut,prec=prec)
    if ("Rz" %in% type)     pval.Rz    <- fregion::get.pval.Rz(x=x,N=N,eigen=e.cov,fpc.cut=cut)
    if ("Rz1" %in% type)    pval.Rz1   <- fregion::get.pval.Rz1(x=x,N=N,eigen=e.cov,fpc.cut=cut)
    if ("Rzs" %in% type)    pval.Rzs   <- fregion::get.pval.Rzs(x=x,N=N,eigen=e.cov,fpc.cut=cut,hat.cov=hat.cov,df=df)
    if ("Rz1s" %in% type)   pval.Rz1s  <- fregion::get.pval.Rz1s(x=x,N=N,eigen=e.cov,fpc.cut=cut,hat.cov=hat.cov,df=df)
    if ("LRt" %in% type)   pval.LRt  <- fregion::get.pval.LRt(x=x, N = N, eigen=e.cov,
      band_pen = band_pen, conf.level=conf.level, fd.eval.grid.size=200, df=NULL)

    names.pval <- ls(pattern=glob2rx("pval.*"))
    pvalues <- sapply(names.pval,get,inherits=FALSE,envir=environment()) #,envir=1
    names(pvalues) <- sub("pval.","",names(pvalues))
    pvalues <- c(pc.cut[i],cut,pvalues)
    result <- rbind(result,pvalues)
    rownames(result)[i] <- i
  }
  colnames(result)[c(1,2)] <- c("pc.cut","pc.used")
  return(result)
}
