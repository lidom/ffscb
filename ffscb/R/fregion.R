#' ffscb
#'
#' This package contains example codes for inferential tools suggested in \url{http://arxiv.org/abs/1607.07771}, titled 'Geometric Approach to Confidence Regions and Bands for Functional Parameters'.
#' The top level functions below take functional estimate and the covariance of it to perform hypothesis testings, to construct confidence bands, and to construct/visuallize hyper-rectangular regions.
#' The inputs can be either vector/matrix or \link{fd}/\link{bifd} object from \link{fda} package.
#' \itemize{
#'   \item \link{fregion.test}: Perform hypothesis testings based on confidence regions.
#'   \item \link{fregion.band}: Construct (simultaneous) confidence bands using hyper-ellipsoid confidence regions.
#' }
#' The example below is given in a mean function estimation context but works for other functional estimates as long as the estimators are Gaussian.
#'
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
#' # Compare different methods for Hypothesis testings.
#' (a1 <- fregion.test(x=hat.mu,x0=mu0,cov=hat.cov.m,tau=hat.tau.v,N=N,type=c("Ec"),pc.cut=c(1,3,4,5,0.99,0.999)))
#'
#' # Make and visualize/compare confidence bands
#' b <- fregion.band(x=hat.mu,cov=hat.cov.m,tau=hat.tau.v,N=N,
#'                   type=c("FFSCB.t", "Bs","BEc","naive.t"),conf.level=c(0.95))
#' plot(b)
#' @docType package
#' @name ffscb
NULL
