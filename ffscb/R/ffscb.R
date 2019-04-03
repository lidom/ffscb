#' ffscb
#'
#' This package contains example codes for inferential tools suggested in Liebl and Reimherr (2019), titled 'Fast and Fair Simultaneous Confidence Bands'.
#' The top level functions allow to construct simultaneous confidence bands and pvalue functions.
#' The inputs can be either vector/matrix or \link{fd}/\link{bifd} object from \link{fda} package.
#' @references Liebl, D. and Reimherr, M. (2019). Fast and fair simultaneous confidence bands.
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
#' # Make and visualize/compare confidence bands
#' b <- confidence_band(x=hat.mu,cov=hat.cov.m,tau=hat.tau.v,N=N,
#'                   type=c("FFSCB.t", "Bs","BEc","naive.t"),conf.level=c(0.95))
#' plot(b)
#' @docType package
#' @name ffscb
NULL
