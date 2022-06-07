#' ffscb
#'
#' This package contains example codes for inferential tools suggested in Liebl and Reimherr (2019), titled 'Fast and Fair Simultaneous Confidence Bands for Functional Parameters'.
#' The top level functions allow to construct simultaneous confidence bands and pvalue functions.
#' The inputs can be either vector/matrix or \link{fd}/\link{bifd} object from \link{fda} package.
#' @references Liebl, D. and Reimherr, M. (2019). Fast and fair simultaneous confidence bands for functional parameters.
#' @examples
#' # Generate a sample
#' p          <- 200 
#' N          <- 80 
#' grid       <- make_grid(p, rangevals=c(0,1))
#' mu0        <- meanf_poly(grid,c(0,1))   
#' names(mu0) <- grid
#' mu         <- meanf_poly(grid,c(0,1.1)) 
#' names(mu)  <- grid
#' cov.m      <- make_cov_m(cov.f = covf_st_matern, grid=grid, cov.f.params=c(2/2,1))
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
#' @docType package
#' @name ffscb
NULL
