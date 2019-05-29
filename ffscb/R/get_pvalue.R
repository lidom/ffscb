#' p-value (ellipsoid region)
#'
#' @param x function argument
#' @param x0 Functional parameter under the null hypothesis. Zero function is assumed if it's not given.
#' @param eigen eigen decomposition of covariance function
#' @param fpc.cut It takes a vector of number of fPC to use in each HT. For integer values, fPC up to those values will be used. If it's a value from 0 to 1, this specifies the proportion of (estimated) variance that should be explained by the fPCs. If it is 0, all the available fPCs will be used as long as the size of eigenvalues are greater than .Machine$double.eps.
#' @param prec This determines the accuracy of \link{imhof}. One may try to modify this if p-value achieved in Ellipsoid form other than Epc gives negative value. It should the the form of c(epsabs, epsrel, limit).
#' @references Choi, H. and Reimherr, M. (2018). A geometric approach to confidence regions and bands for functional parameters. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 80 239-260.
#' @examples 
#' # Generate a sample
#' p <- 200 
#' N <- 80 
#' grid   <- make_grid(p, rangevals=c(0,1))
#' mu0    <- meanf_poly(grid,c(0,1))   ; names(mu0) <- grid
#' mu     <- meanf_poly(grid,c(0,1.1)) ; names(mu)  <- grid
#' cov.m  <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2/2,1,1))
#' sample <- make_sample(mu,cov.m,N)
#'
#' # Compute the estimate and its covariance
#' hat.mu     <- rowMeans(sample)
#' hat.cov    <- crossprod(t(sample - hat.mu)) / N
#' hat.cov.mu <- hat.cov / N
#' 
#' # Eigen decomposition
#' e.cov.mu <- eigen(hat.cov.mu)
#' 
#' # pvalue
#' pval <- get_pval_Ec(x=hat.mu, x0=mu0, eigen=e.cov.mu)
#' pval
#' @export
get_pval_Ec <- function(x, x0=NULL, eigen, fpc.cut=NULL, prec=NULL){
  if (is.null(prec)) {prec <- c(10^(-6),10^(-6),10000)}
  if (is.null(fpc.cut) ) {fpc.cut <- sum(eigen$values > .Machine$double.eps)}
  if (is.null(x0)) { x <- x - x0 }  # this works both for "fd" and vector.
  ##
  c.square <- sqrt(eigen$values[1:fpc.cut])
  if (inherits(x,"fd")) {coef.square <- fda::inprod(x, eigen$harmonics[1:fpc.cut])^2} else
                        {coef.square <- as.vector(crossprod(x, eigen$vectors[,1:fpc.cut])^2)}
  weights <- eigen$values[1:fpc.cut]/c.square
  stat    <- sum(coef.square / c.square)
  CompQuadForm::imhof(stat, weights,epsabs=prec[1],epsrel=prec[2],limit=prec[3])[[1]]
}


# get_pval_Ec <- function(x, x0=NULL, N = 1, eigen, fpc.cut=NULL, prec=NULL){
#   if (is.null(prec)) {prec <- c(10^(-6),10^(-6),10000)}
#   if (is.null(fpc.cut) | fpc.cut==Inf) {fpc.cut <- sum(eigen$values > .Machine$double.eps)}
#   if (is.null(x0)) { x <- x - x0 }  # this works both for "fd" and vector.
#   ##
#   c.square <- sqrt(eigen$values[1:fpc.cut])
#   if (inherits(x,"fd")) {coef.square <- fda::inprod(x, eigen$harmonics[1:fpc.cut])^2} else
#                         {coef.square <- as.vector(crossprod(x, eigen$vectors[,1:fpc.cut])^2)}
#   weights <- eigen$values[1:fpc.cut]/c.square/N
#   stat    <- sum(coef.square / c.square)
#   CompQuadForm::imhof(stat, weights,epsabs=prec[1],epsrel=prec[2],limit=prec[3])[[1]]
# }



#' FFSCB p-value (Gaussian)
#'
#' @param x Functional parameter estimate (for instance, the empirical mean function).
#' @param x0 Functional parameter under the null hypothesis. Default: zero.
#' @param tau Pointwise standard deviation of the standardized and differentiated sample functions. Can be estimated by tau_fun().
#' @param t0 Parameter t0 of the fast and fair simultaneous confidence bands.
#' @param diag.cov The diagonal of Cov(x), in which x is the functional estimator. For instance, the diagonal of the discretized covariance function of the empirical mean function x.  
#' @param n.eval.points Number of evaluation points for the p-value function. Large values (>10) lead to slow computations.
#' @param n_int Number of intervals parameter used by the function make_band_FFSCB_z()
#' @references Liebl, D. and Reimherr, M. (2019). Fast and fair simultaneous confidence bands.
#' @examples 
#' # Generate a sample
#' p <- 200 
#' N <- 80 
#' grid   <- make_grid(p, rangevals=c(0,1))
#' mu0    <- meanf_poly(grid,c(0,1))   ; names(mu0) <- grid
#' mu     <- meanf_poly(grid,c(0,1.1)) ; names(mu)  <- grid
#' cov.m  <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2/2,1,1))
#' sample <- make_sample(mu,cov.m,N)
#'
#' # Compute the estimate and its covariance
#' hat.mu     <- rowMeans(sample)
#' hat.cov    <- crossprod(t(sample - hat.mu)) / N
#' hat.cov.mu <- hat.cov / N
#' hat.tau    <- tau_fun(sample)
#' 
#' # pvalue
#' pval <- get_pvalue_FFSCB_z(x=hat.mu, x0=mu0, tau=hat.tau, 
#'                            diag.cov=diag(hat.cov.mu), 
#'                            n.eval.points=6)
#' plot(y=pval, x=grid, type="l", main="pvalue FFSCB (Gaussian)")
#' @export
get_pvalue_FFSCB_z <- function(x, x0=NULL, tau, t0=NULL, diag.cov, n.eval.points=11, n_int=5){
  if (is.null(x0)) {x0 <- rep(0,times=length(x))}
  t_grid        <- seq(0,1,len=n.eval.points)
  p_grid        <- numeric(n.eval.points)
  x_f           <- stats::splinefun(x = seq(0,1,len=length(tau)), y = x,  method = "natural")
  x0_f          <- stats::splinefun(x = seq(0,1,len=length(tau)), y = x0, method = "natural")
  diag.cov_f    <- stats::splinefun(x = seq(0,1,len=length(tau)), y = diag.cov, method = "natural")
  ##
  myfun <- function(p, t){
    #b    <- ffscb::make_band_FFSCB_z(tau=tau, t0=t0, diag.cov=diag.cov, conf.level=(1-p), n_int=n_int)
    b    <- make_band_FFSCB_z(tau=tau, t0=t0, diag.cov=diag.cov, conf.level=(1-p), n_int=n_int)
    b_f  <- stats::splinefun(x = seq(0,1,len=length(tau)), y = b, method = "natural")
    sgn  <- sign(x_f(t) - x0_f(t))
    tmp  <- x_f(t) - sgn * b_f(t) 
    return((tmp - x0_f(t))^2)
  }
  for(i in 1:length(t_grid)){
    p_grid[i] <- stats::optimize(f=function(p){myfun(p=p,t=t_grid[i])}, interval = c(0,1))$minimum
  }
  p_grid[p_grid>1] <- 1
  p_grid[p_grid<0] <- 0
  pvalue <- stats::spline(x = t_grid, y = p_grid, xout = seq(0,1,len=length(tau)), method = "natural")$y
  return(pvalue)
}


# FFSCB p-value (Gaussian)
#
# @param x function argument
# @param x0 Functional parameter under the null hypothesis. Zero function is assumed if it's not given.
# @param tau Pointwise standard deviation of the standardized and differentiated sample functions. Can be estimated by tau_fun().
# @param t0 Parameter t0 of the fast and fair simultaneous confidence bands.
# @param diag.cov The diagonal of N * Cov(X), in which X is the functional estimator. 
# @param N It should be '1' if 'cov' is the covariance operator for X itself, which is the default value.
# @param n.eval.points Number of evaluation points for the p-value function. Large values (>10) lead to slow computations.
# @param n_int Number of intervals parameter used by the function make_band_FFSCB_z()
# @references Liebl, D. and Reimherr, M. (2019). Fast and fair simultaneous confidence bands.
# @examples 
# # Generate a sample
# p <- 200 ; N <- 80 ; rangeval = c(0,1)
# grid  <- make_grid(p, rangevals=rangeval)
# mu0   <- meanf_poly(grid,c(0,1)) ; names(mu0) = grid
# mu    <- meanf_poly(grid,c(0,1.1)) ; names(mu) = grid
# cov.m <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2/2,1,1))
# dat   <- make_sample(mu,cov.m,N)
#
# # Find the estimate and covariance
# hat.mu       <- rowMeans(dat)
# hat.cov.m    <- crossprod(t(dat - hat.mu)) / (N-1)
# hat.tau.v    <- tau_fun(dat)
# 
# # pvalue
# pval <- get_pvalue_FFSCB_z(x=hat.mu, x0=mu0, tau=hat.tau.v, diag.cov=diag(hat.cov.m), N=N, 
# n.eval.points=6)
# plot(y=pval, x=grid, type="l", main="pvalue FFSCB (Gaussian)")
# @export
# get_pvalue_FFSCB_z <- function(x, x0=NULL, tau, t0=NULL, diag.cov, N, n.eval.points=11, n_int=5){
#   if (is.null(x0)) {x0 <- rep(0,times=length(x))}
#   t_grid        <- seq(0,1,len=n.eval.points)
#   p_grid        <- numeric(n.eval.points)
#   x_f           <- stats::splinefun(x = seq(0,1,len=length(tau)), y = x,  method = "natural")
#   x0_f          <- stats::splinefun(x = seq(0,1,len=length(tau)), y = x0, method = "natural")
#   diag.cov_f    <- stats::splinefun(x = seq(0,1,len=length(tau)), y = diag.cov, method = "natural")
#   ##
#   myfun <- function(p, t){
#     b    <- ffscb::make_band_FFSCB_z(x=x, tau=tau, t0=t0, diag.cov=diag.cov, N=N, conf.level=(1-p), n_int=n_int)
#     b    <- b[,2] - b[,1]
#     b_f  <- stats::splinefun(x = seq(0,1,len=length(tau)), y = b, method = "natural")
#     s    <- sign(x_f(t) - x0_f(t))
#     tmp  <- x_f(t) - s * b_f(t) 
#     ##
#     return((tmp - x0_f(t))^2)
#   }
#   for(i in 1:length(t_grid)){
#     p_grid[i] <- stats::optimize(f=function(p){myfun(p=p,t=t_grid[i])}, interval = c(0,1))$minimum
#   }
#   p_grid[p_grid>1] <- 1
#   p_grid[p_grid<0] <- 0
#   pvalue <- stats::spline(x = t_grid, y = p_grid, xout = seq(0,1,len=length(tau)), method = "natural")$y
#   return(pvalue)
# }




#' FFSCB p-value (t-distr)
#'
#' @param x Functional parameter estimate (for instance, the empirical mean function).
#' @param x0 Functional parameter under the null hypothesis. Default: zero.
#' @param tau Pointwise standard deviation of the standardized and differentiated sample functions. Can be estimated by tau_fun().
#' @param t0 Parameter t0 of the fast and fair simultaneous confidence bands.
#' @param diag.cov The diagonal of Cov(x), in which x is the functional estimator. For instance, the diagonal of the discretized covariance function of the empirical mean function x.  
#' @param df Degrees of freedom parameter for the t-distribution based band 'FFSCB.t'. (Typically, df=N-1)
#' @param n.eval.points Number of evaluation points for the p-value function. Large values (>10) lead to slow computations.
#' @param n_int Number of intervals parameter used by the function make_band_FFSCB_t()
#' @references Liebl, D. and Reimherr, M. (2019). Fast and fair simultaneous confidence bands.
#' @examples 
#' # Generate a sample
#' p <- 200 
#' N <- 80 
#' grid   <- make_grid(p, rangevals=c(0,1))
#' mu0    <- meanf_poly(grid,c(0,1))   ; names(mu0) <- grid
#' mu     <- meanf_poly(grid,c(0,1.1)) ; names(mu)  <- grid
#' cov.m  <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2/2,1,1))
#' sample <- make_sample(mu,cov.m,N)
#'
#' # Compute the estimate and its covariance
#' hat.mu     <- rowMeans(sample)
#' hat.cov    <- crossprod(t(sample - hat.mu)) / N
#' hat.cov.mu <- hat.cov / N
#' hat.tau    <- tau_fun(sample)
#' 
#' # Compute simultaneous pvalue function
#' pval <- get_pvalue_FFSCB_t(x=hat.mu, x0=mu0, tau=hat.tau, 
#'                            diag.cov=diag(hat.cov.mu), df=N-1, 
#'                            n.eval.points=6)
#' plot(y=pval, x=grid, type="l", main="pvalue FFSCB (t-distr)")
#' @export
get_pvalue_FFSCB_t <- function(x, x0=rep(0,times=length(x)), tau, t0=NULL, diag.cov, df, n.eval.points=11, n_int=5){
  t_grid        <- seq(0,1,len=n.eval.points)
  p_grid        <- numeric(n.eval.points)
  x_f           <- stats::splinefun(x = seq(0,1,len=length(tau)), y = x,  method = "natural")
  x0_f          <- stats::splinefun(x = seq(0,1,len=length(tau)), y = x0, method = "natural")
  diag.cov_f    <- stats::splinefun(x = seq(0,1,len=length(tau)), y = diag.cov, method = "natural")
  ##
  myfun <- function(p, t){
    #b    <- ffscb::make_band_FFSCB_t(tau=tau, t0=t0, diag.cov=diag.cov, df=df, conf.level=(1-p), n_int=n_int)
    b    <- make_band_FFSCB_t(tau=tau, t0=t0, diag.cov=diag.cov, df=df, conf.level=(1-p), n_int=n_int)
    b_f  <- stats::splinefun(x = seq(0,1,len=length(tau)), y = b, method = "natural")
    sgn  <- sign(x_f(t) - x0_f(t))
    tmp  <- x_f(t) - sgn * b_f(t) 
    return((tmp - x0_f(t))^2)
  }
  for(i in 1:length(t_grid)){
    p_grid[i] <- stats::optimize(f=function(p){myfun(p=p,t=t_grid[i])}, interval = c(0,1))$minimum
  }
  p_grid[p_grid>1] <- 1
  p_grid[p_grid<0] <- 0
  pvalue <- stats::spline(x = t_grid, y = p_grid, xout = seq(0,1,len=length(tau)), method = "natural")$y
  return(pvalue)
}



# FFSCB p-value (t-distr)
#
# @param x function argument
# @param x0 Functional parameter under the null hypothesis. Zero function is assumed if it's not given.
# @param tau Pointwise standard deviation of the standardized and differentiated sample functions. Can be estimated by tau_fun().
# @param t0 Parameter t0 of the fast and fair simultaneous confidence bands.
# @param diag.cov The diagonal of N * Cov(X), in which X is the functional estimator. 
# @param N It should be '1' if 'cov' is the covariance operator for X itself, which is the default value.
# @param n.eval.points Number of evaluation points for the p-value function. Large values (>10) lead to slow computations.
# @param n_int Number of intervals parameter used by the function make_band_FFSCB_z()
# @references Liebl, D. and Reimherr, M. (2019). Fast and fair simultaneous confidence bands.
# @example 
# # Generate a sample
# p <- 200 ; N <- 80 ; rangeval = c(0,1)
# grid  <- make_grid(p, rangevals=rangeval)
# mu0   <- meanf_poly(grid,c(0,1)) ; names(mu0) = grid
# mu    <- meanf_poly(grid,c(0,1.1)) ; names(mu) = grid
# cov.m <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2/2,1,1))
# dat   <- make_sample(mu,cov.m,N)
#
# # Find the estimate and covariance
# hat.mu       <- rowMeans(dat)
# hat.cov.m    <- crossprod(t(dat - hat.mu)) / (N-1)
# hat.tau.v    <- tau_fun(dat)
# 
# # pvalue
# pval <- get_pvalue_FFSCB_t(x=hat.mu, x0=mu0, tau=hat.tau.v, diag.cov=diag(hat.cov.m), N=N, 
# n.eval.points=6)
# plot(y=pval, x=grid, type="l", main="pvalue FFSCB (t-distr)")
# @export
# get_pvalue_FFSCB_t <- function(x, x0=NULL, tau, t0=NULL, diag.cov, N, n.eval.points=11, n_int=5){
#   if (is.null(x0)) {x0 <- rep(0,times=length(x))}
#   t_grid        <- seq(0,1,len=n.eval.points)
#   p_grid        <- numeric(n.eval.points)
#   x_f           <- stats::splinefun(x = seq(0,1,len=length(tau)), y = x,  method = "natural")
#   x0_f          <- stats::splinefun(x = seq(0,1,len=length(tau)), y = x0, method = "natural")
#   diag.cov_f    <- stats::splinefun(x = seq(0,1,len=length(tau)), y = diag.cov, method = "natural")
#   ##
#   myfun <- function(p, t){
#     b    <- ffscb::make_band_FFSCB_t(x=x, tau=tau, t0=t0, diag.cov=diag.cov, N=N, conf.level=(1-p), n_int=n_int)
#     b    <- b[,2] - b[,1]
#     b_f  <- stats::splinefun(x = seq(0,1,len=length(tau)), y = b, method = "natural")
#     s    <- sign(x_f(t) - x0_f(t))
#     tmp  <- x_f(t) - s * b_f(t) 
#     ##
#     return((tmp - x0_f(t))^2)
#   }
#   for(i in 1:length(t_grid)){
#     p_grid[i] <- stats::optimize(f=function(p){myfun(p=p,t=t_grid[i])}, interval = c(0,1))$minimum
#   }
#   p_grid[p_grid>1] <- 1
#   p_grid[p_grid<0] <- 0
#   pvalue <- stats::spline(x = t_grid, y = p_grid, xout = seq(0,1,len=length(tau)), method = "natural")$y
#   return(pvalue)
# }

