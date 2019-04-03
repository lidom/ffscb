#' @export
get_pval_Ec <- function(x, N = 1, eigen, fpc.cut=NULL, prec=NULL){
  if (is.null(prec)) {prec <- c(10^(-6),10^(-6),10000)}
  if (is.null(fpc.cut) | fpc.cut==Inf) {fpc.cut <- sum(eigen$values > .Machine$double.eps)}
  c.square <- sqrt(eigen$values[1:fpc.cut])
  if (inherits(x,"fd")) {coef.square <- fda::inprod(x, eigen$harmonics[1:fpc.cut])^2} else
                        {coef.square <- as.vector(crossprod(x, eigen$vectors[,1:fpc.cut])^2)}
  weights <- eigen$values[1:fpc.cut]/c.square/N
  stat    <- sum(coef.square / c.square)
  CompQuadForm::imhof(stat, weights,epsabs=prec[1],epsrel=prec[2],limit=prec[3])[[1]]
}

#' @export
get_pvalue_FFSCB_z <- function(x, x0=rep(0,times=length(x)), tau, t0=NULL, diag.cov, N, n_int=10){
  n.eval.points <- n_int + 1
  t_grid        <- seq(0,1,len=n.eval.points)
  p_grid        <- numeric(n.eval.points)
  x_f           <- stats::splinefun(x = seq(0,1,len=length(tau)), y = x,  method = "natural")
  x0_f          <- stats::splinefun(x = seq(0,1,len=length(tau)), y = x0, method = "natural")
  diag.cov_f    <- stats::splinefun(x = seq(0,1,len=length(tau)), y = diag.cov, method = "natural")
  ##
  myfun <- function(p, t){
    b    <- make.band.FFSCB.z(tau=tau, t0=t0, conf.level=(1-p), n_int=5)
    b_f  <- stats::splinefun(x = seq(0,1,len=length(tau)), y = b, method = "natural")
    s    <- sign(x_f(t) - x0_f(t))
    tmp  <- x_f(t) - s * b_f(t) * sqrt(diag.cov_f(t)) / sqrt(N)
    ##
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


#' @export
get_pvalue_FFSCB_t <- function(x, x0=rep(0,times=length(x)), tau, t0=NULL, diag.cov, N, n_int=10){
  n.eval.points <- n_int + 1
  t_grid        <- seq(0,1,len=n.eval.points)
  p_grid        <- numeric(n.eval.points)
  x_f           <- stats::splinefun(x = seq(0,1,len=length(tau)), y = x,  method = "natural")
  x0_f          <- stats::splinefun(x = seq(0,1,len=length(tau)), y = x0, method = "natural")
  diag.cov_f    <- stats::splinefun(x = seq(0,1,len=length(tau)), y = diag.cov, method = "natural")
  ##
  myfun <- function(p, t){
    b    <- ffscb::make.band.FFSCB.t(tau=tau, t0 = t0, conf.level=(1-p), N=N, n_int=5)
    b_f  <- stats::splinefun(x = seq(0,1,len=length(tau)), y = b, method = "natural")
    s    <- sign(x_f(t) - x0_f(t))
    tmp  <- x_f(t) - s * b_f(t) * sqrt(diag.cov_f(t)) / sqrt(N)
    ##
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


# set.seed(1110)
# N             <- 10
# mu            <- meanf_poly(grid, params = c(0,0)) # plot(x=grid,y=mu)
# cov.m         <- make_cov_m(cov.f = covf.st.matern.warp.power_inv, grid=grid, cov.f.params=c(1.25, 1, 1, 2))
# dat           <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
# dat           <- dat + 0.5
# hat.mu        <- rowMeans(dat)
# hat.sd        <- mean(apply(dat, 1,sd))
# hat.cov.m     <- crossprod(t(dat - hat.mu)) / (N-1)
# hat.tau.v     <- tau_fun(dat)
# pvalue.FFSCB.z <- get_pvalue_FFSCB_z(x=hat.mu,tau=hat.tau.v, diag.cov=diag(hat.cov.m), N=N)
# pvalue.FFSCB.t <- get_pvalue_FFSCB_t(x=hat.mu,tau=hat.tau.v, diag.cov=diag(hat.cov.m), N=N)
# pvalue.t.test <- unname(apply(dat, 1, function(x)t.test(x)$p.value))
# 
# par(mfrow=c(2,1))
# matplot(y=cbind(pvalue.FFSCB.z, pvalue.FFSCB.t, pvalue.t.test), 
#         x=seq(0,1,len=nrow(dat)), type="l", pch = c("z","t",NA), 
#         xlab="", ylab="", lty=c(1:3), col=1, log="y", main="p-value")
# legend("bottomright", legend = c("FFSCB.z", "FFSCB.t", "pointwise t"), lty=c(1:3))
# plot(y=hat.mu, x=seq(0,1,len=nrow(dat)), type="l", main="Mean Estimation",
#      ylim = range(hat.sd/2,-hat.sd/2,0), xlab = "", ylab="")
# abline(h=0,lty=2)
# par(mfrow=c(1,1))

