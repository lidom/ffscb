#' @export
get.pval.Enorm <- function(x, N = 1, eigen, fpc.cut=NULL, prec=NULL){ ## prec : precision option for CompQuadForm::imhof
  if (is.null(prec)) {prec <- c(10^(-6),10^(-6),10000)}
  if (is.null(fpc.cut) | fpc.cut==Inf){
    if (inherits(x,"fd")) {Enorm <- N * inprod(x,x) } else  ## fd object
                          {Enorm <- N * crossprod(x) }            ## vector
    validpc <- sum(eigen$values>0)
  } else {
    if (inherits(x,"fd")) {Enorm <- N * sum( inprod(x,eigen$harmonics[1:fpc.cut])^2) } else  ## fd object
                          {Enorm <- N * sum( crossprod(x,eigen$vectors[,1:fpc.cut])^2 ) }            ## vector
    validpc <- fpc.cut
  }
  CompQuadForm::imhof(Enorm,eigen$values[1:validpc],epsabs=prec[1],epsrel=prec[2],limit=prec[3])[[1]] ;
}
#' @export
get.pval.Epc <- function(x, N = 1, eigen, fpc.cut=NULL){
  if (is.null(fpc.cut) | fpc.cut==Inf) {fpc.cut <- sum(eigen$values > .Machine$double.eps)}
  if (inherits(x,"fd")) {Epc <- N * sum( inprod(x, eigen$harmonics[1:fpc.cut])^2 / eigen$values[1:fpc.cut] )} else
                        {Epc <- N * sum( crossprod( x, eigen$vectors[,1:fpc.cut])^2 / eigen$values[1:fpc.cut] )}
  1-pchisq(Epc, fpc.cut)
}
#' @export
get.pval.Ec <- function(x, N = 1, eigen, fpc.cut=NULL, prec=NULL){
  if (is.null(prec)) {prec <- c(10^(-6),10^(-6),10000)}
  if (is.null(fpc.cut) | fpc.cut==Inf) {fpc.cut <- sum(eigen$values > .Machine$double.eps)}
  c.square <- sqrt(eigen$values[1:fpc.cut])
  if (inherits(x,"fd")) {coef.square <- inprod(x, eigen$harmonics[1:fpc.cut])^2} else
                        {coef.square <- as.vector(crossprod(x, eigen$vectors[,1:fpc.cut])^2)}
  weights <- eigen$values[1:fpc.cut]/c.square/N
  stat <- sum(coef.square / c.square)
  CompQuadForm::imhof(stat, weights,epsabs=prec[1],epsrel=prec[2],limit=prec[3])[[1]]
}
#' @export
get.pval.Ec1 <- function(x, N = 1, eigen, fpc.cut=NULL, prec=NULL){
  if (is.null(prec)) {prec <- c(10^(-6),10^(-6),10000)}
  if (is.null(fpc.cut) | fpc.cut==Inf) {fpc.cut <- sum(eigen$values > .Machine$double.eps)}
  c.square <- sqrt(rev(cumsum(rev(eigen$values[1:fpc.cut])))) ## The only difference from Ec
  if (inherits(x,"fd")) {coef.square <- inprod(x, eigen$harmonics[1:fpc.cut])^2} else
                        {coef.square <- as.vector(crossprod(x, eigen$vectors[,1:fpc.cut])^2)}
  weights <- eigen$values[1:fpc.cut]/c.square/N
  stat <- sum(coef.square / c.square)
  CompQuadForm::imhof(stat, weights,epsabs=prec[1],epsrel=prec[2],limit=prec[3])[[1]]
}
#' @export
get.pval.Rz <- function(x, N = 1, eigen, fpc.cut=NULL){
  if (is.null(fpc.cut) | fpc.cut==Inf) {fpc.cut <- sum(eigen$values > .Machine$double.eps)}
  if (inherits(x,"fd")) {Zs.actual <- abs( sqrt(N) * inprod(x, eigen$harmonics[1:fpc.cut])   / sqrt(eigen$values[1:fpc.cut])) } else
                        {Zs.actual <- abs( sqrt(N) * crossprod(x, eigen$vectors[,1:fpc.cut]) / sqrt(eigen$values[1:fpc.cut])) }
  min(1-exp(sum(eigen$values[1:fpc.cut]) / eigen$values[1:fpc.cut] * log(pnorm2(Zs.actual))))
}
#' @export
get.pval.Rz1 <- function(x, N = 1, eigen, fpc.cut=NULL){
  if (is.null(fpc.cut) | fpc.cut==Inf) {fpc.cut <- sum(eigen$values > .Machine$double.eps)}
  if (inherits(x,"fd")) {Zs.actual <- abs( sqrt(N) * inprod(x, eigen$harmonics[1:fpc.cut])   / sqrt(eigen$values[1:fpc.cut])) } else
                        {Zs.actual <- abs( sqrt(N) * crossprod(x, eigen$vectors[,1:fpc.cut]) / sqrt(eigen$values[1:fpc.cut])) }
  M.star <- max(sqrt(2*pi) * eigen$values[1:fpc.cut] * ffscb::sf.f1(Zs.actual))
  Zs.star <- sapply(M.star / (sqrt(2*pi) * eigen$values[1:fpc.cut]), ffscb::sf.f1.inv)
  1 - prod(pnorm2(Zs.star))
}
#' @export
get.pval.Rzs <- function(x, N = 1, eigen, fpc.cut=NULL, hat.cov=NULL, df=NULL){
  if (is.null(hat.cov)) tilde.lambda <- eigen$values else {
    if (inherits(x,"fd")) { tilde.lambda <- diag(inprod(eigen$harmonics$basis, hat.cov$sbasis) %*% hat.cov$coefs %*%
                                                inprod(hat.cov$tbasis, eigen$harmonics$basis))} else
                          { tilde.lambda <- diag(t(eigen$vectors) %*% hat.cov %*% eigen$vectors) } # for vector x
  }
  if (is.null(fpc.cut) | fpc.cut==Inf) fpc.cut <- sum(tilde.lambda > .Machine$double.eps)
  if (inherits(x,"fd")) {Ts.actual <- abs( sqrt(N) * inprod(x, eigen$harmonics[1:fpc.cut] )   / sqrt(tilde.lambda[1:fpc.cut]))} else
                        {Ts.actual <- abs( sqrt(N) * crossprod(x, eigen$vectors[,1:fpc.cut] ) / sqrt(tilde.lambda[1:fpc.cut]))}
  if (is.null(df)) {df=N-1}
  # min(1-exp(sum(eigen$values[1:fpc.cut]) / eigen$values[1:fpc.cut] * log(pt2(Ts.actual,df=df))))
  min(1-exp(sum(tilde.lambda[1:fpc.cut]) / tilde.lambda[1:fpc.cut] * log(pt2(Ts.actual,df=df))))
}
#' @export
get.pval.Rz1s <- function(x, N = 1, eigen, fpc.cut=NULL, hat.cov=NULL, df=NULL){
  if (is.null(hat.cov)) tilde.lambda <- eigen$values else {
    if (inherits(x,"fd")) { tilde.lambda <- diag(inprod(eigen$harmonics$basis, hat.cov$sbasis) %*% hat.cov$coefs %*%
                                                  inprod(hat.cov$tbasis, eigen$harmonics$basis))} else
                          { tilde.lambda <- diag(t(eigen$vectors) %*% hat.cov %*% eigen$vectors) } # for vector x
  }
  if (is.null(fpc.cut) | fpc.cut==Inf) fpc.cut <- sum(tilde.lambda > .Machine$double.eps)
  if (inherits(x,"fd")) {Ts.actual <- abs( sqrt(N) * inprod(x, eigen$harmonics[1:fpc.cut] ) / sqrt(tilde.lambda[1:fpc.cut]))} else
                        {Ts.actual <- abs( sqrt(N) * crossprod(x, eigen$vectors[,1:fpc.cut] ) / sqrt(tilde.lambda[1:fpc.cut]))}
  if (is.null(df)) {df=N-1}
  Zs.actual <- qnorm2(pt2(Ts.actual,df=df))
  # M.star <- max(sqrt(2*pi) * eigen$values[1:fpc.cut] * ffscb::sf.f1(Zs.actual))
  # Zs.star <- sapply(M.star / (sqrt(2*pi) * eigen$values[1:fpc.cut]), ffscb::sf.f1.inv)
  M.star <- max(sqrt(2*pi) * tilde.lambda[1:fpc.cut] * ffscb::sf.f1(Zs.actual))
  Zs.star <- sapply(M.star / (sqrt(2*pi) * tilde.lambda[1:fpc.cut]), ffscb::sf.f1.inv)
  1 - prod(pnorm2(Zs.star))
}


#' @export
get.pvalue.FFSCB.z <- function(x, x0=rep(0,times=length(x)), tau, t0=NULL, diag.cov, N, n_int=10){
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
    p_grid[i] <- optimize(f=function(p){myfun(p=p,t=t_grid[i])}, interval = c(0,1))$minimum
  }
  p_grid[p_grid>1] <- 1
  p_grid[p_grid<0] <- 0
  pvalue <- stats::spline(x = t_grid, y = p_grid, xout = seq(0,1,len=length(tau)), method = "natural")$y
  return(pvalue)
}


#' @export
get.pvalue.FFSCB.t <- function(x, x0=rep(0,times=length(x)), tau, t0=NULL, diag.cov, N, n_int=10){
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
    p_grid[i] <- optimize(f=function(p){myfun(p=p,t=t_grid[i])}, interval = c(0,1))$minimum
  }
  p_grid[p_grid>1] <- 1
  p_grid[p_grid<0] <- 0
  pvalue <- stats::spline(x = t_grid, y = p_grid, xout = seq(0,1,len=length(tau)), method = "natural")$y
  return(pvalue)
}


# set.seed(1110)
# N             <- 10
# mu            <- meanf.poly(grid, params = c(0,0)) # plot(x=grid,y=mu)
# cov.m         <- make.cov.m(cov.f = covf.st.matern.warp.power_inv, grid=grid, cov.f.params=c(1.25, 1, 1, 2))
# dat           <- make.sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
# dat           <- dat + 0.5
# hat.mu        <- rowMeans(dat)
# hat.sd        <- mean(apply(dat, 1,sd))
# hat.cov.m     <- crossprod(t(dat - hat.mu)) / (N-1)
# hat.tau.v     <- tau_fun(dat)
# pvalue.FFSCB.z <- get.pvalue.FFSCB.z(x=hat.mu,tau=hat.tau.v, diag.cov=diag(hat.cov.m), N=N)
# pvalue.FFSCB.t <- get.pvalue.FFSCB.t(x=hat.mu,tau=hat.tau.v, diag.cov=diag(hat.cov.m), N=N)
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

