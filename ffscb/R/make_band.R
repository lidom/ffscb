make_band_Ec <- function(eigen, conf.level, fd.eval.grid.size=200){
  alpha.level <- 1-conf.level
  pc.to.use   <- sum(eigen$values > .Machine$double.eps)
  c.square    <- sqrt(eigen$values[1:pc.to.use])
  weights     <- eigen$values[1:pc.to.use]/c.square
  xi          <- get.schisq.q.gamma(weights,conf.level) ## Approximate Quantile of Weighted Sum of Chi-square by Gamma
  if (inherits(eigen,"pca.fd") | inherits(eigen,"eigen.fd")) {
    evalgrid      <- ffscb::make_grid(p=fd.eval.grid.size, rangevals=eigen$harmonics$basis$rangeval)
    eigen$vectors <- fda::eval.fd(evalgrid,eigen$harmonics)
  }
  band.eval <- sqrt(apply(t(eigen$vectors[,1:pc.to.use]^2) * c.square * xi, 2, sum))
  if (inherits(eigen,"pca.fd") | inherits(eigen,"eigen.fd")) {
    return(fda::Data2fd(evalgrid,band.eval,basisobj=eigen$harmonics$basis)) # return as fd object
  } else return(band.eval)                                                  # return as vector
}


make_band_Bs <- function(cov, conf.level, sim.size=10000, fd.eval.grid.size=200){
  if (inherits(cov,"bifd")) {
    evalgrid <- ffscb::make_grid(p=fd.eval.grid.size, rangevals=cov$sbasis$rangeval)
    cov.m    <- fda::eval.bifd(evalgrid,evalgrid,cov) 
  } else {
    cov.m <- cov
  }
  crit.Bs    <- get.crit.supnorm.simple(cov.m = cov.m, n.sim = sim.size, prob = conf.level)
  band.eval  <- sqrt(diag(cov.m)) * crit.Bs
  if (inherits(cov,"bifd")) {
    return(fda::Data2fd(evalgrid,band.eval,basisobj=cov$sbasis))
  } else return(band.eval)
}


make_band_naive_t <- function(cov, conf.level, df, fd.eval.grid.size=200){
  if (inherits(cov,"bifd")) {
    evalgrid <- ffscb::make_grid(p=fd.eval.grid.size, rangevals=cov$sbasis$rangeval)
    cov.m    <- fda::eval.bifd(evalgrid,evalgrid,cov) 
  } else {
    cov.m <- cov
  }
  band.eval  <- qt2(conf.level,df) * sqrt(diag(cov.m))
  if (inherits(cov,"bifd")) {
    return(fda::Data2fd(evalgrid,band.eval,basisobj=cov$sbasis))
  } else return(band.eval)
}

make_band_naive_t_fragm <- function(diag.cov, conf.level, df){
  ##
  band.eval  <- qt2(conf.level,df) * sqrt(diag.cov)
  ##
  return(band.eval)
}


# Kac-Rice simultaneous confidence band (Gaussian)
# 
# @param x Functional parameter estimate.
# @param tau Pointwise standard deviation of the standardized and differentiated sample functions. Can be estimated by tau_fun().
# @param diag.cov The diagonal of N * Cov(X), in which X is the functional estimator.
# @param N It should be '1' if 'cov' is the covariance operator for X itself, which is the default value.
# @param conf.level confidence level (default: 0.95)
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
# # Result
# band <- make_band_KR_z(x=hat.mu, tau=hat.tau.v, diag.cov=diag(hat.cov.m),N=N)
# matplot(y=band, x=grid, type="l", lty=c(2,1,2), col=1, main="KR-band (Gaussian)")
# @export
# make_band_KR_z <- function(x, tau, diag.cov, N, conf.level=0.95){
#   alpha.level <- 1-conf.level
#   tt          <- seq(0,1,len=length(tau))
#   tau_01      <- sum(tau)*diff(tt)[1] # int_0^1 tau(t) dt
#   myfun       <- function(c){stats::pnorm(c,lower.tail = FALSE)+tau_01*exp(-c^2/2)/(2*pi)-alpha.level/2}
#   cstar       <- stats::uniroot(f = myfun,interval = c(.5,8))$root
#   band        <- rep(cstar, times=length(tau)) * sqrt(diag.cov) / sqrt(N)
#   band_upper  <- x + band
#   band_lower  <- x - band
#   ##
#   result           <- cbind(band_lower,x,band_upper)
#   colnames(result) <- c("KR_z_band_lower","x","KR_z_band_upper")
#   ##
#   return(result)
# }

make_band_KR_z <- function(tau, diag.cov, conf.level=0.95){
  alpha.level <- 1-conf.level
  tt          <- seq(0,1,len=length(tau))
  tau_01      <- sum(tau)*diff(tt)[1] # int_0^1 tau(t) dt
  myfun       <- function(c){stats::pnorm(c,lower.tail = FALSE) + tau_01*exp(-c^2/2)/(2*pi)-alpha.level/2}
  cstar       <- stats::uniroot(f = myfun,interval = c(.5,8))$root
  band        <- rep(cstar, times=length(tau)) * sqrt(diag.cov) 
  ##
  return(band)
}


# Kac-Rice simultaneous confidence band (t-distr)
#
# @param x Functional parameter estimate. 
# @param tau Pointwise standard deviation of the standardized and differentiated sample functions. Can be estimated by tau_fun().
# @param diag.cov The diagonal of N * Cov(X), in which X is the functional estimator. 
# @param N It should be '1' if 'cov' is the covariance operator for X itself, which is the default value.
# @param conf.level confidence level (default: 0.95)
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
# # Result
# band <- make_band_KR_t(x=hat.mu, tau=hat.tau.v, diag.cov=diag(hat.cov.m),N=N)
# matplot(y=band, x=grid, type="l", lty=c(2,1,2), col=1, main="KR-band (t-distr)")
# @export
# make_band_KR_t <- function(x, tau, diag.cov, N, conf.level=0.95){
#   alpha.level <- 1-conf.level
#   nu          <- N-1
#   tt          <- seq(0,1,len=length(tau))
#   tau_01      <- sum(tau)*diff(tt)[1] # int_0^1 tau(t) dt
#   myfun       <- function(c){stats::pt(c, lower.tail = FALSE, df=nu)+tau_01*(1+c^2/nu)^(-nu/2)/(2*pi) - alpha.level/2}
#   cstar       <- stats::uniroot(f = myfun,interval = c(.5,8))$root
#   band        <- rep(cstar, times=length(tau)) * sqrt(diag.cov) / sqrt(N)
#   band_upper  <- x + band
#   band_lower  <- x - band
#   ##
#   result           <- cbind(band_lower,x,band_upper)
#   colnames(result) <- c("KR_z_band_lower","x","KR_z_band_upper")
#   ##
#   return(result)
# }


make_band_KR_t <- function(tau, diag.cov, df, conf.level=0.95){
  alpha.level <- 1-conf.level
  tt          <- seq(0,1,len=length(tau))
  tau_01      <- sum(tau)*diff(tt)[1] # int_0^1 tau(t) dt
  myfun       <- function(c){stats::pt(c, lower.tail = FALSE, df=df) + tau_01*(1+c^2/df)^(-df/2)/(2*pi) - alpha.level/2}
  cstar       <- stats::uniroot(f = myfun,interval = c(.5,8))$root
  band        <- rep(cstar, times=length(tau)) * sqrt(diag.cov)
  ##
  return(band)
}


#' Fast 'n' fair simultaneous confidence band (Gaussian)
#'
#' @param x Functional parameter estimate. 
#' @param diag.cov.x Diagonal of Cov(x), in which x is the functional estimator (for instance, the covariance function of the empirical mean function).
#' @param tau Pointwise standard deviation of the standardized and differentiated sample functions. Can be estimated by tau_fun().
#' @param t0 Parameter t0 of the fast and fair simultaneous confidence bands. If left unspecified (default t0=NULL), t0 is set to the location which maximizes tau.
#' @param conf.level confidence level (default: 0.95)
#' @param n_int Number of equidistant intervals over which the multiple testing component of the type-I error rate (1-conf.level) is distributed uniformly.
#' @param tol tolerance 'tol' parameter used by stats::uniroot(). If tol=NULL, we use tol=.Machine$double.eps^0.32 which increases the default accuracy used by uniroot().
#' @references Liebl, D. and Reimherr, M. (2019). Fast and fair simultaneous confidence bands.
#' @examples 
#' # Generate a sample
#' p          <- 200 
#' N          <- 80 
#' grid       <- make_grid(p, rangevals=c(0,1))
#' mu         <- meanf_poly(grid,c(0,1.1)) 
#' names(mu)  <- grid
#' cov.m      <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2/2,1,1))
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
#' b <- make_band_FFSCB_z(x=hat.mu, diag.cov.x=diag(hat.cov.mu), tau=hat.tau,
#'                        conf.level  = 0.95)
#' matplot(y=b$band[,2:3], x=grid, lty=2)
#' lines(x=grid, y=b$band[,1], lty=1)
#' @export
make_band_FFSCB_z <- function(x, diag.cov.x, tau, t0=0, conf.level=0.95, n_int=4, tol=NULL){
  ##
  if(any(tau < 0.005)){warning("This method may not work if tau(t) is too small.")}
  ##
  result_tmp       <- .make_band_FFSCB_z(tau=tau, t0=t0, diag.cov=diag.cov.x, conf.level=conf.level, n_int=n_int, tol=tol)
  band_m           <- cbind(x, x + result_tmp$band, x - result_tmp$band)
  colnames(band_m) <- c("x", paste0("FFSCB.z.u.", conf.level), paste0("FFSCB.z.l.", conf.level))
  ##
  tt          <- seq(0,1,len=length(tau))
  if(is.null(t0)){t0 <- tt[which.max(tau)]}
  ##
  return(list("band"    = band_m,
              "t0"      = result_tmp$t0,
              "prob_t0" = result_tmp$prob_t0 * 2, # multipled by two, since two-sided SCB
              "a_star"  = result_tmp$a_star  * 2  # multipled by two, since two-sided SCB
              ))
}


.make_band_FFSCB_z <- function(tau, t0=0, diag.cov, conf.level=0.95, n_int=4, tol=NULL){
  ##
  alpha.level <- 1-conf.level
  tt          <- seq(0,1,len=length(tau))
  tau_v       <- tau
  tau_f       <- stats::approxfun(x = seq(0,1,len=length(tau)), y = tau)
  knots       <- seq(0,1,len=(n_int + 1))
  ##
  if(is.null(tol)){
    tol         <- .Machine$double.eps^0.32 # increases the default accuracy for uniroot() (.Machine$double.eps^0.25)
  }
  if(!is.numeric(t0)){
    warning("t0 is not given and set to the default `t0=0`.")
    t0 <- 0
  }
  if(!any(t0==knots) ){ 
    warning("t0 does not equal one of the grid points of the threshold u.\nWe set t0 to the closest grid point.")
    t0 <- knots[which.min(abs(t0-knots))]
  }
  ##
  if(n_int == 1){# Case n_int=1 == constant band == Kac-Rice Band
    tau01      <- sum(tau_v)*diff(tt)[1] # int_0^1 tau(t) dt
    myfun1     <- function(c1){stats::pnorm(q=c1, lower.tail=F)+exp(-c1^2/2)*tau01/(2*pi)-(alpha.level/2)}
    const_band <- stats::uniroot(f = myfun1, interval = c(0,10), extendInt="downX")$root
    const_band <- const_band * sqrt(diag.cov) 
    return(const_band)
  }
  ## Define piecewise linear (pwl) function 'ufun' with derivative=0 at t0.
  c_v         <- numeric(n_int) # coeficients of the pwl-function
  const_int   <- min(findInterval(t0, knots, rightmost.closed = TRUE), n_int)
  fct_body    <- paste0("c_v[",const_int,"]")
  if(const_int >     1){for(j in (const_int-1):1    ){fct_body <- paste0("c_v[",j,"]*pmin(t - knots[",j+1,"],0) +", fct_body)}}
  if(const_int < n_int){for(j in (const_int+1):n_int){fct_body <- paste0(fct_body, "+ c_v[",j,"]*pmax(t - knots[",j,"],0)")}}
  ufun        <- function(t, c_v, knots){}
  body(ufun)  <- parse(text=fct_body)
  ##
  find_u <- function(alpha.aux){# alpha.aux = 0.006
    ##
    ## Determine the initial value of u
    tau_init <- sum(tau_v[knots[const_int] <= tt & tt <= knots[const_int+1]])*diff(tt)[1]
    if(tau_init > 2*pi*(alpha.aux/2)/n_int){# if possible use the analytic solution
      c_v[const_int] <- sqrt(2) * sqrt(log(tau_init/(2*pi*(alpha.aux/2)/n_int))) 
    }else{
      myfun1         <- function(c1){c(exp(-c1^2/2) * tau_init / (2*pi) - (alpha.aux/2)/n_int)}
      # curve(myfun1, 0,5); abline(h=0)
      c_v[const_int] <- stats::uniroot(f = myfun1, interval = c(0,10), extendInt = "downX", tol = tol)$root
    }
    ##
    if(const_int > 1){
      for(j in (const_int-1):1){
        myfunj <- function(cj){
          ##
          if(j==(const_int-1)){c_v_sum <- 0}else{c_v_sum <- c_v[(const_int-1):(j+1)]}# sum(c_v_sum) == u' of the preceeding interval
          ##
          ufun_j <- function(t,cj){ufun(t=t,c_v=c(rep(0,times=(j-1)),cj,c_v[(j+1):const_int],rep(0, times=(n_int-const_int))), knots=knots)}
          ##
          fn1    <- function(t,cj){(tau_f(t)/(2*pi)) * exp(-ufun_j(t,cj)^2/2) * exp(-sum(c(c_v_sum,cj))^2/(2*tau_f(t)^2))}
          fn2    <- function(t,cj){sum(c(c_v_sum,cj))/sqrt(2*pi) * exp(-ufun_j(t,cj)^2/2) * stats::pnorm( sum(c(c_v_sum, cj))/tau_f(t))}
          intgr1 <- sum(fn1(t=tt[knots[j] < tt & tt <= knots[j+1]], cj=cj)) * diff(tt)[1]
          intgr2 <- sum(fn2(t=tt[knots[j] < tt & tt <= knots[j+1]], cj=cj)) * diff(tt)[1]
          ##
          res    <- c(intgr1 + intgr2 - (alpha.aux/2)/n_int)
          return(res)
        }
        # myfunj <- Vectorize(myfunj); curve(myfunj, -10,10); abline(h=0)
        c_v[j] <- stats::uniroot(f = myfunj, interval = c(-10,10), extendInt = "upX", tol = tol)$root
      }
    }
    if(const_int < n_int){
      for(j in (const_int+1):n_int){# j <- (const_int+1)
        myfunj <- function(cj){
          ##
          if(j==(const_int+1)){c_v_sum <- 0}else{c_v_sum <- c_v[(const_int+1):(j-1)]}# sum(c_v_sum) == u'(t) for 0 <= t < (j-1)th interval
          ##
          ufun_j <- function(t,cj){ufun(t=t,c_v=c(rep(0,times=(const_int-1)),c_v[const_int:(j-1)],cj,rep(0, times=(n_int-j))), knots=knots)}
          ##
          fn1    <- function(t,cj){(tau_f(t)/(2*pi)) * exp(-ufun_j(t,cj)^2/2) * exp(-sum(c(c_v_sum,cj))^2/(2*tau_f(t)^2))}
          fn3    <- function(t,cj){sum(c(c_v_sum,cj))/sqrt(2*pi) * exp(-ufun_j(t,cj)^2/2) * stats::pnorm(-sum(c(c_v_sum, cj))/tau_f(t))}
          intgr1 <- sum(fn1(t=tt[knots[j] < tt & tt <= knots[j+1]], cj=cj)) * diff(tt)[1]
          intgr3 <- sum(fn3(t=tt[knots[j] < tt & tt <= knots[j+1]], cj=cj)) * diff(tt)[1]
          ##
          res    <- c(intgr1 - intgr3 - (alpha.aux/2)/n_int)
          return(res)
        }
        # myfunj <- Vectorize(myfunj); curve(myfunj, -10,10); abline(h=0)
        c_v[j] <- stats::uniroot(f = myfunj, interval = c(-10,10), extendInt = "downX", tol = tol)$root
      }
    }
    ## 
    u_star_f  <- function(t){return(ufun(t=t,c_v=c_v,knots=knots))}
    up_star_f <- function(t){
      cp <- c_v; cp[const_int] <- 0
      ##
      if(const_int < n_int){ cp[(const_int+1):n_int] <- cumsum(cp[(const_int+1):n_int]) }
      if(const_int > 1    ){ cp[(const_int-1):1]     <- cumsum(cp[(const_int-1):1])     }
      ##
      return(cp[findInterval(t, knots, rightmost.closed = TRUE)])
    }
    # up_star_f <- function(t){
    #   kn   <- knots[-length(knots)]; cp <- c_v; cp[const_int] <- 0
    #   kn_m <- matrix(kn,                     nrow=length(kn), ncol=length(t))
    #   cp_m <- matrix(cp,                     nrow=length(kn), ncol=length(t))
    #   t_m  <- matrix(rep(t,each=length(kn)), nrow=length(kn), ncol=length(t))
    #   return(colSums(cp_m * (kn_m<t_m)))
    # }
    ##
    fn1_star <- function(t){(tau_f(t)/(2*pi)) * exp(-u_star_f(t)^2/2) * exp(-up_star_f(t)^2/(2*tau_f(t)^2))}
    fn2_star <- function(t){up_star_f(t)/sqrt(2*pi) * exp(-u_star_f(t)^2/2) * stats::pnorm( up_star_f(t)/tau_f(t))}
    fn3_star <- function(t){up_star_f(t)/sqrt(2*pi) * exp(-u_star_f(t)^2/2) * stats::pnorm(-up_star_f(t)/tau_f(t))}
    ##
    intgr1_star <- sum(fn1_star(t=tt)) * diff(tt)[1]
    ##
    if(t0>0){intgr2_star <- sum(fn2_star(t=tt[tt <= t0])) * diff(tt)[1]}else{intgr2_star <- 0}
    if(t0<1){intgr3_star <- sum(fn3_star(t=tt[t0 <  tt])) * diff(tt)[1]}else{intgr3_star <- 0}
    ##
    optim_target <- c( (stats::pnorm(q=u_star_f(t0), lower.tail=FALSE) + intgr1_star + intgr2_star - intgr3_star) - alpha.level/2)
    band.eval    <- u_star_f(t=tt)
    ##
    return(list("optim_target"     = optim_target,
                "band.eval"        = band.eval, 
                "prob_t0"          = stats::pnorm(q=u_star_f(t0), lower.tail=FALSE),
                "a_star"           = (intgr1_star + intgr2_star - intgr3_star) 
                ))
  }
  ##
  fu <- function(x){find_u(x)$optim_target}
  # fu <- Vectorize(fu); curve(fu, .Machine$double.eps, 1); abline(h=0)
  ##
  opt_res <- stats::uniroot(f = fu, interval = c(.Machine$double.eps, (alpha.level/2)), extendInt = "upX", tol=tol)$root
  band    <- find_u(opt_res)$band.eval * sqrt(diag.cov)
  prob_t0 <- find_u(opt_res)$prob_t0
  a_star  <- find_u(opt_res)$a_star
  ##
  if(abs(prob_t0 + a_star - alpha.level/2) > ( (alpha.level/2) / 25 ) ){
    stop("Numeric equation solvers 'uniroot()' give too imprecise results. \nPlease, set the 'tol' argument smaller than the default value (<.Machine$double.eps^0.32).")
  }
  ##
  return(list("band"    = band,
              "t0"      = t0,
              "prob_t0" = prob_t0,
              "a_star"  = a_star))
}

#' Fast 'n' fair simultaneous confidence band (t-distr)
#'
#' @param x Functional parameter estimate. 
#' @param diag.cov.x Diagonal of Cov(x), in which x is the functional estimator (for instance, the covariance function of the empirical mean function).
#' @param tau Pointwise standard deviation of the standardized and differentiated sample functions. Can be estimated by tau_fun().
#' @param t0 Parameter t0 of the fast and fair simultaneous confidence bands. If left unspecified (default t0=NULL), t0 is set to the location which maximizes tau.
#' @param df Degrees of freedom 
#' @param conf.level confidence level (default: 0.95)
#' @param n_int Number of equidistant intervals over which the multiple testing component of the type-I error rate (1-conf.level) is distributed uniformly.
#' @param tol tolerance 'tol' parameter used by stats::uniroot(). If tol=NULL, we use tol=.Machine$double.eps^0.32 which increases the default accuracy used by uniroot().
#' @references Liebl, D. and Reimherr, M. (2019). Fast and fair simultaneous confidence bands.
#' @examples 
#' # Generate a sample
#' p          <- 200 
#' N          <- 80 
#' grid       <- make_grid(p, rangevals=c(0,1))
#' mu         <- meanf_poly(grid,c(0,1.1)) 
#' names(mu)  <- grid
#' cov.m      <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2/2,1,1))
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
#' b <- make_band_FFSCB_t(x=hat.mu, diag.cov.x=diag(hat.cov.mu), tau=hat.tau,
#'                        df = N-1, conf.level  = 0.95)
#' matplot(y=b$band[,2:3], x=grid, lty=2)
#' lines(x=grid, y=b$band[,1], lty=1)
#' @export
make_band_FFSCB_t <- function(x, diag.cov.x, tau, t0=0, df, conf.level=0.95, n_int=4, tol=NULL){
  ##
  if(any(tau < 0.005)){warning("This method may not work if tau(t) is too small.")}
  ##
  if(df <= 100){
    result_tmp       <- .make_band_FFSCB_t(tau=tau, t0=t0, diag.cov=diag.cov.x, df=df, conf.level=conf.level, n_int=n_int, tol=tol)
  }else{
    result_tmp       <- .make_band_FFSCB_z(tau=tau, t0=t0, diag.cov=diag.cov.x,        conf.level=conf.level, n_int=n_int, tol=tol)
  }
  band_m           <- cbind(x, x + result_tmp$band, x - result_tmp$band)
  colnames(band_m) <- c("x", paste0("FFSCB.t.u.", conf.level), paste0("FFSCB.t.l.", conf.level))
  ##
  tt          <- seq(0,1,len=length(tau))
  if(is.null(t0)){t0 <- tt[which.max(tau)]}
  ##
  return(list("band"    = band_m,
              "t0"      = result_tmp$t0,
              "prob_t0" = result_tmp$prob_t0 *2, # multipled by two, since two-sided SCB
              "a_star"  = result_tmp$a_star  *2  # multipled by two, since two-sided SCB
              ))
}


.make_band_FFSCB_t <- function(tau, t0=0, diag.cov, df, conf.level=0.95, n_int=4, tol=NULL){
  ##
  alpha.level <- 1-conf.level
  tt          <- seq(0,1,len=length(tau))
  #tau_v       <- stats::spline(   x = seq(0,1,len=length(tau)), y = tau, method = "natural", xout = tt)$y
  #tau_f       <- stats::splinefun(x = seq(0,1,len=length(tau)), y = tau, method = "natural")
  tau_v       <- tau
  tau_f       <- stats::approxfun(x = seq(0,1,len=length(tau)), y = tau)
  knots       <- seq(0,1,len=(n_int + 1))
  if(!is.numeric(t0)){
    warning("t0 is not given and set to the default `t0=0`.")
    t0 <- 0
  }
  if(!any(t0==knots) ){ 
    warning("t0 does not equal one of the grid points of the threshold u.\nWe set t0 to the closest grid point.")
    t0 <- knots[which.min(abs(t0-knots))]
  }
  if(is.null(tol)){
    tol         <- .Machine$double.eps^0.32 # increases the default accuracy for uniroot() (.Machine$double.eps^0.25)
  }
  nu          <- df
  nup         <- nu+1
  #if(is.null(t0)){t0 <- tt[which.max(tau_v)]}else{t0 <- tt[findInterval(t0, tt, rightmost.closed = TRUE)]}
  ##
  if(n_int == 1){# Case n_int=1 == constant band == Kac-Rice Band
    tau01      <- sum(tau_v)*diff(tt)[1] # int_0^1 tau(t) dt
    myfun1     <- function(c1){c(stats::pt(q=c1, lower.tail=FALSE, df = nu)+(tau01/(2*pi))*(1+c1^2/nu)^(-nu/2)-(alpha.level/2))}
    const_band <- stats::uniroot(f = myfun1, interval = c(0,10), extendInt="downX")$root
    const_band <- const_band * sqrt(diag.cov) 
    return(const_band)
  }
  ## Define piecewise linear (pwl) function 'ufun' with derivative=0 at t0.
  c_v         <- numeric(n_int) # coeficients of the pwl-function
  const_int   <- min(findInterval(t0, knots, rightmost.closed = TRUE), n_int)
  fct_body    <- paste0("c_v[",const_int,"]")
  if(const_int >     1){for(j in (const_int-1):1    ){fct_body <- paste0("c_v[",j,"]*pmin(t - knots[",j+1,"],0) +", fct_body)}}
  if(const_int < n_int){for(j in (const_int+1):n_int){fct_body <- paste0(fct_body, "+ c_v[",j,"]*pmax(t - knots[",j,"],0)")}}
  ufun        <- function(t, c_v, knots){}
  body(ufun)  <- parse(text=fct_body)
  ##
  find_u <- function(alpha.aux){
    ##
    ## Determine the initial value of u
    tau_init <- sum(tau_v[knots[const_int] <= tt & tt <= knots[const_int+1]])*diff(tt)[1]
    if(nu*(( (pi*alpha.aux) / ( tau_init * n_int) )^(-2/nu) -1) >0){# if possible use the analytic solution
      c_v[const_int] <- sqrt(nu*(( (pi*alpha.aux) / ( tau_init *n_int) )^(-2/nu) -1))
    }else{ 
      myfun1         <- function(c1){c((tau_init/(2*pi))*(1+c1^2/nu)^(-nu/2)-(alpha.aux/2)/n_int)}
      # curve(myfun1, 0,1); abline(h=0)
      c_v[const_int] <- stats::uniroot(f = myfun1, interval = c(0,10), extendInt = "downX", tol = tol)$root
    }
    ##
    if(const_int > 1){
      for(j in (const_int-1):1){
        myfunj <- function(cj){
          ##
          if(j==(const_int-1)){c_v_sum <- 0}else{c_v_sum <- c_v[(const_int-1):(j+1)]}# sum(c_v_sum) == u' of the preceeding interval
          ##
          ufun_j <- function(t,cj){ufun(t=t,c_v=c(rep(0,times=(j-1)),cj,c_v[(j+1):const_int],rep(0, times=(n_int-const_int))), knots=knots)}
          afun_j <- function(t,cj){sqrt(nu*tau_f(t)^2*(1+ufun_j(t,cj)^2/nu)/nup)}
          ##
          fn1 <- function(t,cj){tau_f(t) * (1 + ufun_j(t,cj)^2/nu + sum(c(c_v_sum,cj))^2/(nu*tau_f(t)^2))^(-nu/2) / (2*pi)}
          fn2 <- function(t,cj){
            (sum(c(c_v_sum,cj))/(2*pi*tau_f(t))) * (1+ufun_j(t,cj)^2/nu)^(-nu/2 -1)  *
              (gamma(nup/2) * sqrt(nup*pi) * afun_j(t,cj) / gamma((nup+1)/2) ) *
              stats::pt(q = (sum(c(c_v_sum,cj)) / afun_j(t,cj) ), df=nup) 
          }
          res <- ( (sum(fn1(t=tt[knots[j] < tt & tt <= knots[j+1]],cj=cj)) * diff(tt)[1]  +
                      sum(fn2(t=tt[knots[j] < tt & tt <= knots[j+1]],cj=cj)) * diff(tt)[1]) - (alpha.aux/2)/n_int )
          return(res)
        } 
        # myfunj <- Vectorize(myfunj); curve(myfunj, -10,10); abline(h=0)
        c_v[j] <- stats::uniroot(f = myfunj, interval = c(-10,10), extendInt = "upX", tol = tol)$root
      }
    }
    if(const_int < n_int){
      for(j in (const_int+1):n_int){
        myfunj <- function(cj){
          ##
          if(j==(const_int+1)){c_v_sum <- 0}else{c_v_sum <- c_v[(const_int+1):(j-1)]}# sum(c_v_sum) == u'(t) for 0 <= t < (j-1)th interval
          ##
          #ufun_j <- function(t,cj){ufun(t=t,c_v=c(c_v[1:(j-1)],cj,rep(0, times=(n_int-j))), knots=knots)}
          ufun_j <- function(t,cj){ufun(t=t,c_v=c(rep(0,times=(const_int-1)),c_v[const_int:(j-1)],cj,rep(0, times=(n_int-j))), knots=knots)}
          afun_j <- function(t,cj){sqrt(nu*tau_f(t)^2*(1+ufun_j(t,cj)^2/nu)/nup)}
          ##
          fn1 <- function(t,cj){tau_f(t) * (1 + ufun_j(t,cj)^2/nu + sum(c(c_v_sum,cj))^2/(nu*tau_f(t)^2))^(-nu/2) / (2*pi)}
          fn3 <- function(t,cj){
            (sum(c(c_v_sum,cj))/(2*pi*tau_f(t))) * 
              (1+ufun_j(t,cj)^2/nu)^(-nu/2 -1)  *
              (gamma(nup/2) * sqrt(nup*pi) * afun_j(t,cj) / gamma((nup+1)/2) ) *
              stats::pt(q = (-sum(c(c_v_sum,cj)) / afun_j(t,cj) ), df=nup) 
          }
          res <- ( (sum(fn1(t=tt[knots[j] < tt & tt <= knots[j+1]],cj=cj)) * diff(tt)[1]  -
                      sum(fn3(t=tt[knots[j] < tt & tt <= knots[j+1]],cj=cj)) * diff(tt)[1]) - (alpha.aux/2)/n_int )
          return(res)
        } 
        # myfunj <- Vectorize(myfunj); curve(myfunj, -10,10); abline(h=0)
        c_v[j] <- stats::uniroot(f = myfunj, interval = c(-10,10), extendInt = "downX", tol = tol)$root
      }
    }
    ## 
    u_star_f  <- function(t){return(ufun(t=t,c_v=c_v,knots=knots))}
    up_star_f <- function(t){
      cp <- c_v; cp[const_int] <- 0
      ##
      if(const_int < n_int){ cp[(const_int+1):n_int] <- cumsum(cp[(const_int+1):n_int]) }
      if(const_int > 1    ){ cp[(const_int-1):1]     <- cumsum(cp[(const_int-1):1])     }
      ##
      return(cp[findInterval(t, knots, rightmost.closed = TRUE)])
    }
    a_star_f <- function(t){sqrt(nu*tau_f(t)^2*(1+u_star_f(t)^2/nu)/nup)}
    ##
    fn1_star <- function(t){tau_f(t) * (1 + u_star_f(t)^2/nu + up_star_f(t)^2/(nu*tau_f(t)^2))^(-nu/2)/(2*pi)}
    fn2_star <- function(t){
      (up_star_f(t)/(2*pi*tau_f(t))) * (1+u_star_f(t)^2/nu)^(-nu/2 -1)  *
        (gamma(nup/2) * sqrt(nup*pi) * a_star_f(t) / gamma((nup+1)/2) ) *
        stats::pt(q = (up_star_f(t) / a_star_f(t) ), df=nup) 
    }
    fn3_star <- function(t){
      (up_star_f(t)/(2*pi*tau_f(t))) * (1+u_star_f(t)^2/nu)^(-nu/2 -1)  *
        (gamma(nup/2) * sqrt(nup*pi) * a_star_f(t) / gamma((nup+1)/2) ) *
        stats::pt(q = -(up_star_f(t) / a_star_f(t) ), df=nup)
    }
    ##
    intgr1_star <- sum(fn1_star(t=tt)) * diff(tt)[1]
    ##
    if(t0>0){intgr2_star <- sum(fn2_star(t=tt[tt <= t0])) * diff(tt)[1]}else{intgr2_star <- 0}
    if(t0<1){intgr3_star <- sum(fn3_star(t=tt[t0 <  tt])) * diff(tt)[1]}else{intgr3_star <- 0}
    ##
    optim_target <- c( (stats::pt(q=u_star_f(t0), lower.tail=FALSE, df = nu) + (intgr1_star + intgr2_star - intgr3_star)) - alpha.level/2 )
    band.eval    <- u_star_f(t=tt)
    ##
    return(list("optim_target"     = optim_target,
                "band.eval"        = band.eval, 
                "prob_t0"          = stats::pt(q=u_star_f(t0), lower.tail=FALSE, df = nu),
                "a_star"           = (intgr1_star + intgr2_star - intgr3_star) 
                ))
  }
  ##
  fu <- function(x){find_u(x)$optim_target}
  # fu <- Vectorize(fu); curve(fu, .Machine$double.eps, 1); abline(h=0)
  ##
  opt_res <- stats::uniroot(f = fu, interval = c(.Machine$double.eps, (alpha.level/2)), extendInt = "upX", tol=tol)$root
  band    <- find_u(opt_res)$band.eval * sqrt(diag.cov) 
  prob_t0 <- find_u(opt_res)$prob_t0
  a_star  <- find_u(opt_res)$a_star
  ##
  if(abs(prob_t0 + a_star - alpha.level/2) > ( (alpha.level/2) / 25 ) ){
    stop("Numeric equation solvers 'uniroot()' give too imprecise results. \nPlease, set the 'tol' argument smaller than the default value (<.Machine$double.eps^0.32).")
  }
  ##
  return(list("band"    = band,
              "t0"      = t0,
              "prob_t0" = prob_t0,
              "a_star"  = a_star))
}


# make_band_minW_z <- function(x, tau, diag.cov, N, conf.level=0.95, band_pen=1){
#   ##
#     alpha.level <- 1-conf.level
# 
#     lamb_t <- tau
#     ## spline functions u(t) and u'(t)
#     nknots  <- 5
#     knots   <- seq(0,1,length=nknots)
#     times   <- seq(0,1,length=length(lamb_t))
#     bspl    <- splines::bs(times, knots = knots[-c(1,nknots)], intercept = TRUE)
#     bspl_d  <- splines2::dbs(times, knots=knots[-c(1,nknots)], derivs=1, intercept = TRUE)
#     #bspl_d2 <- splines2::dbs(times, knots=knots[-c(1,nknots)], derivs=2, intercept = TRUE)
# 
#     u   <- function(coefs){exp(bspl   %*% coefs)           }
#     up  <- function(coefs){   (bspl_d %*% coefs) * u(coefs)}
#     ## function to minimize
#     eval_f    <- function(coefs){return(mean(u(coefs)^2) + band_pen *  mean(up(coefs)^2) )}
#     ## constraining function
#     eval_g_eq <- function(coefs){
#       return(
#       c(stats::pnorm(u(coefs)[1],lower.tail=FALSE)
#         + mean( lamb_t * exp(-u(coefs)^2/2 - up(coefs)^2/(2*lamb_t^2) )/(2*pi) )
#         + mean( up(coefs)/sqrt(2*pi) * exp(-u(coefs)^2/2 ) * stats::pnorm(-up(coefs)/lamb_t) )
#         - alpha.level/2)
#       )
#     }
#     ## determine starting value = "constant threshold"
#     a       <- mean(lamb_t)# int_0^1 tau
#     myfun_const <-function(c){
#       stats::pnorm(c,lower.tail = FALSE) + a * exp(-c^2/2)/(2*pi) - alpha.level/2
#     }
#     c_m    <- stats::uniroot(f = myfun_const, interval = c(0,10), extendInt="downX")$root
#     qnorm2 <- function(p){-stats::qnorm((1-p)/2)}
#     c_l    <- qnorm2(conf.level)
#     c_u    <- qnorm2(c(1-alpha.level/500))
#     ##
#     opt_result <- nloptr::nloptr(#x0 = runif(n = ncol(bspl), min = log(c_l), max = log(c_u)),#rep(log(c_m), ncol(bspl)),
#       x0        = rep(log(c_m), ncol(bspl)),
#       eval_f    = eval_f,
#       lb        = rep(log(c_l), ncol(bspl)),
#       ub        = rep(log(c_u), ncol(bspl)),
#       #eval_g_eq = eval_g_eq,
#       eval_g_ineq = eval_g_eq,
#       opts      = list("algorithm"         = c("NLOPT_LN_AUGLAG",
#                                                "NLOPT_LN_AUGLAG_EQ",
#                                                "NLOPT_GN_ISRES")[3],
#                        "xtol_rel"          = 1e-20,
#                        #"tol_constraints_eq"= 1e-5,
#                        "maxeval"           = 5000,
#                        "print_level"       = 0,
#                        "print_options_doc" = FALSE)
#     )
#     # nloptr.print.options()
#     # opt_result$objective
#     # opt_result$solution
#     # opt_result$solution - log(c_m)
#     # eval_g_eq(opt_result$solution)
# 
#     u_opt <- exp(bspl  %*% opt_result$solution)
# 
#     # plot(exp(bspl  %*% rep(log(c_m), ncol(bspl))),type="l", ylim=range(0,2.8))
#     # lines(u_opt)
#       band        <- u_opt * sqrt(diag.cov) / sqrt(N)
#       band_upper  <- x + band
#       band_lower  <- x - band
#       ##
#       result           <- cbind(band_lower,x,band_upper)
#       colnames(result) <- c("min_width_z_band_lower","x","min_width_z_band_upper")
#       return(result)
# }

# make.band.min_width <- function(x, tau, diag.cov, N, conf.level=0.95, band_pen=1){
#   alpha.level <- 1-conf.level
# 
#   lamb_t <- tau
# 
#   ## spline functions u(t) and u'(t)
#   nknots  <- 5
#   knots   <- seq(0,1,length=nknots)
#   times   <- seq(0,1,length=length(lamb_t))
#   mybs    <- splines::bs(times, knots = knots[-c(1,nknots)])
#   mybs_d  <- splines2::dbs(times, knots=knots[-c(1,nknots)], derivs=1)
#   mybs_d2 <- splines2::dbs(times, knots=knots[-c(1,nknots)], derivs=2)
# 
#   myfun <- function(xx){# xx: basis coefficients
#     ##
#     u   <- exp(mybs%*%xx)
#     up  <- (mybs_d %*%xx) * u
#     up2 <- (mybs_d2%*%xx) * u + (mybs_d%*%xx) * up
#     ##
#     myfun2 <- function(c, u=u, up=up){# c <- 10
#       c(pnorm(c*u[1],lower.tail=FALSE)
#         + mean( lamb_t * exp(-c^2*u^2/2 - c^2*up^2/(2*lamb_t^2) )/(2*pi) )
#         + mean( c*up/sqrt(2*pi) * exp(- c^2*u^2/2 ) * pnorm(-c*up/lamb_t) )
#         - alpha.level/2)
#     }
#     ##
#     c <- uniroot(f=myfun2,interval=c(.01,10), u=u, up=up, extendInt="downX")$root
#     ##
#     # return(mean(c^2*u^2) + band_pen *  mean((c^2 * up^2)))
#     #return(mean(c^2*u^2) + band_pen *  mean((c^2 * up^2)) )
#     return(mean(c^2*u^2) + band_pen * (mean((c^2 * up^2)) + mean((c^2 * up2^2)))/2 )
#   }
# 
#   # optimal basis coefficients
#   bs_coefs <- optim(par=rep(1,times=dim(mybs)[2]), fn = function(xx){myfun(xx)})$par
#   ##
#   my_u     <- exp(mybs  %*%bs_coefs)
#   my_up    <- exp(mybs_d%*%bs_coefs) * my_u
# 
#   # find the c
#   myfun2 <- function(c, my_u=my_u, my_up=my_up){
#     c(pnorm(c*my_u[1],lower.tail=FALSE)
#       + mean( lamb_t*exp(-c^2*my_u^2/2 - c^2*my_up^2/(2*lamb_t^2) )/(2*pi) )
#       + mean( c*my_up/sqrt(2*pi) * exp(- c^2*my_u^2/2 ) * pnorm(-c*my_up/lamb_t) )
#       - alpha.level/2)
#   }
#   c <- uniroot(f=myfun2,interval=c(.01,10), my_u=my_u, my_up=my_up, extendInt="downX")$root
# 
#   # the optimal threshold fct u
#   my_u        <- c * exp(mybs%*%bs_coefs)
#   band        <- my_u * sqrt(diag.cov) / sqrt(N)
#   band_upper  <- x + band
#   band_lower  <- x - band
#   ##
#   result           <- cbind(band_lower,x,band_upper)
#   colnames(result) <- c("min_width_z_band_lower","x","min_width_z_band_upper")
#   return(result)
# }
# 
