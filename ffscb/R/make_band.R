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



# make_band_KR_z <- function(tau, diag.cov, conf.level=0.95){
#   alpha.level <- 1-conf.level
#   tt          <- seq(0,1,len=length(tau))
#   tau_01      <- sum(tau)*diff(tt)[1] # int_0^1 tau(t) dt
#   myfun       <- function(c){stats::pnorm(c,lower.tail = FALSE) + tau_01*exp(-c^2/2)/(2*pi)-alpha.level/2}
#   cstar       <- stats::uniroot(f = myfun,interval = c(.5,8))$root
#   band        <- rep(cstar, times=length(tau)) * sqrt(diag.cov) 
#   ##
#   return(band)
# }



# ## Confidence band
# make_band_KR_t <- function(tau, diag.cov, df, conf.level=0.95){
#   alpha.level <- 1-conf.level
#   tt          <- seq(0,1,len=length(tau))
#   tau_01      <- sum(tau)*diff(tt)[1] # int_0^1 tau(t) dt
#   myfun       <- function(c){stats::pt(c, lower.tail = FALSE, df=df) + tau_01*(1+c^2/df)^(-df/2)/(2*pi) - alpha.level/2}
#   cstar       <- stats::uniroot(f = myfun,interval = c(.5,8))$root
#   band        <- rep(cstar, times=length(tau)) * sqrt(diag.cov)
#   ##
#   return(band)
# }


# ## Prediction band
# make_prediction_band_KR_t <- function(x, tau, diag.cov, df, conf.level=0.95){
#   alpha.level <- 1-conf.level
#   tt          <- seq(0,1,len=length(tau))
#   tau_01      <- sum(tau)*diff(tt)[1] # int_0^1 tau(t) dt
#   myfun       <- function(c){stats::pt(c, lower.tail = FALSE, df=df) + tau_01*(1+c^2/df)^(-df/2)/(2*pi) - alpha.level/2}
#   cstar       <- stats::uniroot(f = myfun,interval = c(.5,8))$root
#   band        <- rep(cstar, times=length(tau)) * sqrt(diag.cov * df) * sqrt(1 + (1/(df + 1 )))
#   ##
#   band_m           <- cbind(x, x + band, x - band)
#   colnames(band_m) <- c("x", paste0("FFSCB.z.u.", conf.level), paste0("FFSCB.z.l.", conf.level))
#   ##
#   return(band_m)
# }


#' Fast 'n' fair simultaneous confidence band (Gaussian)
#'
#' @param x Functional parameter estimate. 
#' @param diag.cov.x Diagonal of Cov(x), in which x is the functional estimator (for instance, the covariance function of the empirical mean function).
#' @param tau Pointwise standard deviation of the standardized and differentiated sample functions. Can be estimated by tau_fun().
#' @param conf.level confidence level (default: 0.95)
#' @param n_int Number of equidistant intervals over which the multiple testing component of the type-I error rate (1-conf.level) is distributed uniformly.
#' @references Liebl, D. and Reimherr, M. (2021+). Fast and fair simultaneous confidence bands.
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
#' band <- make_band_FFSCB_z(x=hat.mu, diag.cov.x=diag(hat.cov.mu), tau=hat.tau,
#'                        conf.level  = 0.95)
#' matplot(y=band[,2:3], x=grid, lty=2)
#' lines(x=grid, y=band[,1], lty=1)
#' @export
make_band_FFSCB_z <- function(x, diag.cov.x, tau, conf.level=0.95, n_int=4){
  ##
  if(any(tau < 0.005)){warning("This method may not work if tau(t) is too small.")}
  ##
  band             <- .make_band_FFSCB_z(tau=tau, diag.cov=diag.cov.x, conf.level=conf.level, n_int=n_int)
  band_m           <- cbind(x, x + band, x - band)
  colnames(band_m) <- c("x", paste0("FFSCB.z.u.", conf.level), paste0("FFSCB.z.l.", conf.level))
  ##
  return(band_m)
}


.make_band_FFSCB_z <- function(tau, diag.cov, conf.level=0.95, n_int=4){
  ##
  ## location to check whether process started uncrossed
  ##t0          <- 0 
  ## increase the default accuracy for uniroot() (.Machine$double.eps^0.25)
  tol         <- .Machine$double.eps^0.35 
  ##
  alpha       <- 1-conf.level
  tt          <- seq(0,1,len=length(tau))
  tau_v       <- tau
  ## linear interpolation of all tau values for easy integration
  tau_f       <- function(t){stats::approx(x = seq(0,1,len=length(tau)), y = tau, xout=t)$y}
  knots       <- seq(0,1,len=(n_int + 1))
  ##
  if(!is.numeric(n_int)){
    stop("n_int must be a stricitly positive integer value (n_int=1,2,3,...)")
  }
  if(n_int <=0){
    stop("n_int must be a stricitly positive integer value (n_int=1,2,3,...)")
  }
  if(n_int %% 1 != 0){
    stop("n_int must be a stricitly positive integer value (n_int=1,2,3,...)")
  }
  ##
  if(n_int == 1){# Case n_int=1 == constant band == Kac-Rice Band
    tau01      <- sum(tau_v)*diff(tt)[1] # int_0^1 tau(t) dt
    myfun1     <- function(c1){stats::pnorm(q=c1, lower.tail=F)+exp(-c1^2/2)*tau01/(2*pi)-(alpha/2)}
    const_band <- stats::uniroot(f = myfun1, interval = c(0,10), extendInt="downX")$root
    band       <- const_band * sqrt(diag.cov) 
    return(band)
  }
  ##
  ## Remainder part considers the case of n_int >1
  ##
  ## Define piecewise linear (pwl) function 'ufun' with derivative=0 at t0.
  c_v         <- numeric(n_int) # coeficients of the pwl-function
  const_int   <- 1
  fct_body    <- paste0("c_v[",const_int,"]")
  for(j in (const_int+1):n_int){fct_body <- paste0(fct_body, "+ c_v[",j,"]*pmax(t - knots[",j,"],0)")}
  ufun        <- function(t, c_v, knots){}
  body(ufun)  <- parse(text=fct_body)
  ##
  ## Determine the initial value of u
  tau_init <- sum(tau_v[knots[const_int] <= tt & tt <= knots[const_int+1]])*diff(tt)[1]
  myfun1   <- function(c1){
    stats::pnorm(-c1) +
      c(exp(-c1^2/2) * tau_init / (2*pi) - (alpha/2)/n_int)
  }
  ## curve(myfun1, 0,5); abline(h=0)
  ##
  c_v[const_int] <- stats::uniroot(f = myfun1, interval = c(0,10), extendInt = "downX", tol = tol)$root
  ##
  for(j in (const_int+1):n_int){# j <- (const_int+1)
    myfunj <- function(cj){
      ## sum(c_v_sum) == u'(t) for all t in the (j-1)th interval
      ## sum(c(c_v_sum,c_j)) == u'(t) for all t in the jth interval
      if(j==(const_int+1)){c_v_sum <- 0}else{c_v_sum <- c_v[(const_int+1):(j-1)]}
      ##
      ufun_j <- function(t,cj){
        ufun(t=t,
             c_v=c(#rep(0,times=(const_int-1)),
                   c_v[const_int:(j-1)], cj, rep(0, times=(n_int-j))),
             knots=knots)
      }
      fn1    <- function(t,cj){
        (tau_f(t)/(2*pi)) *
          exp(-ufun_j(t,cj)^2/2) *
          exp(-sum(c(c_v_sum,cj))^2/(2*tau_f(t)^2))
      }
      fn2    <- function(t,cj){
        sum(c(c_v_sum,cj))/sqrt(2*pi) * 
          exp(-ufun_j(t,cj)^2/2) * 
          stats::pnorm( sum(c(c_v_sum, cj))/tau_f(t))
      }
      fn3    <- function(t,cj){
        sum(c(c_v_sum,cj))/sqrt(2*pi) *
          exp(-ufun_j(t,cj)^2/2) *
          stats::pnorm(-sum(c(c_v_sum, cj))/tau_f(t))
      }
      intgr1 <- sum(fn1(t=tt[knots[j] < tt & tt <= knots[j+1]], cj=cj)) * diff(tt)[1]
      if(j %% 2 != 0){# odd j
        intgr2 <- sum(fn2(t=tt[knots[j] < tt & tt <= knots[j+1]], cj=cj)) * diff(tt)[1]
      }else{intgr2 <- 0}
      if(j %% 2 == 0){# even j
        intgr3 <- sum(fn3(t=tt[knots[j] < tt & tt <= knots[j+1]], cj=cj)) * diff(tt)[1]
      }else{intgr3 <- 0}
      ##
      ## res    <- c(stats::pnorm(-c_v[1])/n_int + intgr1 - intgr3 - (alpha/2)/n_int)
      if(j %% 2 == 0){# even j
        res <- c(stats::pnorm(-ufun(knots[j], c_v = c_v, knots = knots)) + intgr1 + intgr2 - intgr3 - (alpha/2)/n_int) # '-ufun' since we need upper-tail and can use symmetriy of normal
      }
      if(j %% 2 != 0){# odd j
        res <- c(stats::pnorm(-ufun_j(t=knots[j+1], cj))                 + intgr1 + intgr2 - intgr3 - (alpha/2)/n_int) # '-ufun' since we need upper-tail and can use symmetriy of normal
      }
      ##
      return(res)
    }
    # myfunj <- Vectorize(myfunj); curve(myfunj, -50,50); abline(h=0); myfunj(c_v[j])
    c_v[j] <- stats::uniroot(f = myfunj, interval = c(-10,10), extendInt = "downX", tol = tol)$root
  }
  ##
  band.eval <- ufun(t=tt, c_v=c_v, knots=knots) # plot(x=tt,y=band.eval, type="l", main="z")
  band      <- band.eval * sqrt(diag.cov) # plot(x=tt,y=band, type="l", main="z")
  ##
  return(band)
}



#' Fast 'n' fair simultaneous confidence band (t-distr)
#'
#' @param x Functional parameter estimate. 
#' @param diag.cov.x Diagonal of Cov(x), in which x is the functional estimator (for instance, the covariance function of the empirical mean function).
#' @param tau Pointwise standard deviation of the standardized and differentiated sample functions. Can be estimated by tau_fun().
#' @param df Degrees of freedom 
#' @param conf.level confidence level (default: 0.95)
#' @param n_int Number of equidistant intervals over which the multiple testing component of the type-I error rate (1-conf.level) is distributed uniformly.
#' @references Liebl, D. and Reimherr, M. (2021+). Fast and fair simultaneous confidence bands.
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
#' band <- make_band_FFSCB_t(x=hat.mu, diag.cov.x=diag(hat.cov.mu), tau=hat.tau,
#'                        df = N-1, conf.level  = 0.95)
#' matplot(y=band[,2:3], x=grid, lty=2)
#' lines(x=grid, y=band[,1], lty=1)
#' @export
make_band_FFSCB_t <- function(x, diag.cov.x, tau, df, conf.level=0.95, n_int=4){
  ##
  if(any(tau < 0.005)){warning("This method may not work if tau(t) is too small.")}
  ##
  if(df <= 101){
    band       <- .make_band_FFSCB_t(tau=tau, diag.cov=diag.cov.x, df=df, conf.level=conf.level, n_int=n_int)
  }else{
    band       <- .make_band_FFSCB_z(tau=tau, diag.cov=diag.cov.x,        conf.level=conf.level, n_int=n_int)
  }
  band_m       <- cbind(x, x + band, x - band)
  colnames(band_m) <- c("x", paste0("FFSCB.t.u.", conf.level), paste0("FFSCB.t.l.", conf.level))
  ##
  return(band_m)
}



.make_band_FFSCB_t <- function(tau, diag.cov, df, conf.level=0.95, n_int=4){
  ##
  ## location to check whether process started uncrossed
  #t0          <- 0 
  ## increases the default accuracy for uniroot() (.Machine$double.eps^0.25)
  tol         <- .Machine$double.eps^0.35 
  ##
  alpha       <- 1-conf.level
  tt          <- seq(0,1,len=length(tau))
  tau_v       <- tau
  nu          <- df
  nup         <- nu+1
  tau_f       <- function(t){stats::approx(x = seq(0,1,len=length(tau)), y = tau, xout=t)$y}
  knots       <- seq(0,1,len=(n_int + 1))
  ##
  if(!is.numeric(n_int)){
    stop("n_int must be a stricitly positive integer value (n_int=1,2,3,...)")
  }
  if(n_int <=0){
    stop("n_int must be a stricitly positive integer value (n_int=1,2,3,...)")
  }
  if(n_int %% 1 != 0){
    stop("n_int must be a stricitly positive integer value (n_int=1,2,3,...)")
  }
  ##
  if(n_int == 1){# Case n_int=1 == constant band == Kac-Rice Band
    tau01      <- sum(tau_v)*diff(tt)[1] # int_0^1 tau(t) dt
    myfun1     <- function(c1){c(stats::pt(q=c1, lower.tail=FALSE, df = nu)+(tau01/(2*pi))*(1+c1^2/nu)^(-nu/2)-(alpha/2))}
    const_band <- stats::uniroot(f = myfun1, interval = c(0,10), extendInt="downX")$root
    const_band <- const_band * sqrt(diag.cov) 
    return(const_band)
  }
  ##
  ## Remainder part considers the case of n_int >1
  ##
  ## Define piecewise linear (pwl) function 'ufun' with derivative=0 at t0.
  c_v         <- numeric(n_int) # coeficients of the pwl-function
  const_int   <- 1
  fct_body    <- paste0("c_v[",const_int,"]")
  for(j in (const_int+1):n_int){fct_body <- paste0(fct_body, "+ c_v[",j,"]*pmax(t - knots[",j,"],0)")}
  ufun        <- function(t, c_v, knots){}
  body(ufun)  <- parse(text=fct_body)
  ##
  ## Determine the initial value of u
  tau_init <- sum(tau_v[knots[const_int] <= tt & tt <= knots[const_int+1]])*diff(tt)[1]
  myfun1         <- function(c1){
    stats::pt(q=-c1, df = nu) +
      c((tau_init/(2*pi)) * (1+c1^2/nu)^(-nu/2) - (alpha/2)/n_int)
  }
  ## curve(myfun1, 0,1); abline(h=0)
  ##
  c_v[const_int] <- stats::uniroot(f = myfun1, interval = c(0,10), extendInt = "downX", tol = tol)$root
  ##
  for(j in (const_int+1):n_int){# j <- const_int+1
    myfunj <- function(cj){
      ##
      ## sum(c_v_sum) == u'(t) for 0 <= t < (j-1)th interval
      if(j==(const_int+1)){c_v_sum <- 0}else{c_v_sum <- c_v[(const_int+1):(j-1)]}
      ##
      ufun_j <- function(t,cj){
        ufun(t=t,
             c_v=c(c_v[const_int:(j-1)],cj,rep(0, times=(n_int-j))),
             knots=knots)
      }
      afun_j <- function(t,cj){sqrt(nu*tau_f(t)^2*(1+ufun_j(t,cj)^2/nu)/nup)}
      ##
      fn1    <- function(t,cj){
        tau_f(t) * (1 + ufun_j(t,cj)^2/nu +
                      sum(c(c_v_sum,cj))^2/(nu*tau_f(t)^2))^(-nu/2) / (2*pi)
      }
      fn2 <- function(t,cj){
                   (sum(c(c_v_sum,cj))/(2*pi*tau_f(t))) * (1+ufun_j(t,cj)^2/nu)^(-nu/2 -1)  *
                     (gamma(nup/2) * sqrt(nup*pi) * afun_j(t,cj) / gamma((nup+1)/2) ) *
                     stats::pt(q = (sum(c(c_v_sum,cj)) / afun_j(t,cj) ), df=nup) 
      }
      fn3    <- function(t,cj){
        (sum(c(c_v_sum,cj))/(2*pi*tau_f(t))) * 
          (1+ufun_j(t,cj)^2/nu)^(-nu/2 -1)  *
          (gamma(nup/2) * sqrt(nup*pi) * afun_j(t,cj) / gamma((nup+1)/2) ) *
          stats::pt(q = (-sum(c(c_v_sum,cj)) / afun_j(t,cj) ), df=nup) 
      }
      intgr1 <- sum(fn1(t=tt[knots[j] < tt & tt <= knots[j+1]],cj=cj)) * diff(tt)[1]
      if(j %% 2 != 0){# odd j
      intgr2 <- sum(fn2(t=tt[knots[j] < tt & tt <= knots[j+1]],cj=cj)) * diff(tt)[1]
      }else{intgr2 <- 0}
      if(j %% 2 == 0){# even j
      intgr3 <- sum(fn3(t=tt[knots[j] < tt & tt <= knots[j+1]],cj=cj)) * diff(tt)[1]
      }else{intgr3 <- 0}
      ##
      #res    <- c(stats::pt(q=-c_v[1], df = nu) + intgr1 - intgr3 - (alpha/2)/n_int)
      if(j %% 2 == 0){# even j
        #res <- c(stats::pnorm(-ufun(knots[j], c_v = c_v, knots = knots)) + intgr1 + intgr2 - intgr3 - (alpha/2)/n_int)
        res <- c(stats::pt(-ufun(knots[j], c_v = c_v, knots = knots), df=df) + intgr1 + intgr2 - intgr3 - (alpha/2)/n_int) # '-ufun' since we need upper-tail and can use symmetriy of normal
      }
      if(j %% 2 != 0){# odd j
        res <- c(stats::pt(-ufun_j(t=knots[j+1], cj), df=df)                 + intgr1 + intgr2 - intgr3 - (alpha/2)/n_int) # '-ufun' since we need upper-tail and can use symmetriy of normal
      }
      return(res)
    } 
    # myfunj <- Vectorize(myfunj); curve(myfunj, -10,10); abline(h=0)
    # 2.7206646  0.5871968  0.8107724 -0.7136544
    c_v[j] <- stats::uniroot(f = myfunj, interval = c(-10,10), extendInt = "downX", tol = tol)$root
  }
  ##
  band.eval <- ufun(t=tt, c_v=c_v, knots=knots) # plot(y=band.eval,x=tt, type="l")
  band      <- band.eval * sqrt(diag.cov)       # plot(y=band,x=tt, type="l", main="t")
  ##
  return(band)
}
