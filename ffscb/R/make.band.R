#' @export
make.band.BEc <- function(eigen, conf.level, fd.eval.grid.size=200){
  alpha.level <- 1-conf.level
  pc.to.use   <- sum(eigen$values > .Machine$double.eps)
  c.square    <- sqrt(eigen$values[1:pc.to.use])
  weights     <- eigen$values[1:pc.to.use]/c.square
  xi          <- ffscb:::get.schisq.q.gamma(weights,conf.level) ## Approximate Quantile of Weighted Sum of Chi-square by Gamma
  if (inherits(eigen,"pca.fd") | inherits(eigen,"eigen.fd")) {
    evalgrid      <- ffscb::make.grid(p=fd.eval.grid.size, rangevals=eigen$harmonics$basis$rangeval)
    eigen$vectors <- fda::eval.fd(evalgrid,eigen$harmonics)
  }
  band.eval <- sqrt(apply(t(eigen$vectors[,1:pc.to.use]^2) * c.square * xi,2,sum))
  if (inherits(eigen,"pca.fd") | inherits(eigen,"eigen.fd")) {
    return(fda::Data2fd(evalgrid,band.eval,basisobj=eigen$harmonics$basis)) # return as fd object
  } else return(band.eval)                                             # return as vector
}


#' @export
make.band.Bs <- function(cov, conf.level, sim.size=10000, fd.eval.grid.size=200){
  if (inherits(cov,"bifd")) {
    evalgrid <- ffscb::make.grid(p=fd.eval.grid.size, rangevals=cov$sbasis$rangeval)
    cov.m    <- fda::eval.bifd(evalgrid,evalbrid,cov) 
  } else {
    cov.m <- cov
  }
  crit.Bs    <- ffscb:::get.crit.supnorm.simple(cov.m,n.sim=sim.size,p=conf.level)
  band.eval  <- sqrt(diag(cov.m)) * crit.Bs
  if (inherits(cov,"bifd")) {
    return(fda::Data2fd(evalgrid,band.eval,basisobj=cov$sbasis))
  } else return(band.eval)
}

#' @export
make.band.naive.t <- function(cov, conf.level, df, fd.eval.grid.size=200){
  if (inherits(cov,"bifd")) {
    evalgrid <- ffscb::make.grid(p=fd.eval.grid.size, rangevals=cov$sbasis$rangeval)
    cov.m    <- fda::eval.bifd(evalgrid,evalgrid,cov) 
  } else {
    cov.m <- cov
  }
  band.eval  <- ffscb:::qt2(conf.level,df) * sqrt(diag(cov.m))
  if (inherits(cov,"bifd")) {
    return(fda::Data2fd(evalgrid,band.eval,basisobj=cov$sbasis))
  } else return(band.eval)
}

#' @export
make.band.KR.z <- function(tau, conf.level){
  alpha.level <- 1-conf.level
  tt          <- seq(0,1,len=length(tau))
  tau_01      <- sum(tau)*diff(tt)[1] # int_0^1 tau(t) dt
  myfun       <- function(c){stats::pnorm(c,lower.tail = FALSE)+tau_01*exp(-c^2/2)/(2*pi)-alpha.level/2}
  cstar       <- stats::uniroot(f = myfun,interval = c(.5,8))$root
  band.eval   <- rep(cstar, times=length(tau))
  return(band.eval)
}

#' @export
make.band.KR.t <- function(tau, conf.level, N){
  alpha.level <- 1-conf.level
  nu          <- N-1
  tt          <- seq(0,1,len=length(tau))
  tau_01      <- sum(tau)*diff(tt)[1] # int_0^1 tau(t) dt
  myfun       <- function(c){stats::pt(c, lower.tail = FALSE, df=nu)+tau_01*(1+c^2/nu)^(-nu/2)/(2*pi) - alpha.level/2}
  cstar       <- stats::uniroot(f = myfun,interval = c(.5,8))$root
  band.eval   <- rep(cstar, times=length(tau))
  return(band.eval)
}


#' @export
make.band.FFSCB.z <- function(tau, t0, conf.level, n_int){
  ##
  alpha.level <- 1-conf.level
  tt          <- seq(0,1,len=length(tau))
  tau_v       <- stats::spline(   x = seq(0,1,len=length(tau)), y = tau, method = "natural", xout = tt)$y
  tau_f       <- stats::splinefun(x = seq(0,1,len=length(tau)), y = tau, method = "natural")
  knots       <- seq(0,1,len=(n_int + 1))
  tol         <- .Machine$double.eps^0.25 # default accuracy for optimize() 
  if(is.null(t0)){t0 <- tt[which.max(tau_v)]}else{t0 <- tt[findInterval(t0,tt)]}
  ##
  if(n_int == 1){# Case n_int=1 == constant band == Kac-Rice Band
    tau01      <- sum(tau_v)*diff(tt)[1] # int_0^1 tau(t) dt
    myfun1     <- function(c1){stats::pnorm(q=c1, lower.tail=F)+exp(-c1^2/2)*tau01/(2*pi)-(alpha.level/2)}
    const_band <- stats::uniroot(f = myfun1, interval = c(0,10), extendInt="downX")$root
    return(rep(const_band, times=length(tau)))
  }
  ## Define piecewise linear (pwl) function 'ufun' with derivative=0 at t0.
  c_v         <- numeric(n_int) # coeficients of the pwl-function
  const_int   <- min(findInterval(t0, knots), n_int)
  fct_body    <- paste0("c_v[",const_int,"]")
  if(const_int >     1){for(j in (const_int-1):1    ){fct_body <- paste0("c_v[",j,"]*pmin(t - knots[",j,"],0) +", fct_body)}}
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
      myfun1         <- function(c){c(exp(-c^2/2) * tau_init / (2*pi) - (alpha.aux/2)/n_int)^2}
      c_v[const_int] <- stats::optimize(f = myfun1, interval = c(0,10), tol = tol)$minimum
    }
    ##
    if(const_int > 1){
      for(j in (const_int-1):1){
        myfunj <- function(cj){
          ##
          if(j==(const_int-1)){c_v_sum <- 0}else{c_v_sum <- c_v[(const_int-1):(j+1)]}# c_v_sum == u' of the preceeding interval
          ##
          ufun_j <- function(t,cj){ufun(t=t,c_v=c(rep(0,times=(j-1)),cj,c_v[(j+1):const_int],rep(0, times=(n_int-const_int))), knots=knots)}
          ##
          fn1    <- function(t,cj){(tau_f(t)/(2*pi)) * exp(-ufun_j(t,cj)^2/2) * exp(-sum(c(c_v_sum,cj))^2/(2*tau_f(t)^2))}
          fn2    <- function(t,cj){sum(c(c_v_sum,cj))/sqrt(2*pi) * exp(-ufun_j(t,cj)^2/2) * stats::pnorm( sum(c(c_v_sum, cj))/tau_f(t))}
          intgr1 <- sum(fn1(t=tt[knots[j] < tt & tt <= knots[j+1]], cj=cj)) * diff(tt)[1]
          intgr2 <- sum(fn2(t=tt[knots[j] < tt & tt <= knots[j+1]], cj=cj)) * diff(tt)[1]
          ##
          res    <- c(intgr1 + intgr2 - (alpha.aux/2)/n_int)#^2
          return(res)
        }
        # c_v[j] <- stats::optimize(f = myfunj, interval = c(-4,4), tol = tol)$minimum
        c_v[j] <- stats::uniroot(f = myfunj, interval = c(-4,4), extendInt = "upX", tol = tol)$root
      }
    }
    if(const_int < n_int){
      for(j in (const_int+1):n_int){# j <- (const_int+1)
        myfunj <- function(cj){
          ##
          if(j==(const_int+1)){c_v_sum <- 0}else{c_v_sum <- c_v[(const_int+1):(j-1)]}# c_v_sum == u'(t) for 0 <= t < (j-1)th interval
          ##
          ufun_j <- function(t,cj){ufun(t=t,c_v=c(c_v[1:(j-1)],cj,rep(0, times=(n_int-j))), knots=knots)}
          ##
          fn1    <- function(t,cj){(tau_f(t)/(2*pi)) * exp(-ufun_j(t,cj)^2/2) * exp(-sum(c(c_v_sum,cj))^2/(2*tau_f(t)^2))}
          fn3    <- function(t,cj){sum(c(c_v_sum,cj))/sqrt(2*pi) * exp(-ufun_j(t,cj)^2/2) * stats::pnorm(-sum(c(c_v_sum, cj))/tau_f(t))}
          intgr1 <- sum(fn1(t=tt[knots[j] < tt & tt <= knots[j+1]], cj=cj)) * diff(tt)[1]
          intgr3 <- sum(fn3(t=tt[knots[j] < tt & tt <= knots[j+1]], cj=cj)) * diff(tt)[1]
          ##
          res    <- c(intgr1 - intgr3 - (alpha.aux/2)/n_int)#^2
          return(res)
        }
        #c_v[j] <- stats::optimize(f = myfunj, interval = c(-4,4), tol = tol)$minimum
        c_v[j] <- stats::uniroot(f = myfunj, interval = c(-4,4), extendInt = "downX", tol = tol)$root
      }
    }
    ## 
    u_star_f  <- function(t){return(ufun(t=t,c_v=c_v,knots=knots))}
    up_star_f <- function(t){
      kn   <- knots[-length(knots)]; cp <- c_v; cp[1] <- 0
      kn_m <- matrix(kn,                     nrow=length(kn), ncol=length(t))
      cp_m <- matrix(cp,                     nrow=length(kn), ncol=length(t))
      t_m  <- matrix(rep(t,each=length(kn)), nrow=length(kn), ncol=length(t))
      return(colSums(cp_m * (kn_m<t_m)))
    }
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
    optim_target <- c( (stats::pnorm(q=u_star_f(t0), lower.tail=F) + intgr1_star + intgr2_star - intgr3_star) - alpha.level/2)^2
    band.eval    <- u_star_f(t=seq(0,1,len=length(tau)))
    ##
    return(list("optim_target" = optim_target,
                "band.eval"    = band.eval))
  }
  ##
  opt_res <- stats::optimize(f = function(x){find_u(alpha.aux = x)$optim_target}, interval = c(0,alpha.level), tol=tol)$minimum
  ##
  return(find_u(opt_res)$band.eval)
}


#' @export
make.band.FFSCB.t <- function(tau, t0, conf.level, N, n_int){
  ##
  alpha.level <- 1-conf.level
  tt          <- seq(0,1,len=length(tau))
  tau_v       <- stats::spline(   x = seq(0,1,len=length(tau)), y = tau, method = "natural", xout = tt)$y
  tau_f       <- stats::splinefun(x = seq(0,1,len=length(tau)), y = tau, method = "natural")
  knots       <- seq(0,1,len=(n_int + 1))
  tol         <- .Machine$double.eps^0.25 # default accuracy for optimize() 
  nu          <- N-1
  nup         <- nu+1
  if(is.null(t0)){t0 <- tt[which.max(tau_v)]}else{t0 <- tt[findInterval(t0,tt)]}
  ##
  if(n_int == 1){# Case n_int=1 == constant band == Kac-Rice Band
    tau01      <- sum(tau_v)*diff(tt)[1] # int_0^1 tau(t) dt
    myfun1     <- function(c1){stats::pnorm(q=c1, lower.tail=F)+exp(-c1^2/2)*tau01/(2*pi)-(alpha.level/2)}
    const_band <- stats::uniroot(f = myfun1, interval = c(0,10), extendInt="downX")$root
    return(rep(const_band, times=length(tau)))
  }
  ## Define piecewise linear (pwl) function 'ufun' with derivative=0 at t0.
  c_v         <- numeric(n_int) # coeficients of the pwl-function
  const_int   <- min(findInterval(t0, knots), n_int)
  fct_body    <- paste0("c_v[",const_int,"]")
  if(const_int >     1){for(j in (const_int-1):1    ){fct_body <- paste0("c_v[",j,"]*pmin(t - knots[",j,"],0) +", fct_body)}}
  if(const_int < n_int){for(j in (const_int+1):n_int){fct_body <- paste0(fct_body, "+ c_v[",j,"]*pmax(t - knots[",j,"],0)")}}
  ufun        <- function(t, c_v, knots){}
  body(ufun)  <- parse(text=fct_body)
  ##
  find_u <- function(alpha.aux){
    ##
    ## Determine the initial value of u
    tau_init <- sum(tau_v[knots[const_int] <= tt & tt <= knots[const_int+1]])*diff(tt)[1]
    if(nu*(( (pi*alpha.aux) / ( tau_init *n_int) )^(-2/nu) -1) >0){# if possible use the analytic solution
      c_v[const_int] <- sqrt(nu*(( (pi*alpha.aux) / ( tau_init *n_int) )^(-2/nu) -1))
    }else{
      myfun1         <- function(c1){c((tau_init/(2*pi))*(1+c1^2/nu)^(-nu/2)-(alpha.aux/2)/n_int)^2}
      c_v[const_int] <- stats::optimize(f = myfun1, interval = c(0,10), tol = tol)$minimum
    }
    ##
    if(const_int > 1){
      for(j in (const_int-1):1){
        myfunj <- function(cj){
          ##
          if(j==(const_int-1)){c_v_sum <- 0}else{c_v_sum <- c_v[(const_int-1):(j+1)]}# c_v_sum == u' of the preceeding interval
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
                    sum(fn2(t=tt[knots[j] < tt & tt <= knots[j+1]],cj=cj)) * diff(tt)[1]) - (alpha.aux/2)/n_int )#^2
          return(res)
        } 
        #c_v[j] <- stats::optimize(f = myfunj, interval = c(-4,4), tol = tol)$minimum
        c_v[j] <- stats::uniroot(f = myfunj, interval = c(-4,4), extendInt = "upX", tol = tol)$root
      }
    }
    if(const_int < n_int){
      for(j in (const_int+1):n_int){
        myfunj <- function(cj){
          ##
          if(j==(const_int+1)){c_v_sum <- 0}else{c_v_sum <- c_v[(const_int+1):(j-1)]}# c_v_sum == u'(t) for 0 <= t < (j-1)th interval
          ##
          ufun_j <- function(t,cj){ufun(t=t,c_v=c(c_v[1:(j-1)],cj,rep(0, times=(n_int-j))), knots=knots)}
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
                    sum(fn3(t=tt[knots[j] < tt & tt <= knots[j+1]],cj=cj)) * diff(tt)[1]) - (alpha.aux/2)/n_int )#^2
          return(res)
        } 
        #c_v[j] <- stats::optimize(f = myfunj, interval = c(-4,4), tol = tol)$minimum
        c_v[j] <- stats::uniroot(f = myfunj, interval = c(-4,4), extendInt = "downX", tol = tol)$root
      }
    }
    ## 
    u_star_f  <- function(t){return(ufun(t=t,c_v=c_v,knots=knots))}
    up_star_f <- function(t){
      kn   <- knots[-length(knots)]; cp <- c_v; cp[1] <- 0
      kn_m <- matrix(kn,                     nrow=length(kn), ncol=length(t))
      cp_m <- matrix(cp,                     nrow=length(kn), ncol=length(t))
      t_m  <- matrix(rep(t,each=length(kn)), nrow=length(kn), ncol=length(t))
      return(colSums(cp_m * (kn_m<t_m)))
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
    optim_target <- c( (stats::pt(q=u_star_f(t0), lower.tail=F, df = nu) + intgr1_star + intgr2_star - intgr3_star) - alpha.level/2)^2
    band.eval    <- u_star_f(t=seq(0,1,len=length(tau)))
    ##
    return(list("optim_target" = optim_target,
                "band.eval"    = band.eval))
  }
  ##
  opt_res <- stats::optimize(f = function(x){find_u(alpha.aux = x)$optim_target}, interval = c(0,alpha.level), tol=tol)$minimum
  ##
  return(find_u(opt_res)$band.eval)
}

