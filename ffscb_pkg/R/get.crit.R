#' @export
get.crit.supnorm.simple <- function(cov.m,n.sim,prob){ ### Finding Maximum of Normal - Supnorm criteria ----
  J <- length(diag(cov.m)) # number of points to estimate  ### Simple Maximum Simulation (Degras) ----
  X <- array(0,dim=c(J,n.sim))
  sd.v <- sqrt(diag(cov.m))
  cor.m <- cov.m / outer(sd.v,sd.v)
  #diag(cor.m) <- 1         # Set diagonal elements to 1
  eigen <- eigen(cor.m)
  eigen$values[eigen$values<.Machine$double.eps] <- 0 # handle numerical error
  sd.m.2 <- crossprod(t(eigen$vectors), diag(sqrt(eigen$values)))
  X <- crossprod(t(sd.m.2), matrix(rnorm(J*n.sim),ncol=n.sim))
  MaxX <- apply(abs(X),2,max)
  quantile(MaxX,prob)
}

#' @export
get.crit.Rz <- function(lambdas,prob){ ## Assigning log(1-alpha) proportional to lambda.
  if (!is.null(dim(lambdas)[2])) {
    lambdas <- eigen(lambdas)$values   # May Feed either eigenvalues or cov.matrix
    lambdas <- lambdas[lambdas>.Machine$double.eps]
  }
  qnorm2(exp(lambdas / sum(lambdas) * log(prob)))
}

#' @export
get.crit.Rz1 <- function(lambdas,prob){  ## Minimize Size ##
  if (!is.null(dim(lambdas)[2])) {
    lambdas=eigen(lambdas)$values
    lambdas=lambdas[lambdas>.Machine$double.eps]
  }
  lowerb <- sqrt(2*pi) * lambdas[1] * sf.f1(qnorm2(prob))
  upperb <- sqrt(2*pi) * lambdas[1] * sf.f1(qnorm2(prob^(1/length(lambdas))))
  M <- optim(c((lowerb+upperb)/2),fn=sf.diffP,lambdas=lambdas,prob=prob,
            method="Brent",lower=lowerb,upper=upperb)$par
  sapply(M / (sqrt(2*pi) * lambdas), sf.f1.inv)
}

#' @export
get.crit.Rz0 <- function(lambdas,prob){  ## Minimize Size ##
  if (!is.null(dim(lambdas)[2])) {
    lambdas <- eigen(lambdas)$values   # May Feed either eigenvalues or cov.matrix
    lambdas <- lambdas[lambdas>.Machine$double.eps]
  }
  #l <- length(lambdas)
  #qnorm2(exp(rep(1/l, l) * log(prob)))
  qnorm2(exp(lambdas^{0.37} / sum(lambdas^{0.37}) * log(prob)))
}

#' @export
get.crit.Rzs <- function(lambdas,prob,df){ qt2(pnorm2(get.crit.Rz(lambdas,prob)),df) }

#' @export
get.crit.Rz1s <- function(lambdas,prob,df){ qt2(pnorm2(get.crit.Rz1(lambdas,prob)),df) }
