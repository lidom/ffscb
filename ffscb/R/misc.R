#' Make grid
#'
#' @param p Length of the grid 
#' @param rangevals Endpoints of the grid 
#' @export
make_grid <- function(p=100,rangevals=c(0,1)){seq(rangevals[1],rangevals[2],len=p)}


#' Make discretized covariance function
#'
#' @param cov.f Covariance function
#' @param grid Evaluation-grid 
#' @param cov.f.params Parameters of the covariance function
#' @export
make_cov_m <- function(cov.f=covf.st.matern, grid, cov.f.params=NULL){  ### Make cov. matrix from cov. function.
  grid.size <- length(grid)
  cov.m <- matrix(0,nrow=grid.size,ncol=grid.size)
  if (is.null(cov.f.params)) {
    for (i in (1:grid.size)){
      cov.m[i,]=sapply(grid, cov.f, x1=grid[i])
    }
  }
  else{
    for (i in (1:grid.size)){
      cov.m[i,]=sapply(grid, cov.f, x1=grid[i], params=cov.f.params)
    }
  }
  return(cov.m)
}

#' Make sample (for simulation)
#'
#' @param mean.v Mean-Vector (discretized mean function)
#' @param cov.m Covariance-Matrix (discretized covariance function)
#' @param N Number of functions
#' @param dist Distribution
#' @param ... further parameters for dist
#' @export
make_sample <- function(mean.v,cov.m,N,dist="rnorm",...){
  p <- length(mean.v)
  if (p != dim(cov.m)[1] | p != dim(cov.m)[2]) stop("Dimensions of mean vector and cov. matrix do not match")
  dist        <- get(dist, mode <- "function", envir <- parent.frame())
  Z           <- matrix(dist(N*p,...),nrow=p,ncol=N)
  eigen.cov.m <- eigen(cov.m);
  eigen.cov.m$values[eigen.cov.m$values<0] <- 0
  X           <- crossprod( t(eigen.cov.m$vectors), crossprod(diag(sqrt(eigen.cov.m$values)), Z)) + mean.v
  rownames(X) <- names(mean.v)
  colnames(X) <- paste0("x",c(1:N))
  return(X)
}


#' Locate crossings
#'
#' @param x_vec Vector which might cross the given threshold
#' @param threshold Threshold-vector
#' @param type locate 'up' or 'down' crossings?
#' @export
locate_crossings <- function(x_vec, threshold, type=c("up", "down")){
  ##
  if(!any(type==c("up", "down"))){"'type' must be 'up' or 'down'."}
  if(type=="up"){  check <-  x_vec >  threshold}
  if(type=="down"){check <- -x_vec > -threshold}
  ##
  names(check) <- 1:length(check)
  check        <- diff(check)
  res          <- as.numeric(names(check)[check==1])
  return(res)
}


qt2 <- function(p,df){-stats::qt((1-p)/2,df)}

get.req.n.pc <- function(proportion, lambdas){ # Finds how many pcs are required to achieve desired variance proportion(prop. as 1-10^(-x))
  satisfied  <- (proportion <= cumsum(lambdas) / sum(lambdas))
  min(c(1:length(lambdas))[satisfied])
}


get.schisq.q.gamma <- function(weights,prob){  # prob as 1-alpha
  meanvalue  <- sum(weights)
  varvalue   <- sum(2*weights^2)
  gammascale <- varvalue / meanvalue   # 's'
  gammashape <- meanvalue / gammascale # 'a'
  stats::qgamma(p=prob,shape=gammashape,scale=gammascale,lower.tail=TRUE)
}


get.crit.supnorm.simple <- function(cov.m,n.sim,prob){ ### Finding Maximum of Normal - Supnorm criteria ----
  J <- length(diag(cov.m)) # number of points to estimate  ### Simple Maximum Simulation (Degras) ----
  X <- array(0,dim=c(J,n.sim))
  sd.v <- sqrt(diag(cov.m))
  cor.m <- cov.m / outer(sd.v,sd.v)
  #diag(cor.m) <- 1         # Set diagonal elements to 1
  eigen <- eigen(cor.m)
  eigen$values[eigen$values<.Machine$double.eps] <- 0 # handle numerical error
  sd.m.2 <- crossprod(t(eigen$vectors), diag(sqrt(eigen$values)))
  X <- crossprod(t(sd.m.2), matrix(stats::rnorm(J*n.sim),ncol=n.sim))
  MaxX <- apply(abs(X),2,max)
  stats::quantile(MaxX,prob)
}

# eigen.fd2 <- function(cov.fd){
#   W2inv <- inprod(cov.fd$sbasis,cov.fd$tbasis) ## sbisis and tbasis should be the same.
#   coefs <- cov.fd$coefs
#   A <- coefs %*% W2inv
#   e.A <- eigen(A)
#   e.A$values <- Re(e.A$values)
#   e.A$vectors <- Re(e.A$vectors)
#   cut <- sum(e.A$values > .Machine$double.eps)
#   vcoefs <- e.A$vectors[,1:cut]
#   normv <- diag(t(vcoefs) %*% W2inv %*% vcoefs)
#   vcoefs <- vcoefs %*% diag(1/sqrt(normv))
#   harmonics <- fd(vcoefs,cov.fd$sbasis)
#   rtnobj <- list(values=e.A$values,harmonics=harmonics)
#   class(rtnobj) <- "eigen.fd"
#   return(rtnobj)
# }