#' @export
pnorm2 <- function(z){1-2*(1-pnorm(z))} #### make two-sided version of pnorm with only positive z
#' @export
qnorm2 <- function(p){-qnorm((1-p)/2)}
#' @export
pt2 <- function(t,df){1-2*(1-pt(t,df))}
#' @export
qt2 <- function(p,df){-qt((1-p)/2,df)}
#' @export
make.grid <- function(p=100,rangevals=c(0,1)){seq(rangevals[1],rangevals[2],len=p)}
#' @export
make.cov.m <- function(cov.f=covf.st.matern, grid=100, cov.f.params=NULL){  ### Make cov. matrix from cov. function.
  if (length(grid)==1) {grid=make.grid(p=grid)} ## input grid as a single number (as grid size), or as vector (actual grid)
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
#' @export
make.sample <- function(mean.v,cov.m,N,dist="rnorm",...){
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
#' @export
get.req.n.pc <- function(proportion, lambdas){ # Finds how many pcs are required to achieve desired variance proportion(prop. as 1-10^(-x))
  satisfied  <- (proportion <= cumsum(lambdas) / sum(lambdas))
  min(c(1:length(lambdas))[satisfied])
}
#' @export
get.schisq.q.gamma <- function(weights,prob){  # prob as 1-alpha
  meanvalue  <- sum(weights)
  varvalue   <- sum(2*weights^2)
  gammascale <- varvalue / meanvalue   # 's'
  gammashape <- meanvalue / gammascale # 'a'
  qgamma(p=prob,shape=gammashape,scale=gammascale,lower.tail=TRUE)
}

#' @export
eigen.fd2 <- function(cov.fd){
  W2inv <- inprod(cov.fd$sbasis,cov.fd$tbasis) ## sbisis and tbasis should be the same.
  coefs <- cov.fd$coefs
  A <- coefs %*% W2inv
  e.A <- eigen(A)
  e.A$values <- Re(e.A$values)
  e.A$vectors <- Re(e.A$vectors)
  cut <- sum(e.A$values > .Machine$double.eps)
  vcoefs <- e.A$vectors[,1:cut]
  normv <- diag(t(vcoefs) %*% W2inv %*% vcoefs)
  vcoefs <- vcoefs %*% diag(1/sqrt(normv))
  harmonics <- fd(vcoefs,cov.fd$sbasis)
  rtnobj <- list(values=e.A$values,harmonics=harmonics)
  class(rtnobj) <- "eigen.fd"
  return(rtnobj)
}


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


make.grid <- function(p=100,rangevals=c(0,1),type=c("close","open")[1]){
  if(!any(type==c("close","open"))){"'type' must be 'close' or 'open'."}
  if(type=="open"){
    seq(0.5/p, (p-0.5)/p, by=1/p) * (rangevals[2] - rangevals[1]) + rangevals[1]
  }
  if(type=="close"){
    seq(rangevals[1], rangevals[2], len=p)
  }
}

#' #' @export
#' counting_upcrosses <- function(x_vec, threshold) { 
#' ##
#' check        <- x_vec > threshold
#' names(check) <- 1:length(check)
#' ## Force all TRUE-values to be interior (for diff())
#' check1 <- c(FALSE,check,FALSE)
#' names(check1)[length(check1)] <- as.numeric(names(check1[length(check1)-1]))+1
#' check2 <- diff(check1)
#' check2 <- check2[check2!=0]
#' ##
#' if(sum(check2)!=0){stop("Error: Each upcrossing must have a start (1) and an end (-1)")}
#' ##
#' n_upcross <- length(check2)/2
#' ##
#' if(n_upcross>0){
#'   start_locations <- rep(NA,n_upcross)
#'   extents         <- rep(NA,n_upcross)
#'   for(i in 1:n_upcross){ # i <- 3
#'     # locate the upcrossing at its beginning
#'     start_locations[i] <- as.numeric(names(check2)[(2*(i-1)+1)])
#'     extents[i]         <- as.numeric(names(check2)[(2*(i-1)+2)]) - as.numeric(names(check2)[(2*(i-1)+1)])
#'   }
#' }else{
#'   start_locations <- NULL
#'   extents         <- 0
#' }
#' ##
#' return(list("n_upcross"       = n_upcross, 
#'             "start_locations" = start_locations,
#'             "extents"         = extents))
#' }
