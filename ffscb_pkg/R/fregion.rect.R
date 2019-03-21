#' Builds a Rectangular Confidence Region
#'
#' @param x Functional parameter estimate. It can be either a vector or \link{fd} object from \link{fda}.
#' @param cov N * Cov(X), in which X is the functional estimator. It can be either matrix or \link{bifd} object from \link{fda}. The eigen decomposition of Cov(X) can be used instead.
#' @param N It should be '1' if 'cov' is the covariance operator for X itself, which is the default value.
#' @param type This specifies which rectangular region to be constructed. It should be either one of "Rz", "Rz1", "Rzs", or "Rz1s".
#' @param conf.level Confidence level of the region.
#' @param pc.cut A numeric value. For integer values, fPC up to those values will be used. If it's a value from 0 to 1, this specifies the proportion of (estimated) variance that should be explained by the fPCs. If it is 0, all the available fPCs will be used as long as the size of eigenvalues are greater than .Machine$double.eps.
#' @param df Degree of freedom to use in small sample versions.
#' @return fregion.rect Use \link{plot.fregion.rect} to visualize.
#' @examples
#' # 1. Vector/matrix version
#'
#' # Generate a sample
#' p = 200 ; N = 80 ; rangeval = c(0,1)
#' grid = make.grid(p, rangevals=rangeval)
#' mu0 = meanf.poly(grid,c(0,1)) ; names(mu0) = grid
#' mu = meanf.poly(grid,c(0,1.1)) ; names(mu) = grid
#' cov.m = make.cov.m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2/2,1,1))
#' x = make.sample(mu,cov.m,N)
#'
#' # Find the estimate and covariance
#' hat.mu = rowMeans(x)
#' hat.cov.m = crossprod(t(x - hat.mu)) / (N-1)
#' e.hat.cov.m = eigen(hat.cov.m)   # <- This is optional and can be provide into the functions instead of hat.cov.m below.
#'
#' # Make rectangular region and visulize
#' c <- fregion.rect(hat.mu-mu0,hat.cov.m,N=N)
#' plot(c)
#'
#'
#' # 2. fd/bifd version
#'
#' # create basis, convert vector/matrix into fd/bifd objects.
#' require(fda)
#' nbasis <- round(p*0.9)
#' fd.basis <- create.bspline.basis(rangeval,nbasis)
#' mu0.fd <- Data2fd(names(mu0),mu0,fd.basis)
#' mu.fd <- Data2fd(names(mu),mu,fd.basis)
#' x.fd <- Data2fd(rownames(x),x,fd.basis)
#' hat.mu.fd <- mean.fd(x.fd)
#' hat.cov.fd <- var.fd(x.fd)
#' e.hat.cov.fd <- eigen.fd(hat.cov.fd)   # <- This is optional and can be provide into the functions instead of hat.cov.fd below.
#'
#' # Make rectangular region and visulize
#' c.fd <- fregion.rect(hat.mu.fd-mu0.fd,hat.cov.fd,N=N)
#' plot(c.fd)

#' @export

fregion.rect <- function(x, cov, N=1, type="Rz", conf.level=0.95, pc.cut=0.999, df=NULL){

  ### 1. Check the data type ###
  if (inherits(x,"fd") & (inherits(cov,"bifd") | inherits(cov,"pca.fd") | inherits(cov,"eigen.fd"))) datatype="fd" else if
  ((inherits(x,"numeric") | inherits(x,"matrix"))  & (inherits(cov,"matrix") | inherits(cov,"list") | inherits(cov,"eigen") )) datatype="vector" else stop ("The format of data is unknown")

  ### 2. If covariance is given as it is (not eigen decomposition), take eigen decomposition.
  e.cov <- cov
  if (inherits(cov,"bifd")) e.cov <- eigen.fd(cov)
  if (inherits(cov,"matrix")) e.cov <- eigen(cov)

  # Trim cov to have all non-negative eigenvalues (precautionary) and scale it according to N
  e.cov$values[e.cov$values<0] <- 0
  e.cov$values <- e.cov$values / N

  if (is.null(df)) df <- N - 1
  cut=pc.cut
  if (cut < 0) stop("fpc.cut should be some positive fraction from 0 to 1, or integer greater than or equal to 1, or just 0 (to use all available pcs)")
  if (cut == 0 ) cut <- sum(e.cov$values > .Machine$double.eps) else { # if fpc.cut is 0, use all PCs
    if (cut < 1) cut <- fregion::get.req.n.pc(cut,e.cov$values) else {
      if (cut != round(cut)) {cut=round(cut) ; print("fpc.cut[",i,"] was rounded to the closest integer")}
    }
  }
  pc.range <- c(1:cut)

  eigenvalues <- e.cov$values[pc.range]

  harmonics <- eigenfunctions <- NULL
  if (datatype=="fd"){
    harmonics <- e.cov$harmonics[pc.range]
  } else {
    eigenfunctions <- e.cov$vectors[,pc.range]
    rownames(eigenfunctions) <- names(x)
    colnames(eigenfunctions) <- paste0("v",pc.range)
  }

  if (datatype=="fd") {coef <- as.vector(inprod(x, harmonics))} else
                     {coef <- as.vector(crossprod(x, eigenfunctions))}
  score <- coef / sqrt(eigenvalues)
  if ("Rz" %in% type) {mult <- fregion::get.crit.Rz(eigenvalues,conf.level) ; marginal.coverage <- pnorm2(mult) }
  if ("Rz1" %in% type) {mult <- fregion::get.crit.Rz1(eigenvalues,conf.level) ; marginal.coverage <- pnorm2(mult) }
  if ("Rzs" %in% type) {mult <- fregion::get.crit.Rzs(eigenvalues,conf.level,df) ; marginal.coverage <- pt2(mult,df) }
  if ("Rz1s" %in% type) {mult <- fregion::get.crit.Rz1s(eigenvalues,conf.level,df) ; marginal.coverage <- pt2(mult,df) }
  ME <- mult * sqrt(eigenvalues)
  cum.var.explained <- cumsum(eigenvalues) / sum(e.cov$values[e.cov$values>0])
  var.explained <- eigenvalues / sum(e.cov$values[e.cov$values>0])
  table <- data.frame(PC=factor(pc.range), coefficient=coef, MoE=ME, score=score, abs.score=abs(score),
                     score.bound=mult, var.explained=var.explained, cum.var.explained=cum.var.explained, marginal.coverage=marginal.coverage)
  result <- list(table=table, eigenfunctions=eigenfunctions, harmonics=harmonics, type=type)
  class(result) <- "fregion.rect"
  return(result)
}
