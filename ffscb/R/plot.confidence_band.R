
#' Visualizes confidence bands constructed from 
#' \link{confidence_band}.
#'
#' @param x A 'confidence_band' object, the output from 
#' \link{confidence_band} funciton.
#' @param center Whether to include the functional estimate or not.
#' @param legendx position `x' of the legend. If NULL is passed, the 
#' legend will not be drawn (However, it may be added manually)
#' @param legendy position `y' of the legend.
#' @param ... Graphical parameters to be passed/overrided. If 
#' 'center' is TRUE, the first elements of 'col', 'lwd', 'lty' will
#'  be used for the estimate and the next ones will be used for the 
#'  bands, but using the same values for one pair, i.e., lower and 
#'  upper bounds.
#' @method plot confidence_band
#' @export
plot.confidence_band <- function(
  x, center=TRUE, legendx="topleft", legendy=NULL, ...
){
  ##
  band <- x
  # if(is.null(rownames(band))){rownames(band) <- seq(from=0,to=1,len=nrow(band))}
  
  ##
  if (is.null(dim(band))) {
    type <- "fd"
    class(band) <- "fd"
  } else {
    type <- "vector"
  }
  
  if (type=="vector") { 
    bandnames <- colnames(band)
    rownames(band) <- seq(from=0,to=1,len=nrow(band))
  } else {
    bandnames <- colnames(band$coefs)
  }
  
  bandnames <- sub(".u","",bandnames)
  bandnames <- sub(".l","",bandnames)
  bandnames <- bandnames[c(T,F)]
  
  ## extract graphical parameters from '...'
  gp       <- list(...)
  gpi_col  <- which(names(gp)=="col")
  gpi_lwd  <- which(names(gp)=="lwd")
  gpi_lty  <- which(names(gp)=="lty")
  gpi_xlab <- which(names(gp)=="xlab")
  gpi_ylab <- which(names(gp)=="ylab")
  
  if (length(gpi_col )==1) {
    col  <- gp[[gpi_col]]  
  } else { 
    col  <- c('black', '#e41a1c','#377eb8','#4daf4a',
              '#984ea3','#ff7f00','#ffff33','#a65628')
  }
  if (length(gpi_lty )==1) {
    lty  <- gp[[gpi_lty]]  
  } else {
    lty  <- c(1,5,3,4,2,6)
  }
  if (length(gpi_lwd )==1) {
    lwd  <- gp[[gpi_lwd]]  
  } else {
    lwd  <- c(1,1,1,1,1,1)
  }
  if (length(gpi_xlab)==1) {
    xlab <- gp[[gpi_xlab]] 
  } else { 
    xlab <- "T"
  }
  if (length(gpi_ylab)==1) {
    ylab <- gp[[gpi_ylab]] 
  } else {
    ylab <- "bands"
  }
  
  # reserve first value for the center (if center==TRUE)
  if (center) {
    colb <- c(col[-1],col[1])
    ltyb <- c(lty[-1],lty[1])
    lwdb <- c(lwd[-1],lwd[1])
  } else {
    colb <- col              
    ltyb <- lty               
    lwdb <- lwd    
  }
  
  gp$col  <- rep(colb,each=2)
  gp$lwd <- rep(lwdb,each=2)
  gp$lty <- rep(ltyb,each=2)
  gp$xlab <- xlab
  gp$ylab <- ylab
  
  if (type=="vector") {
    gp$x <- rownames(band)
    gp$y <- band[,-1]
    gp$type <- "l"
    do.call(graphics::matplot, gp)
  } else {
    gp$x <- band[-1]
    do.call(fda::plot.fd, gp)
  }
  
  if (center) {
    if (type=="vector") {
      graphics::lines(
        rownames(band),band[,1],col=col[1],lwd=lwd[1],lty=lty[1]
      )
    } else {
      fda::plot.fd(band,col=col,lwd=lwd[1],lty=lty[1],add=TRUE)
    }
    
    if (!is.null(legendx)) {
      graphics::legend(
        x=legendx,y=legendy,legend=c("Estimate",bandnames[-1]),
        col=col,lty=lty,lwd=lwd
      )
    } 
    
  } else {
    if (!is.null(legendx)) {
      graphics::legend(
        x=legendx,y=legendy,legend=bandnames[-1],
        col=col,lty=lty,lwd=lwd
      )
    }
  }
}
