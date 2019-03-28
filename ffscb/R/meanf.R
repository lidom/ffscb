#' Meanfunction (polynomial)
#'
#' @param x function argument
#' @param params Parameters: params=c(shift,scale). 
#' @example 
#' curve(meanf.poly(x, c(0,2)), from=0, to=1, 
#' main="Meanfct Poly", ylab="",xlab="")
#' curve(meanf.poly(x, c(0,1)), from=0, to=1, 
#' lty=2, add=TRUE)
#' @export
meanf.poly <- function(x,params=c(0,1)){ f <- params[1] + params[2]*(10*x^3 - 15*x^4 + 6*x^5) ; names(f) <- x ; return(f)} #params=c(shift,scale)



#' Meanfunction with local rectangle
#'
#' @param x function argument
#' @param params parameters params=c(start, end, height) 
#' @example 
#' curve(meanf.rect(x, c(0, 1/4, 1/5)), from=0, to=1,
#' main="Meanfct Rect", ylab="",xlab="")
#' curve(meanf.rect(x, c(2/4, 3/4, 1/5)), from=0, to=1, 
#' lty=2, add=TRUE)
#' @export
meanf.rect <- function(x, params=c(0, 1/4, 1/5)){
  tmp <- numeric(length(x))
  tmp[params[1]<=x&x<=params[2]] <- params[3] 
  return(tmp)
  }


#' Meanfunction with local peak
#'
#' @param x function argument
#' @param params Parameters
#' @example 
#' curve(meanf.peak(x, c(0,1,2,16,1)), from=0, to=1, 
#' main="Meanfct Peak", ylab="",xlab="")
#' curve(meanf.peak(x, c(-1,1,2,16,1)), from=0, to=1, 
#' lty=2, add=TRUE)
#' curve(meanf.peak(x, c(0,2,2,16,1)), from=0, to=1, 
#' lty=2, add=TRUE)
#' curve(meanf.peak(x, c(0,1,1+1,10,1)), from=0, to=1, 
#' lty=2, add=TRUE)
#' @export
meanf.peak <- function(x,params=c(0,1,2,16,1)){
  peak <- params[3] - params[4]*abs(x-0.5)^(params[5])
  peak[peak<0] <- 0   #params=c(shift,scale,peak-top,peak-sharpness,peak-power )
  f <- params[2] * (params[1] + peak)
  names(f) <- x
  return(f)
}


#' Meanfunction (polynomial simple scaling)
#'
#' @param x function argument
#' @param delta scaling parameter 
#' @example 
#' curve(meanf.scale(x, 1), from=0, to=1)
#' curve(meanf.scale(x, 0), from=0, to=1, lty=2, add=TRUE)
#' @export
meanf.scale <- function(x, delta=0){ meanf.poly(x,c(0,1+delta)) }


#' Meanfunction (polynomial simple shifting)
#'
#' @param x function argument
#' @param delta shifting parameter 
#' @example 
#' curve(meanf.shift(x, 0), from=0, to=1)
#' curve(meanf.shift(x,.1), from=0, to=1, lty=2, add=TRUE)
#' @export
meanf.shift <- function(x, delta=0){ meanf.poly(x,c(delta,1)) }


#' Meanfunction with local peak (imple shifting)
#'
#' @param x function argument
#' @param delta shifting parameter 
#' @example 
#' curve(meanf.localshift(x, 0), from=0, to=1)
#' curve(meanf.localshift(x,.1), from=0, to=1, lty=2, add=TRUE)
#' @export
meanf.localshift <- function(x, delta=0){meanf.peak(x, c(0,1,1+delta,10,1)) }
