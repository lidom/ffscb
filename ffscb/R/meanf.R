#' @export
meanf.poly <- function(x,params=c(0,1)){ f <- params[1] + params[2]*(10*x^3 - 15*x^4 + 6*x^5) ; names(f) <- x ; return(f)} #params=c(shift,scale)

#' @export
meanf.peak <- function(x,params=c(0,1,2,16,1)){
  peak <- params[3] - params[4]*abs(x-0.5)^(params[5])
  peak[peak<0] <- 0   #params=c(shift,scale,peak-top,peak-sharpness,peak-power )
  f <- params[2] * (params[1] + peak)
  names(f) <- x
  return(f)
}

#' @export
meanf.scale <- function(x, delta=0){ meanf.poly(x,c(0,1+delta)) }

#' @export
meanf.shift <- function(x, delta=0){ meanf.poly(x,c(delta,1)) }

#' @export
meanf.localshift <- function(x, delta=0){ meanf.peak(x, c(0,1,1+delta,10,1)) }
