# Supplementary functions of other functions.

#' @export
sf.f1 <- function(z){ z * exp(z^2/2) * pnorm2(z) } ###### Make f1 function
#' @export
sf.difff1 <- function(y,x){ (y-sf.f1(x))^2 } ###### Make f1-inv using optimizaiton (not using table)
#' @export
sf.f1.inv <- function(y){ optim(2,fn=sf.difff1,y=y,method="Brent",lower=0,upper=30)$par} # initial 2, Min 0, Max 30
#' @export
sf.diffP <- function(m,lambdas,prob){ ### define the function to minimize #lambdas ;eigen values #prob : target probability
  tmpZ <- sapply(m / (sqrt(2*pi) * lambdas), sf.f1.inv)
  abs(prod(pnorm2(tmpZ)) - prob)
}
