## ################################
##
## Make Package Codes
##
## ################################

## Remove pkg 
remove.packages("ffscb")

## Create/update documentation and (re-)write NAMESPACE
devtools::document("ffscb")

## CRAN-check pkg
# devtools::check("ffscb")       # check the package

## Install
devtools::install_local("ffscb", force = TRUE)
##
library("ffscb")
help("ffscb")
## #################################



                                        # Generate a sample
p          <- 200 
N          <- 80 
grid       <- make_grid(p, rangevals=c(0,1))
mu0        <- meanf_poly(grid,c(0,1))   
names(mu0) <- grid
mu         <- meanf_poly(grid,c(0,1.1)) 
names(mu)  <- grid
cov.m      <- make_cov_m(cov.f = covf.nonst.matern, grid=grid, cov.f.params=c(2, 1/4, 1/4))
sample     <- make_sample(mu,cov.m,N)

                                        # Compute the estimate and its covariance
hat.mu     <- rowMeans(sample)
hat.cov    <- crossprod(t(sample - hat.mu)) / N
hat.cov.mu <- hat.cov / N

                                        # Compute the tau-parameter (for the KR- and FFSCB-bands)
hat.tau    <- tau_fun(sample)
## Make and plot confidence bands
b <- confidence_band(x=hat.mu, cov.x=hat.cov.mu, tau=hat.tau, df=N-1,
                     type=c("FFSCB.t", "Bs","BEc","naive.t"),
                     conf.level  = 0.95)
plot(b)
lines(x=grid, y=mu0, lty=2)


