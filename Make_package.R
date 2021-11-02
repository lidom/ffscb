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



## Generate a sample
p          <- 200 
N          <- 80 
grid       <- make_grid(p, rangevals=c(0,1))
mu0        <- meanf_poly(grid,c(0,0))   
names(mu0) <- grid
mu         <- meanf_poly(grid,c(0,0)) 
names(mu)  <- grid
cov.m      <- make_cov_m(cov.f = covf.nonst.matern, grid=grid, cov.f.params=c(2, 1/4, 1/4))
sample     <- make_sample(mu,cov.m,N)

                                        # Compute the estimate and its covariance
hat.mu     <- rowMeans(sample)
hat.cov    <- crossprod(t(sample - hat.mu)) / N
hat.cov.mu <- hat.cov / N

                                        # Compute the tau-parameter (for the KR- and FFSCB-bands)
hat.tau    <- tau_fun(sample)

tau = hat.tau; diag.cov = diag(hat.cov.mu);  conf.level=0.95; n_int=4; df=N-1

x = hat.mu; cov.x = hat.cov.mu; tau = hat.tau; df=N-1; conf.level  = 0.95; 


FFSCB_t1 <- make_band_FFSCB_t(x = hat.mu, diag.cov.x = diag(hat.cov.mu), tau = hat.tau, df = N-1, conf.level = 0.95, n_int = 1)
FFSCB_t2 <- make_band_FFSCB_t(x = hat.mu, diag.cov.x = diag(hat.cov.mu), tau = hat.tau, df = N-1, conf.level = 0.95, n_int = 8)
FFSCB_z <- make_band_FFSCB_z(x = hat.mu, diag.cov.x = diag(hat.cov.mu), tau = hat.tau, conf.level = 0.95)

matplot(FFSCB_t1[,-1], type="l", lty=1, col=1)
matlines(FFSCB_t2[,-1], type="l", lty=1, col=2)
lines(y=hat.mu, x=1:p)
lines(y=mu0, x=1:p, col=3)

## Make and plot confidence bands
b <- confidence_band(x=hat.mu, cov.x=hat.cov.mu, tau=hat.tau, df=N-1,
                     type=c("FFSCB.z", "FFSCB.t", "Bs","BEc","naive.t")[c(1,2,5)],
                     conf.level  = 0.95)
plot(b)
lines(x=grid, y=mu0, lty=2)

yy <- c(0, 2., 4., 3.)
xx <- seq(0,1,len=length(y0))

sfun0  <- approxfun(xx, yy, method = "const")

sfun0(1/3)

plot(sfun0)#,xlim = c(0,1))
