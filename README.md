# R-package: ffscb
The R-package `ffscb` contains the fast 'n' fair simultaneous confidence bands for functional parameters as introduced in the paper [**Fast and Fair Simultaneous Confidence Bands for Functional Parameters**](http://arxiv.org/abs/1910.00131) by [Dominik Liebl](www.dliebl.com) and [Matthew Reimherr](http://www.personal.psu.edu/mlr36/).

## Installation 
```r
devtools::install_github("lidom/ffscb")
```

## Getting Started
```r
library("ffscb")
help("ffscb")
```

## Example Codes

```r
# Generate a sample
p          <- 200 
N          <- 80 
grid       <- make_grid(p, rangevals=c(0,1))
mu0        <- meanf_poly(grid,c(0,1))   
names(mu0) <- grid
mu         <- meanf_poly(grid,c(0,1.1)) 
names(mu)  <- grid
cov.m      <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2/2,1))
sample     <- make_sample(mu,cov.m,N)

# Compute the estimate, hat.mu, and its covariance, hat.cov.mu
hat.mu     <- rowMeans(sample)
hat.cov    <- crossprod(t(sample - hat.mu)) / N
hat.cov.mu <- hat.cov / N

# Compute the tau-parameter 
# I.e., the 'roughness parameter function' needed for the KR- and FFSCB-bands
hat.tau    <- tau_fun(sample)

# Make and plot confidence bands
b <- confidence_band(x=hat.mu, cov.x=hat.cov.mu, tau=hat.tau, df=N-1,
                     type=c("FFSCB.t", "Bs","BEc","naive.t"),
                     conf.level  = 0.95)
plot(b)
lines(x=grid, y=mu0, lty=2)
```

