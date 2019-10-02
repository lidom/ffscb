# R-package: ffscb
The R-package `ffscb` contains the fast 'n' fair simultaneous confidence bands for functional parameters as introduced in the paper [**Fast and Fair Simultaneous Confidence Bands for Functional Parameters** (arXiv:1910.00131)](http://arxiv.org/abs/1910.00131) by [Dominik Liebl](www.dliebl.com) and [Matthew Reimherr](http://www.personal.psu.edu/mlr36/).

The folder `Applications` contains the `r-scripts` to reproduce the applications in the manuscript.

## Installation 
```r
devtools::install_github("lidom/ffscb/ffscb")
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
mu         <- meanf_poly(grid,c(0,.25)) 
names(mu)  <- grid
cov.m      <- make_cov_m(cov.f = covf.nonst.matern, grid=grid, cov.f.params=c(2, 1/4, 1/4))
sample     <- make_sample(mu,cov.m,N)

# Compute the estimate, hat.mu, and its covariance, hat.cov.mu
hat.mu     <- rowMeans(sample)
hat.cov    <- crossprod(t(sample - hat.mu)) / N
hat.cov.mu <- hat.cov / N

# Compute the tau-parameter 
# I.e., the 'roughness parameter function' needed for the KR- and FFSCB-bands
hat.tau    <- tau_fun(sample)

# Make and plot confidence bands
b <- confidence_band(x = hat.mu, cov.x = hat.cov.mu, tau = hat.tau, df = N-1,
                     n_int = 5, t0 = 1, conf.level  = 0.95,
                     type=c("FFSCB.t", "Bs","BEc","naive.t"))
plot(b)
```

Legend for band-`type`:

|Argument   | Description
|:-----------|:--------------
|`FFSCB.t`  | Our fast 'n' fair simultaneous confidence band
|`Bs`       | Bootstrap based simultaneous confidence band similar to [Degras (2011)](http://www3.stat.sinica.edu.tw/statistica/j21n4/j21n412/j21n412.html)
|`BEc`      | Simultaneous confidence band of [Choi and Reimherr (2018)](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12239)
|`naive.t`  | Naive pointwise confidence band