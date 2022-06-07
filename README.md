
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R-package: ffscb

<!-- badges: start -->

[![R-CMD-check](https://github.com/lidom/ffscb/workflows/R-CMD-check/badge.svg)](https://github.com/lidom/ffscb/actions)

<!-- badges: end -->

## Description

Fast ‘n’ fair simultaneous confidence bands for functional parameters as
introduced in the paper [**Fast and Fair Simultaneous Confidence Bands
for Functional Parameters**
(arXiv:1910.00131)](https://arxiv.org/abs/1910.00131) by [Dominik
Liebl](www.dliebl.com) and [Matthew
Reimherr](http://www.personal.psu.edu/mlr36/).

The folders `Simulations` and `Applications` contain the `r-scripts` for
reproducing the simulations and real-data applications of the
manuscript.

## Installation

``` r
devtools::install_github("lidom/ffscb")
```

## Small Example (Artifical Data)

``` r
library("ffscb")
# Generate a sample
p          <- 200 
N          <- 80 
grid       <- make_grid(p, rangevals=c(0,1))
mu         <- meanf_poly(grid,c(0,.25)) 
names(mu)  <- grid
cov.m      <- make_cov_m(cov.f = covf_nonst_matern, grid=grid, cov.f.params=c(2, 1/4, 1/4))
sample     <- make_sample(mu,cov.m,N)

# Compute the estimate, hat.mu, and its covariance, hat.cov.mu
hat.mu     <- rowMeans(sample)
hat.cov    <- crossprod(t(sample - hat.mu)) / N
hat.cov.mu <- hat.cov / N

# Compute the tau-parameter 
# I.e., the 'roughness parameter function' needed for the KR- and FFSCB-bands
hat.tau    <- tau_fun(sample)

# Make and plot confidence bands
b <- confidence_band(x          = hat.mu, 
                     cov.x      = hat.cov.mu, 
                     tau        = hat.tau, 
                     df         = N-1,
                     type       = c("FFSCB.t", "Bs", "BEc", "naive.t"),
                     conf.level = 0.95, 
                     n_int      = 4)
plot(b)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

<!-- Legend for band-`type`: -->
<!-- |Argument   | Description -->
<!-- |:-----------|:-------------- -->
<!-- |`FFSCB.t`  | Our fast 'n' fair simultaneous confidence band based on the $t$-distribution -->
<!-- |`FFSCB.z`  | Our fast 'n' fair simultaneous confidence band based on the normal distribution -->
<!-- |`Bs`       | Bootstrap based simultaneous confidence band similar to [Degras (2011)](http://www3.stat.sinica.edu.tw/statistica/j21n4/j21n412/j21n412.html) -->
<!-- |`BEc`      | Simultaneous confidence band of [Choi and Reimherr (2018)](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12239) -->
<!-- |`naive.t`  | Naive pointwise confidence band (not a simultaneous confidence band) -->

## Replicating the Sports Biomechanics Example

The following code replicates the sports biomechanics real data example
of our paper. First, we read in and setup the data.

``` r
library("ffscb")
data("Biomechanics")
## Data
slct_EC      <- grepl(pattern = "Extra_Cush_",  x = colnames(Biomechanics))
slct_NC      <- grepl(pattern = "Normal_Cush_", x = colnames(Biomechanics))
##
grid         <- Biomechanics[,1]/100   # Grid in [0,1]  
EC_mat       <- Biomechanics[,slct_EC] # Torque curves when running with extra  cushioned running shoes
NC_mat       <- Biomechanics[,slct_NC] # Torque curves when running with normal cushioned running shoes
Dff_mat      <- EC_mat - NC_mat        # Difference curves
N            <- ncol(Dff_mat)          # Sample size

## Plot
width     <- 7
height    <- 2.5
mar       <- c(4.1, 4.1, 4.1, 1.1)
##
layout(mat = matrix(c(1:3), nrow=1, ncol=3), widths = c(1, 1, 1))  
##
par(family = "serif", ps=13, cex.main=1.3, cex.lab=1.25, cex.axis=1.3, font.main = 1, mar=mar)
matplot(y = EC_mat,  x = grid*100, lwd=.5, col=gray(.25), type="l", lty=1, 
        ylab="N m / kg", xlab = "% of Stance Phase", main="Extra Cushioned\nRunning Shoe")
matplot(y = NC_mat,  x = grid*100, lwd=.5, col=gray(.25), type="l", lty=1, 
        ylab="N m / kg", xlab = "% of Stance Phase", main="Normal Cushioned\nRunning Shoe")
matplot(y = Dff_mat, x = grid*100, lwd=.5, col=gray(.5), type="l", lty=1, ylim = c(-0.25, max(Dff_mat)), 
        ylab="N m / kg", xlab = "% of Stance Phase", main="Difference Curves\n(Extra - Normal)")
lines(y = rowMeans(Dff_mat), x = grid*100)
legend("bottomright", legend = expression(paste("Estimated mean")), col="black", lty=1, bty="n", lwd = 1.5, cex=1.3)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

Next, we need to compute the estimate of the mean function `hat.mu` and
the estimate of its covariance function `hat.cov.mu`. From the latter,
we can compute the estimate of the roughness parameter `hat.tau`.

``` r
## Computing the estimates
hat_mu       <- rowMeans(Dff_mat)
hat.cov      <- crossprod(t(Dff_mat - hat_mu)) / (N-1)
hat.cov.mu   <- hat.cov / N
hat.tau      <- tau_fun(Dff_mat) 
```

Now, we have everything in place to compute the fair simultaneous 95%
confidence band that balances the false positive rate *α* = 0.05 by
allocating the equal shares *α*(1/`n_int`) to each of the here `n_int=8`
sub-intervalls \[0,12.5%\], \[12.5%,25\], …, \[87.5%,100%\].

Note: The paper treats also the case of balanced false positive rates
over arbitrary partitions of the domain, but this feature is not yet
implemented in our package.

``` r
alpha.level  <- 0.05
n_int        <- 8

## Over each subinterval [0, 12.5%], [12.5%, 25%], ... , [87.5%, 100%]  
## we compute a 99.375% simultanouse confidence interval:
(1 - alpha.level*(1/n_int)) * 100 
#> [1] 99.375

## Computing the confidence bands
b                 <- confidence_band(x          = hat_mu, 
                                     cov        = hat.cov.mu, 
                                     tau        = hat.tau, 
                                     df         = N-1, 
                                     type       ="FFSCB.t", 
                                     conf.level = 1-alpha.level, 
                                     n_int      = n_int)
```

A quick and dirty plot can be generated by `plot(b)`. The following
code, however, replicates the plot in our manuscript.

``` r
FF8_t_band        <- b[,-1] 
FF8_crit_value_u  <- (FF8_t_band[,1] - hat_mu)/sqrt(diag(hat.cov.mu))

## Plots
width     <- 7
height    <- 3.125
cex       <- .9
cexs      <- 0.95
##
layout(mat = matrix(c(1,2,3,3), nrow=2, ncol=2),
       heights = c(1, 1),     # Heights of the two rows
       widths =  c(1, 1.25))  # Widths of the two columns

par(family = "serif", ps=13, cex.axis=1.05, font.main = 1)
par(mar=c(3.1, 2.1, 2, 1.1))
plot(y=hat.tau,x=grid*100, type="l", main = "", xlab = "", ylab="")
mtext(text = expression(paste("Roughness estimate ",hat(tau))), 3, line = 0.4, adj = 0, cex=cex)
##
par(mar=c(3.1, 2.1, 2, 1.1))
matplot(y=cbind(FF8_crit_value_u), x=grid*100, col="black", lty=c(2,1), type = "l", 
        main="", xlab = "", ylab="", lwd=c(1.5,1), ylim=c(min(FF8_crit_value_u),4.1))
abline(v=c( 0:8 * 100/8), col="gray", lwd=1.3)
mtext(text = expression(paste("Fair adaptive critical value function ", hat(u)[alpha/2]^"*")), 3, line = 0.4, adj = 0, cex=cex)
mtext(text = "% of Stance Phase", 1, line = 2.25, cex=cexs)
## 
par(mar=c(3.1, 4.1, 3.1, 2.1))
matplot( y = b, x = grid*100, type="n", ylim=c(min(FF8_t_band),0.28), 
         ylab="", xlab = "", main="")
polygon(x=c(grid*100,rev(grid*100)), y=c(FF8_t_band[,2],rev(FF8_t_band[,1])), col = gray(.75), border = gray(.75))
abline(  h = 0, lwd=0.7)
abline(v=c( 0:8 * 100/8), col="gray")
lines(   y = hat_mu,  x = grid*100, col=1, lty=1)
axis(4, at = 0, labels = expression(H[0]:~theta==0))
legend(x=50, y=0.3, legend = c(expression(paste("Estimated mean")), 
                               expression(FF[t]^8)), y.intersp=1.05, bg="white", box.col = "white",
       lty=c(1,1), lwd = c(1.5,10), col=c("black", gray(.75)), cex =cexs, seg.len=2)
mtext(text = "% of Stance Phase", 1, line = 2.25, cex=cexs)
mtext(text = "N m / kg", 2, line = 2.5, cex=cex)
mtext(text = "95% Simultaneous Confidence Band (SCB)\n99.375% SCB over each subinterval", 3, line = 0.4, cex=cex)
text(x = 12.5/2, y = .2, labels="99.375% SCB", srt=90)
box()
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />
