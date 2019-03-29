# Load packages 
library("devtools")
library("fregion")
#devtools::install_github("alessiapini/fdatest")
library("fdatest")
library("fda")

p         <- 101 
N         <- c(10,50,100)[1]
rangeval  <- c(0,1)
grid      <- make_grid(p, rangevals=rangeval)#, type == "close")
mu        <- meanf_poly(grid, params = c(0,0)) # plot(x=grid,y=mu)
cov.m     <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2, 1, 1))
t0        <- grid[1]


##
x         <-  make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
matplot(grid, x, type="l", lty=1); lines(grid, mu, lwd=2); confint(lm(x[1,]~1)) # naive.t



IWT_res  <- fdatest::IWT1(data = t(x), mu = 0)
plot(IWT_res)
par(mfrow=c(2,1))
plot(rowMeans(x),type="l");abline(h=0,lty=3)
plot(IWT_res$adjusted_pval,type="l")

