library("ffscb")
# library("tidyverse")
# 
# ## ################################################
# load("Simulation_Results/DGP3_N=10")
# 
# type       <- c("naive.t", "Bs", "BEc", "KR.t", "FFSCB.t")
# type_order <- pmatch(type, names(widths))
# type       <- type[type_order]
# crossings_df <- dplyr::tibble(
#   !!type[1] :=  c(crossings_loc[,,1]),
#   !!type[2] :=  c(crossings_loc[,,2]),
#   !!type[3] :=  c(crossings_loc[,,3]),
#   !!type[4] :=  c(crossings_loc[,,4]),
#   !!type[5] :=  c(crossings_loc[,,5])
# ) %>% 
#   tidyr::gather(band, cr_loc, factor_key=TRUE) %>% 
#   dplyr::mutate(intervals = as.factor(findInterval(x   = cr_loc, 
#                                                    vec = seq(0,1,len=5), 
#                                                    rightmost.closed = T)))
# 
# 
# crossings_df %>% 
#   tidyr::drop_na() %>% 
#   dplyr::group_by(band) %>% 
#   dplyr::mutate(n_cr = n()) %>% 
#   dplyr::group_by(band, intervals) %>% 
#   dplyr::summarise(rfrq = n()/n_cr[1]) %>% 
#   mutate_if(is.numeric, round, 2) %>% 
#   print(n=Inf)
# ## ################################################



set.seed(1110)
p           <- 200 
N           <- 15
rangeval    <- c(0,1)
grid        <- make.grid(p, rangevals=rangeval)
mu          <- meanf.poly(grid, params = c(0,0)) # plot(x=grid,y=mu)
names(mu)   <- grid
cov.m       <- make.cov.m(cov.f = covf.st.matern.warp.power, 
                          grid=grid, cov.f.params=c(1.25, 1, 1, 2.5))
n_int       <- 10
type        <- c("naive.t", "Bs", "FFSCB.t")
alpha.level <- 0.10
dat         <-  make.sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
## Estimate mean, covariance, and tau
hat.mu      <- rowMeans(dat)
hat.cov.m   <- crossprod(t(dat - hat.mu)) / (N-1)
hat.tau.v   <- tau_fun(dat)# plot(hat.tau.v)
## Confidence bands
b           <- fregion.band(x=hat.mu,cov=hat.cov.m,tau=hat.tau.v,N=N,type=type,conf.level=(1-alpha.level),n_int=n_int)# plot(b)
##

## PLOT
pdf(file = "Fig1.pdf", width = 9, height = 4.5)
par(mfrow=c(1,2), mar=c(3.1, 3.1, 2.1, 1.1), family="serif")
matplot(y=dat, x=grid, axes=F, xlab="", ylab="", type = "l", lty=1, col=gray(.5,.5))
mtext("Functions with Inhomogeneous Roughness", side = 3, line = .5, adj = 0)
axis(1, at=seq(0,1,len=5)); axis(2); box()
lines(y=hat.mu, x=grid,   lty=4, lwd=1.5)
lines(y=c(0,0), x=c(0,1), lty=1, lwd=1.5)
legend("topleft", legend = c("Estimated Mean", "True Mean"), lty=c(4,1), bty="n", lwd = 1.3)
##
matplot(y=b[,-1], x=grid, type = "n", axes=F, xlab="", ylab="", xlim=c(-.025,1), ylim=c(-1.25,1.5))
mtext("Simultaneous 90% Confidence Bands", side = 3, line = .5, adj = 0)
axis(1, at=seq(0,1,len=5)); axis(2); box()
matlines(y=b[,-1], x=grid, lty=rep(c(2,3,1),each=2), col="1")
lines(y=hat.mu, x=grid, lty=4, lwd=1.5)
legend(x = par("usr")[1], y = par("usr")[4], 
       legend = c("FFSCB",
                  "Bootstrap",
                  "naive t"),lty = c(1,2,3), bty="n")
legend(x = par("usr")[1]+.25, y = par("usr")[4], 
       legend = c("(emp. coverage: 0.90)", 
                  "(emp. coverage: 0.91)",
                  "(emp. coverage: 0.72)"), lty = c(1), col=gray(.5,0), bty="n")
textloc1 <- -0.6
textloc2 <- -1
textloc3 <- -.8
text(x = par("usr")[1], y= textloc1, labels = "Intervalwise exceedances (rel. frequencies):", pos = 4, cex = .9)
#lines(x=c(0,0), y=c(par("usr")[3],-1.25))
#text(x = -0.1,    y=textloc2, labels = c("FFSCG\nBootstr\nnaive t"), cex = .9)
lines(x = c(par("usr")[1]+.05, par("usr")[1]+.125), y=rep(-0.85,2), lty=1, lwd=1.1)
lines(x = c(par("usr")[1]+.05, par("usr")[1]+.125), y=rep(-1.00,2), lty=2, lwd=1.1)
lines(x = c(par("usr")[1]+.05, par("usr")[1]+.125), y=rep(-1.15,2), lty=3, lwd=1.1)
text(x = .125*1,    y=textloc2, labels = c("0.24\n0.06\n0.02"), cex = .9)
lines(x=c(.25,.25), y=c(par("usr")[3],textloc3))
text(x = .125*3,    y=textloc2, labels = c("0.25\n0.19\n0.22"), cex = .9)
lines(x=c(.5,.5),   y=c(par("usr")[3],textloc3))
text(x = .125*5,    y=textloc2, labels = c("0.25\n0.31\n0.35"), cex = .9)
lines(x=c(.75,.75), y=c(par("usr")[3],textloc3))
text(x = .125*7,    y=textloc2, labels = c("0.26\n0.44\n0.42"), cex = .9)
#lines(x=c(1,1), y=c(par("usr")[3],-1.25))
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

