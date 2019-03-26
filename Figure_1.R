library("fregion")
set.seed(1110)
p           <- 200 
N           <- 15
rangeval    <- c(0,1)
grid        <- make.grid(p, rangevals=rangeval)
mu          <- meanf.poly(grid, params = c(0,0)) # plot(x=grid,y=mu)
names(mu)   <- grid
cov.m       <-  make.cov.m(cov.f = covf.st.matern.warp.power, grid=grid, cov.f.params=c(1.25, 1, 1, 2))
n_int       <- 4
type        <- c("naive.t", "KR.t", "ESCB.t")
alpha.level <- 0.10
dat         <-  make.sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
## Estimate mean, covariance, and tau
hat.mu      <- rowMeans(dat)
hat.cov.m   <- crossprod(t(dat - hat.mu)) / (N-1)
hat.tau.v   <- tau_fun(dat)# plot(hat.tau.v)
## Confidence bands
b           <- fregion.band(x=hat.mu,cov=hat.cov.m,tau=hat.tau.v,N=N,type=type,conf.level=(1-alpha.level),n_int=n_int)# plot(b)
##
par(mfrow=c(1,2), mar=c(3.1, 3.1, 2.1, 1.1), family="serif")
matplot(y=dat, x=grid, axes=F, xlab="", ylab="", type = "l", lty=1, col=gray(.5,.5))
axis(1, at=seq(0,1,len=(n_int+1))); axis(2); box()
lines(y=hat.mu, x=grid, lty=4, lwd=1.5)
legend("topleft", legend = c("Estimated Mean", "True Mean = 0"), lty=c(4,NA), bty="n", lwd = 1.5)
##
matplot(y=b[,-1], x=grid, type = "n", axes=F, xlab="", ylab="", xlim=c(0,1), ylim=c(-1.75,max(b[,-1])))
axis(1, at=seq(0,1,len=(n_int+1))); axis(2); box()
matlines(y=b[,-1], x=grid, lty=rep(length(type):1,each=2), col="1")
lines(y=hat.mu, x=grid, lty=4, lwd=1.5)
legend("topleft", legend = c("Fair SCB", "Kac-Rice SCB", "Pointwise t"),lty = c(1,2,3), bty="n")
text(x = 0, y= -1, labels = "Intervalwise exceedances:", pos = 4)
lines(x=c(0,0), y=c(par("usr")[3],-1.25))
text(x = .125*1, y=-1.5, labels = c("0.00\n0.00\n0.00"))
lines(x=c(.25,.25), y=c(par("usr")[3],-1.25))
text(x = .125*3, y=-1.5, labels = c("0.00\n0.00\n0.00"))
lines(x=c(.5,.5), y=c(par("usr")[3],-1.25))
text(x = .125*5, y=-1.5, labels = c("0.00\n0.00\n0.00"))
lines(x=c(.75,.75), y=c(par("usr")[3],-1.25))
text(x = .125*7, y=-1.5, labels = c("0.00\n0.00\n0.00"))
lines(x=c(1,1), y=c(par("usr")[3],-1.25))
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))


