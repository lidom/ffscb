## #########################################
## Plots of the Data Generating Processes
## #########################################
## Load packages 
library("ffscb")

p            <- 101
grid         <- make_grid(p, rangevals=c(0,1))
delta        <- 0
mu           <- meanf_shift(grid, delta)
##
cov.m_s      <- make_cov_m(cov.f = covf.st.matern,    grid=grid, cov.f.params=c(3/2, 1/4))
cov.m_r      <- make_cov_m(cov.f = covf.st.matern,    grid=grid, cov.f.params=c(1/2, 1/4))
cov.m_s2r    <- make_cov_m(cov.f = covf.nonst.matern, grid=grid, cov.f.params=c(4/2, 1/4, 1/4))


## Plots:
N            <- 100
n_int        <- 3
seq(0,1,len=n_int+1)
alpha.level  <- 0.05
type         <- c("Bs", "BEc", "naive.t", "FFSCB.t")
seed         <- 1
##
## Smooth
set.seed(seed)
sim.dat_s      <- make_sample(mean.v = mu, cov.m = cov.m_s,   N = N, dist = "rnorm")
t0_s           <- grid[1]
hat_mu_s       <- rowMeans(sim.dat_s)
hat.cov_s      <- crossprod(t(sim.dat_s - hat_mu_s)) / (N-1)
hat.cov.mu_s   <- hat.cov_s / N
hat.tau_s      <- tau_fun(sim.dat_s) 
b_s            <- confidence_band(x=hat_mu_s, cov=hat.cov.mu_s, tau=hat.tau_s, t0=t0_s, df=N-1, 
                                  type=type, conf.level=(1-alpha.level), n_int=n_int)
FFSCB_t_band_s <- b_s[, grepl("FFSCB.t", colnames(b_s))]
naive_t_band_s <- b_s[, grepl("naive.t", colnames(b_s))]
Bs_band_s      <- b_s[, grepl("Bs",      colnames(b_s))]
BEc_band_s     <- b_s[, grepl("BEc",     colnames(b_s))]
FFSCB_diff_s   <- FFSCB_t_band_s[,1] - FFSCB_t_band_s[,2]
naive_diff_s   <- naive_t_band_s[,1] - naive_t_band_s[,2]
Bs_diff_s      <- Bs_band_s[,1]      - Bs_band_s[,2]
BEc_diff_s     <- BEc_band_s[,1]     - BEc_band_s[,2]
## true cov
true.cov.mu_s  <- cov.m_s/N
true.cov.tau_s <- cov2tau_fun(cov.m_s)
true.cov.b_s   <- confidence_band(x=hat_mu_s, cov=true.cov.mu_s, tau=true.cov.tau_s, t0=t0_s, df=N-1, 
                                  type=type, conf.level=(1-alpha.level), n_int=n_int)
true.cov.FFSCB_t_band_s <- true.cov.b_s[, grepl("FFSCB.t", colnames(true.cov.b_s))]
true.cov.naive_t_band_s <- true.cov.b_s[, grepl("naive.t", colnames(true.cov.b_s))]
true.cov.Bs_band_s      <- true.cov.b_s[, grepl("Bs",      colnames(true.cov.b_s))]
true.cov.BEc_band_s     <- true.cov.b_s[, grepl("BEc",     colnames(true.cov.b_s))]
true.cov.FFSCB_diff_s   <- true.cov.FFSCB_t_band_s[,1] - true.cov.FFSCB_t_band_s[,2]
true.cov.naive_diff_s   <- true.cov.naive_t_band_s[,1] - true.cov.naive_t_band_s[,2]
true.cov.Bs_diff_s      <- true.cov.Bs_band_s[,1]      - true.cov.Bs_band_s[,2]
true.cov.BEc_diff_s     <- true.cov.BEc_band_s[,1]     - true.cov.BEc_band_s[,2]



##
## Rough  
set.seed(seed)
sim.dat_r      <- make_sample(mean.v = mu, cov.m = cov.m_r,   N = N, dist = "rnorm")
t0_r           <- grid[1]
hat_mu_r       <- rowMeans(sim.dat_r)
hat.cov_r      <- crossprod(t(sim.dat_r - hat_mu_r)) / (N-1)
hat.cov.mu_r   <- hat.cov_r / N
hat.tau_r      <- tau_fun(sim.dat_r) 
b_r            <- confidence_band(x=hat_mu_r, cov=hat.cov.mu_r, tau=hat.tau_r, t0=t0_r, df=N-1, 
                                  type=type, conf.level=(1-alpha.level), n_int=n_int)
FFSCB_t_band_r <- b_r[, grepl("FFSCB.t", colnames(b_r))]
naive_t_band_r <- b_r[, grepl("naive.t", colnames(b_r))]
Bs_band_r      <- b_r[, grepl("Bs",      colnames(b_r))]
BEc_band_r     <- b_r[, grepl("BEc",     colnames(b_r))]
FFSCB_diff_r   <- FFSCB_t_band_r[,1] - FFSCB_t_band_r[,2]
naive_diff_r   <- naive_t_band_r[,1] - naive_t_band_r[,2]
Bs_diff_r      <- Bs_band_r[,1]      - Bs_band_r[,2]
BEc_diff_r     <- BEc_band_r[,1]     - BEc_band_r[,2]
## true cov
true.cov.mu_r  <- cov.m_r/N
true.cov.tau_r <- cov2tau_fun(cov.m_r)
true.cov.b_r   <- confidence_band(x=hat_mu_r, cov=true.cov.mu_r, tau=true.cov.tau_r, t0=t0_r, df=N-1, 
                                  type=type, conf.level=(1-alpha.level), n_int=n_int)
true.cov.FFSCB_t_band_r <- true.cov.b_r[, grepl("FFSCB.t", colnames(true.cov.b_r))]
true.cov.naive_t_band_r <- true.cov.b_r[, grepl("naive.t", colnames(true.cov.b_r))]
true.cov.Bs_band_r      <- true.cov.b_r[, grepl("Bs",      colnames(true.cov.b_r))]
true.cov.BEc_band_r     <- true.cov.b_r[, grepl("BEc",     colnames(true.cov.b_r))]
true.cov.FFSCB_diff_r   <- true.cov.FFSCB_t_band_r[,1] - true.cov.FFSCB_t_band_r[,2]
true.cov.naive_diff_r   <- true.cov.naive_t_band_r[,1] - true.cov.naive_t_band_r[,2]
true.cov.Bs_diff_r      <- true.cov.Bs_band_r[,1]      - true.cov.Bs_band_r[,2]
true.cov.BEc_diff_r     <- true.cov.BEc_band_r[,1]     - true.cov.BEc_band_r[,2]
##
## From smooth to rough
set.seed(seed)
sim.dat_s2r    <- make_sample(mean.v = mu, cov.m = cov.m_s2r, N = N, dist = "rnorm")
t0_s2r         <- grid[1]
hat_mu_s2r     <- rowMeans(sim.dat_s2r)
hat.cov_s2r    <- crossprod(t(sim.dat_s2r - hat_mu_s2r)) / (N-1)
hat.cov.mu_s2r <- hat.cov_s2r / N
hat.tau_s2r    <- tau_fun(sim.dat_s2r) 
b_s2r          <- confidence_band(x=hat_mu_s2r, cov=hat.cov.mu_s2r, tau=hat.tau_s2r, t0=t0_s2r, df=N-1, 
                                  type=type, conf.level=(1-alpha.level), n_int=n_int)
FFSCB_t_band_s2r <- b_s2r[, grepl("FFSCB.t", colnames(b_s2r))]
naive_t_band_s2r <- b_s2r[, grepl("naive.t", colnames(b_s2r))]
Bs_band_s2r      <- b_s2r[, grepl("Bs",      colnames(b_s2r))]
BEc_band_s2r     <- b_s2r[, grepl("BEc",     colnames(b_s2r))]
FFSCB_diff_s2r   <- FFSCB_t_band_s2r[,1] - FFSCB_t_band_s2r[,2]
naive_diff_s2r   <- naive_t_band_s2r[,1] - naive_t_band_s2r[,2]
Bs_diff_s2r      <- Bs_band_s2r[,1]      - Bs_band_s2r[,2]
BEc_diff_s2r     <- BEc_band_s2r[,1]     - BEc_band_s2r[,2]
## true cov
true.cov.mu_s2r  <- cov.m_s2r/N
true.cov.tau_s2r <- cov2tau_fun(cov.m_s2r)
true.cov.b_s2r   <- confidence_band(x=hat_mu_s2r, cov=true.cov.mu_s2r, tau=true.cov.tau_s2r, t0=t0_s2r, df=N-1, 
                                    type=type, conf.level=(1-alpha.level), n_int=n_int)
true.cov.FFSCB_t_band_s2r <- true.cov.b_s2r[, grepl("FFSCB.t", colnames(true.cov.b_s2r))]
true.cov.naive_t_band_s2r <- true.cov.b_s2r[, grepl("naive.t", colnames(true.cov.b_s2r))]
true.cov.Bs_band_s2r      <- true.cov.b_s2r[, grepl("Bs",      colnames(true.cov.b_s2r))]
true.cov.BEc_band_s2r     <- true.cov.b_s2r[, grepl("BEc",     colnames(true.cov.b_s2r))]
true.cov.FFSCB_diff_s2r   <- true.cov.FFSCB_t_band_s2r[,1] - true.cov.FFSCB_t_band_s2r[,2]
true.cov.naive_diff_s2r   <- true.cov.naive_t_band_s2r[,1] - true.cov.naive_t_band_s2r[,2]
true.cov.Bs_diff_s2r      <- true.cov.Bs_band_s2r[,1]      - true.cov.Bs_band_s2r[,2]
true.cov.BEc_diff_s2r     <- true.cov.BEc_band_s2r[,1]     - true.cov.BEc_band_s2r[,2]


res_true_s2r <- make_band_FFSCB_t(x=hat_mu_s2r, diag.cov.x = diag(true.cov.mu_s2r), tau=true.cov.tau_s2r, t0=t0_s2r, df=N-1,
                                  conf.level=(1-alpha.level), n_int=n_int)

res_true_s2r$prob_t0 + res_true_s2r$a_star


res_est_s2r <- make_band_FFSCB_t(x=hat_mu_s2r, diag.cov.x = diag(hat.cov.mu_s2r), tau=hat.tau_s2r, t0=t0_s2r, df=N-1,
                                 conf.level=(1-alpha.level), n_int=n_int)

res_est_s2r$prob_t0 + res_est_s2r$a_star



res_true_s <- make_band_FFSCB_t(x=hat_mu_s, diag.cov.x = diag(true.cov.mu_s), tau=true.cov.tau_s, t0=t0_s, df=N-1,
                                conf.level=(1-alpha.level), n_int=n_int)

res_true_s$prob_t0 + res_true_s$a_star


res_true_r <- make_band_FFSCB_t(x=hat_mu_r, diag.cov.x = diag(true.cov.mu_r), tau=true.cov.tau_r, t0=t0_r, df=N-1,
                                conf.level=(1-alpha.level), n_int=n_int)

res_true_r$prob_t0 + res_true_r$a_star


# 
# x=hat_mu_s2r
# diag.cov= diag(hat.cov.mu_s2r)
# tau=hat.tau_s2r
# t0=t0_s2r
# df=N-1
# conf.level=(1-alpha.level)
# n_int=n_int


## Plots
width     <- 3.2
height    <- 6.2
mar_u1    <- c(0.5, 1.9, 2.5, 0.05)
mar_l1    <- c(3.1, 1.9, 0.1, 0.05)
mar_u2    <- c(0.5, 1.1, 2.5, 0.10)
mar_l2    <- c(3.1, 1.1, 0.1, 0.10)
mar_u3    <- c(0.5, 1.0, 2.5, 0.12)
mar_l3    <- c(3.1, 1.0, 0.1, 0.12)

plot_max  <- 50

ylimu     <- range(sim.dat_s,sim.dat_r,sim.dat_s2r)
yliml     <- range(b_r - hat_mu_r,b_s - hat_mu_s,b_s2r - hat_mu_s2r)
path_plot <- "/home/dom/Dropbox/Forschung/PRJ_OPEN/PRJ_Inference4_FDA_using_RFT/Manuscript/"
## 
pdf(file = paste0(path_plot, "Fig_SIM_DGPs.pdf"), width = width, height = height)
par(mfrow=c(2,1), family = "serif", ps=13, cex.main=1, font.main = 1, mar=mar_u1)
matplot(y = sim.dat_s[,1:plot_max],  x = grid, lwd=.5, col=gray(.75), type="l", lty=1, 
        ylab="", xlab = "t", main="", axes=F, ylim = ylimu)
lines(y=hat_mu_s,x=grid)
legend("bottomright", legend=expression(paste("Estimated mean")), lty=1,bty="n",cex = .95)
axis(2);box()
mtext(text = "Smooth (Cov1)", side = 3, line = .75)
text(x = -.05, y=ylimu[2]*0.95, labels = "Sample Paths (n=100)", cex = .95, pos = 4)
text(x = -.05, y=ylimu[2]*0.8, labels = "(first 50 plotted)", cex = .85, pos = 4)
par(mar=mar_l1)
matplot( y = 0, x = 0, type="n", ylab="", xlab = "", main="", ylim = yliml, xlim = c(0,1))
matlines(x=grid, y=FFSCB_t_band_s - hat_mu_s, col=1, lty=1, lwd=.85)
matlines(x=grid, y=Bs_band_s      - hat_mu_s, col=1, lty=2, lwd=.85)
matlines(x=grid, y=BEc_band_s     - hat_mu_s, col=1, lty=3, lwd=.85)
polygon(x=c(grid,rev(grid)), y=c(naive_t_band_s[,1]-hat_mu_s,rev(naive_t_band_s[,2]-hat_mu_s)), 
        col=gray(.75), border = gray(.75))
legend("center", legend = c("FFSCB-t", "Bootstrap", expression(hat(B)[Ec])), 
       lty=c(1,2,3), bty="n", lwd = c(1), col=c("black"), cex = 1)
legend(x = 0.35, y=yliml[1]*0.65, legend = "naive-t", pch=22, col=gray(.78), pt.bg = gray(.75), pt.cex = 1.75, cex = 1, bty="n")
#legend(x = 0.35, y=-0.12, legend = "naive-t", lty = 1, col=gray(.78), lwd=10, cex = .9, bty="n")
mtext(text = "t", side = 1, line = 1.75)
text(x = -0.05, y=yliml[2]*0.95, labels = "Centered Confidence Bands", cex = .95, pos = 4)
dev.off()
##
pdf(file = paste0(path_plot, "Fig_SIM_DGPr.pdf"), width = width, height = height)
par(mfrow=c(2,1), family = "serif", ps=13, cex.main=1, font.main = 1, mar=mar_u2)
matplot(y = sim.dat_r[,1:plot_max],  x = grid, lwd=.5, col=gray(.75), type="l", lty=1, 
        ylab="", xlab = "t", main="", axes=F, ylim = ylimu)
lines(y=hat_mu_r,x=grid)
#legend("bottomright", legend=expression(paste("Estimated mean")), lty=1,bty="n",cex = .87)
box()
mtext(text = "Rough (Cov2)", side = 3, line = .75)
par(mar=mar_l2)
matplot( y = 0, x = 0, type="n", ylab="", xlab = "", main="", ylim = yliml, xlim = c(0,1), axes=F)
axis(1);box()
matlines(x=grid, y=FFSCB_t_band_r - hat_mu_r, col=1, lty=1, lwd=.85)
matlines(x=grid, y=Bs_band_r      - hat_mu_r, col=1, lty=2, lwd=.85)
matlines(x=grid, y=BEc_band_r     - hat_mu_r, col=1, lty=3, lwd=.85)
polygon(x=c(grid,rev(grid)), y=c(naive_t_band_r[,1]-hat_mu_r,rev(naive_t_band_r[,2]-hat_mu_r)), 
        col=gray(.75), border = gray(.75))
# legend("left", legend = c("FFSCB (t distr.)", "Bootstrap", expression(hat(B)[Ec])), 
#        lty=c(1,2,3), bty="n", lwd = c(1), col=c("black"), cex = .87)
#legend(x = 0, y=-0.12, legend = "naive (t distr.)", pch=22, col=gray(.75), pt.bg = gray(.75), cex = .87, bty="n")
mtext(text = "t", side = 1, line = 1.75)
dev.off()
##
pdf(file = paste0(path_plot, "Fig_SIM_DGPs2r.pdf"), width = width, height = height)
par(mfrow=c(2,1), family = "serif", ps=13, cex.main=1, font.main = 1, mar=mar_u3)
matplot(y = sim.dat_s2r[,1:plot_max],  x = grid, lwd=.5, col=gray(.75), type="l", lty=1, 
        ylab="", xlab = "t", main="", axes=F, ylim = ylimu)
lines(y=hat_mu_s2r,x=grid)
#legend("bottomright", legend=expression(paste("Estimated mean")), lty=1,bty="n",cex = .87)
box()
mtext(text = "Smooth to Rough (Cov3)", side = 3, line = .75)
par(mar=mar_l3)
matplot( y = 0, x = 0, type="n", ylab="", xlab = "", main="", ylim = yliml, xlim = c(0,1), axes=F)
axis(1);box()
matlines(x=grid, y=FFSCB_t_band_s2r- hat_mu_s2r, col=1, lty=1, lwd=.85)
matlines(x=grid, y=Bs_band_s2r     - hat_mu_s2r, col=1, lty=2, lwd=.85)
matlines(x=grid, y=BEc_band_s2r    - hat_mu_s2r, col=1, lty=3, lwd=.85)
polygon(x=c(grid,rev(grid)), y=c(naive_t_band_s2r[,1]-hat_mu_s2r,rev(naive_t_band_s2r[,2]-hat_mu_s2r)), 
        col=gray(.75), border = gray(.75))
# legend("left", legend = c("FFSCB (t distr.)", "Bootstrap", expression(hat(B)[Ec])), 
#        lty=c(1,2,3), bty="n", lwd = c(1), col=c("black"), cex = .87)
#legend(x = 0, y=-0.12, legend = "naive (t distr.)", pch=22, col=gray(.75), pt.bg = gray(.75), cex = .87, bty="n")
mtext(text = "t", side = 1, line = 1.75)
dev.off()





