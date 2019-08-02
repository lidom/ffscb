## #########################################
## Plots of the Data Generating Processes
## #########################################
## Load packages 
library("ffscb")

p            <- 101
grid         <- make_grid(p, rangevals=c(0,1))
delta        <- 0
mu           <- meanf_shift(grid, delta)
cov.m_s      <- make_cov_m(cov.f = covf.st.matern,    grid=grid, cov.f.params=c(3/2,        1/4))
cov.m_r      <- make_cov_m(cov.f = covf.st.matern,    grid=grid, cov.f.params=c(1/2,        1/4))
cov.m_s2r    <- make_cov_m(cov.f = covf.nonst.matern, grid=grid, cov.f.params=c(3/2,   1/2, 1/4))


## Plots:
N            <- 100
n_int        <- 3
seq(0,1,len=n_int+1)
alpha.level  <- 0.05
type         <- c("Bs", "BEc", "naive.t", "FFSCB.t", "KR.t")
seed         <- 2
##
## Smooth
set.seed(seed)
sim.dat_s      <- make_sample(mean.v = mu, cov.m = cov.m_s,   N = N, dist = "rnorm")
t0_s           <- grid[1]
hat_mu_s       <- rowMeans(sim.dat_s)
hat.cov_s      <- crossprod(t(sim.dat_s - hat_mu_s)) / (N-1)
hat.cov.mu_s   <- hat.cov_s / N
hat.tau_s      <- cov2tau_fun(hat.cov_s)# tau_fun(sim.dat_s) 
b_s            <- confidence_band(x=hat_mu_s, cov=hat.cov.mu_s, tau=hat.tau_s, t0=t0_s, df=N-1, 
                                  type=type, conf.level=(1-alpha.level), n_int=n_int)
FFSCB_t_band_s <- b_s[, grepl("FFSCB.t", colnames(b_s))]
naive_t_band_s <- b_s[, grepl("naive.t", colnames(b_s))]
Bs_band_s      <- b_s[, grepl("Bs",      colnames(b_s))]
BEc_band_s     <- b_s[, grepl("BEc",     colnames(b_s))]
KR_t_band_s    <- b_s[, grepl("KR.t",    colnames(b_s))]
##
## Rough  
set.seed(seed)
sim.dat_r      <- make_sample(mean.v = mu, cov.m = cov.m_r,   N = N, dist = "rnorm")
t0_r           <- grid[1]
hat_mu_r       <- rowMeans(sim.dat_r)
hat.cov_r      <- crossprod(t(sim.dat_r - hat_mu_r)) / (N-1)
hat.cov.mu_r   <- hat.cov_r / N
hat.tau_r      <- cov2tau_fun(hat.cov_r)# tau_fun(sim.dat_r) 
b_r            <- confidence_band(x=hat_mu_r, cov=hat.cov.mu_r, tau=hat.tau_r, t0=t0_r, df=N-1, 
                                  type=type, conf.level=(1-alpha.level), n_int=n_int)
FFSCB_t_band_r <- b_r[, grepl("FFSCB.t", colnames(b_r))]
naive_t_band_r <- b_r[, grepl("naive.t", colnames(b_r))]
Bs_band_r      <- b_r[, grepl("Bs",      colnames(b_r))]
BEc_band_r     <- b_r[, grepl("BEc",     colnames(b_r))]
KR_t_band_r    <- b_r[, grepl("KR.t",    colnames(b_s))]
##
## From smooth to rough
set.seed(seed)
sim.dat_s2r    <- make_sample(mean.v = mu, cov.m = cov.m_s2r, N = N, dist = "rnorm")
t0_s2r         <- grid[1]
hat_mu_s2r     <- rowMeans(sim.dat_s2r)
hat.cov_s2r    <- crossprod(t(sim.dat_s2r - hat_mu_s2r)) / (N-1)
hat.cov.mu_s2r <- hat.cov_s2r / N
hat.tau_s2r    <- cov2tau_fun(hat.cov_s2r)# tau_fun(sim.dat_s2r) 
b_s2r          <- confidence_band(x=hat_mu_s2r, cov=hat.cov.mu_s2r, tau=hat.tau_s2r, t0=t0_s2r, df=N-1, 
                                  type=type, conf.level=(1-alpha.level), n_int=n_int)
FFSCB_t_band_s2r <- b_s2r[, grepl("FFSCB.t", colnames(b_s2r))]
naive_t_band_s2r <- b_s2r[, grepl("naive.t", colnames(b_s2r))]
Bs_band_s2r      <- b_s2r[, grepl("Bs",      colnames(b_s2r))]
BEc_band_s2r     <- b_s2r[, grepl("BEc",     colnames(b_s2r))]
KR_t_band_s2r    <- b_s2r[, grepl("KR.t",    colnames(b_s2r))]



## Plots
width     <- 8.5
height    <- 6
mar_u     <- c(0.5, 2, 2.5, 0.05)
mar_l     <- c(4.1, 2, 0.5, 0.05)

plot_max  <- 50

ylimu     <- range(sim.dat_s,sim.dat_r,sim.dat_s2r)
yliml     <- range(b_r - hat_mu_r,b_s - hat_mu_s,b_s2r - hat_mu_s2r)
path_plot <- "/home/dom/Dropbox/Forschung/PRJ_OPEN/PRJ_Inference4_FDA_using_RFT/Manuscript/"
## 
pdf(file = paste0(path_plot, "Fig_SIM_DGP.pdf"), width = width, height = height)
layout(mat = matrix(c(1:6), nrow=2, ncol=3),
       heights = c(1, 1.75), # Heights of the two rows
       widths = c(1, 1, 1))  # Widths of the two columns

#layout.show(6)
par(family = "serif", ps=13, cex.main=1, cex.axis=1.3, font.main = 1, mar=mar_u)
matplot(y = sim.dat_s[,1:plot_max],  x = grid, lwd=.5, col=gray(.75), type="l", lty=1, 
        ylab="", xlab = "t", main="", axes=F, ylim = ylimu)
lines(y=hat_mu_s,x=grid)
legend("bottomright", legend=expression(paste("Estimated mean")), lty=1,bty="n",cex = 1.4)
axis(2);box()
mtext(text = "Smooth (Cov1)", side = 3, line = .95)
text(x = 0, y=ylimu[2]*0.9, labels = "Sample Paths (n=100)", cex =1.5, pos = 4)
##
par(mar=mar_l)
matplot( y = 0, x = 0, type="n", ylab="", xlab = "", main="", ylim = yliml, xlim = c(0,1))
polygon(x=c(grid,rev(grid)), y=c(Bs_band_s[,1]-hat_mu_s,rev(Bs_band_s[,2]-hat_mu_s)), 
        col=gray(.75), border = gray(.75))
polygon(x=c(grid,rev(grid)), y=c(naive_t_band_s[,1]-hat_mu_s,rev(naive_t_band_s[,2]-hat_mu_s)), 
        col=gray(.60), border = gray(.60))
matlines(x=grid, y=FFSCB_t_band_s - hat_mu_s, col=1, lty=1, lwd=.85)
matlines(x=grid, y=KR_t_band_s    - hat_mu_s, col=1, lty=5, lwd=.85)
matlines(x=grid, y=BEc_band_s     - hat_mu_s, col=1, lty=3, lwd=.85)
legend(x = 0,    y=yliml[1]*0.6, legend = c(expression(B[Ec]), "FF-t", "KR-t"), 
       lty=c(3,1,5), bty="n", lwd = c(1), col=c("black"), cex = 1.5)
legend(x = 0.55, y=yliml[1]*0.6, legend = c("Bootstrap", "naive-t"), pch=22, bty="n",
       col=c(gray(.75),gray(.60)), pt.bg = c(gray(.75),gray(.60)), pt.cex = 2.5, cex = 1.5)
mtext(text = "t", side = 1, line = 2.5)
text(x = 0, y=yliml[2]*0.9, labels = "Centered Confidence Bands", cex = 1.5, pos = 4)
##
par(mar=mar_u)
matplot(y = sim.dat_r[,1:plot_max],  x = grid, lwd=.5, col=gray(.75), type="l", lty=1, 
        ylab="", xlab = "t", main="", axes=F, ylim = ylimu)
lines(y=hat_mu_r,x=grid)
box()
mtext(text = "Rough (Cov2)", side = 3, line = .75)
##
par(mar=mar_l)
matplot( y = 0, x = 0, type="n", ylab="", xlab = "", main="", ylim = yliml, xlim = c(0,1), axes=F)
axis(1);box()
polygon(x=c(grid,rev(grid)), y=c(Bs_band_r[,1]-hat_mu_r,rev(Bs_band_r[,2]-hat_mu_r)), 
        col=gray(.75), border = gray(.75))
polygon(x=c(grid,rev(grid)), y=c(naive_t_band_r[,1]-hat_mu_r,rev(naive_t_band_r[,2]-hat_mu_r)), 
        col=gray(.60), border = gray(.60))
matlines(x=grid, y=FFSCB_t_band_r - hat_mu_r, col=1, lty=1, lwd=.85)
matlines(x=grid, y=KR_t_band_r    - hat_mu_r, col=1, lty=5, lwd=.85)
matlines(x=grid, y=BEc_band_r     - hat_mu_r, col=1, lty=3, lwd=.85)
mtext(text = "t", side = 1, line = 2.5)
##
par(mar=mar_u)
matplot(y = sim.dat_s2r[,1:plot_max],  x = grid, lwd=.5, col=gray(.75), type="l", lty=1, 
        ylab="", xlab = "t", main="", axes=F, ylim = ylimu)
lines(y=hat_mu_s2r,x=grid)
box()
mtext(text = "Smooth to Rough (Cov3)", side = 3, line = .75)
##
par(mar=mar_l)
matplot( y = 0, x = 0, type="n", ylab="", xlab = "", main="", ylim = yliml, xlim = c(0,1), axes=F)
axis(1);box()
polygon(x=c(grid,rev(grid)), y=c(Bs_band_s2r[,1]-hat_mu_s2r,rev(Bs_band_s2r[,2]-hat_mu_s2r)), 
        col=gray(.75), border = gray(.75))
polygon(x=c(grid,rev(grid)), y=c(naive_t_band_s2r[,1]-hat_mu_s2r,rev(naive_t_band_s2r[,2]-hat_mu_s2r)), 
        col=gray(.60), border = gray(.60))
matlines(x=grid, y=FFSCB_t_band_s2r - hat_mu_s2r, col=1, lty=1, lwd=.85)
matlines(x=grid, y=KR_t_band_s2r    - hat_mu_s2r, col=1, lty=5, lwd=.85)
matlines(x=grid, y=BEc_band_s2r     - hat_mu_s2r, col=1, lty=3, lwd=.85)
mtext(text = "t", side = 1, line = 2.5)
dev.off()



## ################################## ## ################################## ## ##################################
## Fragments                          ## Fragments                          ## Fragments 
## ################################## ## ################################## ## ##################################


## Plots:
N            <- 500
n_int        <- 3
seq(0,1,len=n_int+1)
alpha.level  <- 0.05
type         <- c("naive.t", "FFSCB.t", "KR.t")
##

# N_frag_mat     <- matrix(NA,p,1000)
# for(r in 1:1000){
# sim.dat_s      <- make_fragm_sample(mean.v = mu, cov.m = cov.m_s, fragm_len = 40,  N = N, dist = "rnorm")
# N_frag_mat[,r] <- apply(sim.dat_s$X_frag_mat, 1, function(x) length(c(na.omit(x))))
# }
# 
# range(N_frag_mat)

## Smooth
seed         <- 2
set.seed(seed)
sim.dat_s      <- make_fragm_sample(mean.v = mu, cov.m = cov.m_s, fragm_len = 40,  N = N, dist = "rnorm")
t0_s           <- grid[1]
hat_mu_s       <- rowMeans(sim.dat_s$X_frag_mat, na.rm = TRUE)
hat.cov_s      <- ffscb:::cov_partial_fd(sim.dat_s$X_frag_mat)
N_frag_s       <- apply(sim.dat_s$X_frag_mat, 1, function(x) length(c(na.omit(x))))
hat.cov.mu_s   <- hat.cov_s / N_frag_s
hat.tau_s      <- cov2tau_fun(hat.cov_s)
hat.tau_s      <- tau_fun(sim.dat_s) 
b_s            <- ffscb:::confidence_band_fragm(x=hat_mu_s, diag.cov.x = diag(hat.cov.mu_s), tau=hat.tau_s, t0=t0_s, df=min(N_frag_s)-1, 
                                  type=type, conf.level=(1-alpha.level), n_int=n_int)
FFSCB_t_band_s <- b_s[, grepl("FFSCB.t", colnames(b_s))]
naive_t_band_s <- b_s[, grepl("naive.t", colnames(b_s))]
KR_t_band_s    <- b_s[, grepl("KR.t",    colnames(b_s))]
##
## Rough  
set.seed(seed)
sim.dat_r      <- make_fragm_sample(mean.v = mu, cov.m = cov.m_r, fragm_len = 40,  N = N, dist = "rnorm")
t0_r           <- grid[1]
hat_mu_r       <- rowMeans(sim.dat_r$X_frag_mat, na.rm = TRUE)
hat.cov_r      <- ffscb:::cov_partial_fd(sim.dat_r$X_frag_mat)
N_frag_r       <- apply(sim.dat_r$X_frag_mat, 1, function(x) length(c(na.omit(x))))
hat.cov.mu_r   <- hat.cov_r / N_frag_r
hat.tau_r      <- cov2tau_fun(hat.cov_r)
hat.tau_r      <- tau_fun(sim.dat_r) 
b_r            <- ffscb:::confidence_band_fragm(x=hat_mu_r, diag.cov.x = diag(hat.cov.mu_r), tau=hat.tau_r, t0=t0_r, df=min(N_frag_r)-1, 
                                                type=type, conf.level=(1-alpha.level), n_int=n_int)
FFSCB_t_band_r <- b_r[, grepl("FFSCB.t", colnames(b_r))]
naive_t_band_r <- b_r[, grepl("naive.t", colnames(b_r))]
KR_t_band_r    <- b_r[, grepl("KR.t",    colnames(b_r))]
##
## From smooth to rough
set.seed(seed)
sim.dat_s2r      <- make_fragm_sample(mean.v = mu, cov.m = cov.m_s2r, fragm_len = 40,  N = N, dist = "rnorm")
t0_s2r           <- grid[1]
hat_mu_s2r       <- rowMeans(sim.dat_s2r$X_frag_mat, na.rm = TRUE)
hat.cov_s2r      <- ffscb:::cov_partial_fd(sim.dat_s2r$X_frag_mat)
N_frag_s2r       <- apply(sim.dat_s2r$X_frag_mat, 1, function(x) length(c(na.omit(x))))
hat.cov.mu_s2r   <- hat.cov_s2r / N_frag_s2r
hat.tau_s2r      <- cov2tau_fun(hat.cov_s2r)
hat.tau_s2r      <- tau_fun(sim.dat_s2r) 
b_s2r            <- ffscb:::confidence_band_fragm(x=hat_mu_s2r, diag.cov.x = diag(hat.cov.mu_s2r), tau=hat.tau_s2r, t0=t0_s2r, df=min(N_frag_s2r)-1, 
                                                type=type, conf.level=(1-alpha.level), n_int=n_int)
FFSCB_t_band_s2r <- b_s2r[, grepl("FFSCB.t", colnames(b_s2r))]
naive_t_band_s2r <- b_s2r[, grepl("naive.t", colnames(b_s2r))]
KR_t_band_s2r    <- b_s2r[, grepl("KR.t",    colnames(b_s2r))]



## Plots
width     <- 8.5
height    <- 6
mar_u     <- c(0.5, 2, 2.5, 0.05)
mar_l     <- c(4.1, 2, 0.5, 0.05)

plot_max  <- 50

ylimu     <- range(sim.dat_s$X_frag_mat,sim.dat_r$X_frag_mat,sim.dat_s2r$X_frag_mat,na.rm = T)
yliml     <- range(b_r - hat_mu_r,b_s - hat_mu_s,b_s2r - hat_mu_s2r)
path_plot <- "/home/dom/Dropbox/Forschung/PRJ_OPEN/PRJ_Inference4_FDA_using_RFT/Manuscript/"
## 
pdf(file = paste0(path_plot, "Fig_SIM_DGP_FRGM.pdf"), width = width, height = height)
layout(mat = matrix(c(1:6), nrow=2, ncol=3),
       heights = c(1, 1.75), # Heights of the two rows
       widths = c(1, 1, 1))  # Widths of the two columns

#layout.show(6)
par(family = "serif", ps=13, cex.main=1, cex.axis=1.3, font.main = 1, mar=mar_u)
matplot(y = sim.dat_s$X_frag_mat[,1:plot_max],  x = sim.dat_s$grid_frag_mat[,1:plot_max], lwd=.5, col=gray(.75), type="l", lty=1, 
        ylab="", xlab = "t", main="", axes=F, ylim = ylimu)
lines(y=hat_mu_s,x=grid)
legend("bottomright", legend=expression(paste("Estimated mean")), lty=1,bty="n",cex = 1.4)
axis(2);box()
mtext(text = "Smooth (Cov1)", side = 3, line = .95)
text(x = 0, y=ylimu[2]*0.9, labels = "Sample Paths (n=100)", cex =1.5, pos = 4)
##
par(mar=mar_l)
matplot( y = 0, x = 0, type="n", ylab="", xlab = "", main="", ylim = yliml, xlim = c(0,1))
polygon(x=c(grid,rev(grid)), y=c(naive_t_band_s[,1]-hat_mu_s,rev(naive_t_band_s[,2]-hat_mu_s)), 
        col=gray(.60), border = gray(.60))
matlines(x=grid, y=FFSCB_t_band_s - hat_mu_s, col=1, lty=1, lwd=.85)
matlines(x=grid, y=KR_t_band_s    - hat_mu_s, col=1, lty=5, lwd=.85)
legend(x = 0,    y=yliml[1]*0.65, legend = c("FF-t", "KR-t"), 
       lty=c(1,5), bty="n", lwd = c(1), col=c("black"), cex = 1.5)
legend(x = 0.55, y=yliml[1]*0.65, legend = c("naive-t"), pch=22, bty="n",
       col=c(gray(.60)), pt.bg = c(gray(.60)), pt.cex = 2.5, cex = 1.5)
mtext(text = "t", side = 1, line = 2.5)
text(x = 0, y=yliml[2]*0.9, labels = "Centered Confidence Bands", cex = 1.5, pos = 4)
##
par(mar=mar_u)
matplot(y = sim.dat_r$X_frag_mat[,1:plot_max],  x = sim.dat_r$grid_frag_mat[,1:plot_max], lwd=.5, col=gray(.75), type="l", lty=1, 
        ylab="", xlab = "t", main="", axes=F, ylim = ylimu)
lines(y=hat_mu_r,x=grid)
box()
mtext(text = "Rough (Cov2)", side = 3, line = .75)
##
par(mar=mar_l)
matplot( y = 0, x = 0, type="n", ylab="", xlab = "", main="", ylim = yliml, xlim = c(0,1), axes=F)
axis(1);box()
polygon(x=c(grid,rev(grid)), y=c(naive_t_band_r[,1]-hat_mu_r,rev(naive_t_band_r[,2]-hat_mu_r)), 
        col=gray(.60), border = gray(.60))
matlines(x=grid, y=FFSCB_t_band_r - hat_mu_r, col=1, lty=1, lwd=.85)
matlines(x=grid, y=KR_t_band_r    - hat_mu_r, col=1, lty=5, lwd=.85)
mtext(text = "t", side = 1, line = 2.5)
##
par(mar=mar_u)
matplot(y = sim.dat_s2r$X_frag_mat[,1:plot_max],  x = sim.dat_s2r$grid_frag_mat[,1:plot_max], lwd=.5, col=gray(.75), type="l", lty=1, 
        ylab="", xlab = "t", main="", axes=F, ylim = ylimu)
lines(y=hat_mu_s2r,x=grid)
box()
mtext(text = "Smooth to Rough (Cov3)", side = 3, line = .75)
##
par(mar=mar_l)
matplot( y = 0, x = 0, type="n", ylab="", xlab = "", main="", ylim = yliml, xlim = c(0,1), axes=F)
axis(1);box()
polygon(x=c(grid,rev(grid)), y=c(naive_t_band_s2r[,1]-hat_mu_s2r,rev(naive_t_band_s2r[,2]-hat_mu_s2r)), 
        col=gray(.60), border = gray(.60))
matlines(x=grid, y=FFSCB_t_band_s2r - hat_mu_s2r, col=1, lty=1, lwd=.85)
matlines(x=grid, y=KR_t_band_s2r    - hat_mu_s2r, col=1, lty=5, lwd=.85)
mtext(text = "t", side = 1, line = 2.5)
dev.off()

