## #######################################################
##
## Summarizing the simulation results
##
## #######################################################

library("ffscb")
library("tidyverse")
## 

## path to simulation results:
my_path       <- "/home/dom/Dropbox/Forschung/PRJ_OPEN/PRJ_Inference4_FDA_using_RFT/"

## Looping-Variables
DGP_seq       <- c("DGP1_shift","DGP1_scale","DGP1_local",
                   "DGP2_shift","DGP2_scale","DGP2_local", 
                   "DGP3_shift","DGP3_scale","DGP3_local")
##
delta_Nsmall  <- c(0, seq(from = 0.05, to = 0.45, len = 5))
delta_Nlarge  <- c(0, seq(from = 0.02, to = 0.1,  len = 5))
alpha.level   <- 0.05
N_seq         <- c(10,100)
p             <- 101
grid          <- make_grid(p, rangevals=c(0,1))
## 
IWT_bool      <- c(FALSE, TRUE)


## Wrangling
for(IWT in IWT_bool){
  if(IWT){ IWT_SimResults_df <- NULL }else{ SimResults_df <- NULL}
  for(DGP in DGP_seq){# IWT <- TRUE
    for(N in N_seq) {
      if ( N==min(N_seq) ) delta_seq <- delta_Nsmall else delta_seq <- delta_Nlarge
      for(delta in delta_seq) {# DGP <- "DGP3_shift"; N <- 10; delta <- 0
        ## Load sim_df
        if(IWT){
          load(file = paste0(my_path, "Simulation_Results/IWT_", DGP, "_N=", N, "_Delta=", delta, ".RData"))
        }else{
          load(file = paste0(my_path, "Simulation_Results/",     DGP, "_N=", N, "_Delta=", delta, ".RData"))
        }
        ## Compute which share of the difference between mu and mu0 was correctly found
        if(grepl("shift", DGP)) { mu0 <- meanf_shift(grid, 0);  mu <- meanf_shift(grid, delta) }
        if(grepl("scale", DGP)) { mu0 <- meanf_scale(grid, 0);  mu <- meanf_scale(grid, delta) }
        if(grepl("local", DGP)) { mu0 <- meanf_rect( grid, 0);  mu <- meanf_rect( grid, delta) }
        ##
        SimResults_tmp <- sim_df %>% 
          dplyr::group_by(band) %>% 
          dplyr::summarise(rfrq_excd    = mean(excd),
                           #rfrq_excd_i1 = mean(excd_i1),
                           #rfrq_excd_i2 = mean(excd_i2),
                           avg_width    = mean(wdth),
                           n_rep        = unique(sim_df$n_rep),
                           DGP          = unique(sim_df$DGP),
                           delta        = unique(sim_df$delta),
                           N            = unique(sim_df$N),
                           alpha        = alpha.level) 
        ##
        ## Row-Binding all 'SimResults_tmp' data frames:
        if(IWT){
          IWT_SimResults_df <- SimResults_tmp %>% 
            dplyr::select(band, DGP, N, delta, n_rep, alpha, rfrq_excd) %>% 
            dplyr::bind_rows(IWT_SimResults_df, .)
        }else{
          SimResults_df <- SimResults_tmp %>% 
            dplyr::select(band, DGP, N, delta, n_rep, alpha, rfrq_excd) %>% 
            dplyr::bind_rows(SimResults_df, .)
        }
      }
    }
  }
}

## ###############################################
## Joining the aggregated simulation results
## ###############################################

## 'band'-factor to 'band'-character (needed for the rowbinding)
SimResults_df     <- SimResults_df     %>% mutate(band = as.character(band))
IWT_SimResults_df <- IWT_SimResults_df %>% mutate(band = as.character(band))

Size_and_Power_df <- dplyr::bind_rows(SimResults_df, IWT_SimResults_df) %>% 
  ## back to 'band'-factor
  mutate(band = as.factor(band)) %>% 
  dplyr::arrange(DGP, N, delta)

## ###############################################
## Building data frames  
## ###############################################

## n=10
Size_and_Power_n10_df <- Size_and_Power_df %>% 
  filter(band!="naive.t" & N==10) %>% 
  select(-alpha, -n_rep, -N) %>% 
  spread(delta, rfrq_excd) %>% 
  select(DGP, band, `0`:`0.45`) %>%
  arrange(DGP, band) 

Size_and_Power_n10_df %>% print(n=Inf)

Size_n10_df <- Size_and_Power_n10_df %>% 
  filter(grepl("_local", Size_and_Power_n10_df$DGP)) %>% 
  select(DGP, band, `0`) %>% 
  spread(band, `0`) %>% 
  mutate(DGP=c("DGP1","DGP2","DGP3")) %>% 
  select(DGP, FFSCB.t, FFSCB.z, KR.t, KR.z, IWT, BEc, Bs)

## n=100
Size_and_Power_n100_df <- Size_and_Power_df %>% 
  filter(band!="naive.t" & N==100) %>% 
  select(-alpha, -n_rep, -N) %>% 
  spread(delta, rfrq_excd) %>% 
  select(DGP, band, `0`:`0.1`) %>% 
  arrange(DGP, band) 

Size_and_Power_n100_df %>% print(n=Inf)

Size_n100_df <- Size_and_Power_n100_df %>% 
  filter(grepl("_local", Size_and_Power_n100_df$DGP)) %>% 
  select(DGP, band, `0`) %>% 
  spread(band, `0`) %>% 
  mutate(DGP=c("DGP1","DGP2","DGP3")) %>% 
  select(DGP, FFSCB.t, FFSCB.z, KR.t, KR.z, IWT, BEc, Bs)

## Plots

Size_and_Power_n100_df




## #########################################
## Plots of the Data Generating Processes
## #########################################

cov2tau_fun <- function(cov_mat, df=NULL){
  ## 'cov_mat' denotes the sample covariance function computed from the sample functions X_1(t),...X_n(t)
  ## Caution: assumed grid is in [0,1]
  ##
  if(is.null(df)           ){scl <- 1}
  if(is.numeric(df) & df >2){scl <- (df-2)/df}else{scl <- 1}
  ##
  p        <- ncol(cov_mat)
  grid     <- seq(0,1,len=p)
  corr_mat <- cov2cor(cov_mat) * scl
  ## computing the numeric approximation to c_12(t,s) := \partial^2 c(t,s)/(\partial t \partial t) 
  ## with c_12(t,s) evaluated at t=t and s=t
  a1  <- corr_mat[cbind(2:p,2:p)]        # corr(t+h, t+h) with h=0.005, and t\in{0.005,0.015,0.025,...,0.995}
  a2  <- corr_mat[cbind(1:(p-1),2:p)]    # corr(t-h, t+h) 
  a3  <- corr_mat[cbind(1:(p-1),1:(p-1))]# corr(t-h, t-h) 
  ##
  h   <- diff(grid)[1]/2
  tau <- sqrt( c(a1- 2*a2 + a3) /( 4* h^2 ) )
  ## tau has length p-1, so we interpolate to get length p
  xx  <- c(0,      grid[-p] + h, 1               )
  yy  <- c(tau[1], tau,          tau[length(tau)])
  tau <- spline(y=yy, x=xx, xout = grid, method = "natural")$y
  ##
  return(tau)
}

## Load packages 
library("ffscb")
##
p            <- 101
grid         <- make_grid(p, rangevals=c(0,1))
delta        <- 0
mu           <- meanf_shift(grid, delta)
##
cov.m_s      <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(3/2, 1, 1/4))
cov.m_r      <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(1/2, 1, 1/4))
cov.m_s2r    <- make_cov_m(cov.f = covf.st.matern.warp.power, grid=grid, cov.f.params=c(4/5, 1, 1/4, 4))

## Plots:
N            <- 50
n_int        <- 10
alpha.level  <- 0.05
type         <- c("Bs", "BEc", "naive.t", "FFSCB.t")

##
## Smooth
set.seed(1)
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
true.cov.tau_s <- cov2tau_fun(cov.m_s, df = (N-1) )
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


matplot(y=cbind(FFSCB_diff_s,true.cov.FFSCB_diff_s_OLD, true.cov.FFSCB_diff_s), x=grid, type="l")


##
## Rough  
set.seed(1)
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
set.seed(1)
sim.dat_s2r    <- make_sample(mean.v = mu, cov.m = cov.m_s2r, N = N, dist = "rnorm")
t0_s2r         <- grid[p]
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







## Plots
width     <- 3.2
height    <- 6.2
mar_u1    <- c(0.5, 1.9, 2.5, 0.05)
mar_l1    <- c(3.1, 1.9, 0.1, 0.05)
mar_u2    <- c(0.5, 1.1, 2.5, 0.10)
mar_l2    <- c(3.1, 1.1, 0.1, 0.10)
mar_u3    <- c(0.5, 1.0, 2.5, 0.12)
mar_l3    <- c(3.1, 1.0, 0.1, 0.12)

ylimu     <- range(sim.dat_s,sim.dat_r,sim.dat_s2r)
yliml     <- range(b_r - hat_mu_r,b_s - hat_mu_s,b_s2r - hat_mu_s2r)
path_plot <- "/home/dom/Dropbox/Forschung/PRJ_OPEN/PRJ_Inference4_FDA_using_RFT/Manuscript/"
## 
pdf(file = paste0(path_plot, "Fig_SIM_DGPs.pdf"), width = width, height = height)
par(mfrow=c(2,1), family = "serif", ps=13, cex.main=.99, font.main = 1, mar=mar_u1)
matplot(y = sim.dat_s,  x = grid, lwd=.5, col=gray(.25), type="l", lty=1, 
        ylab="", xlab = "t", main="", axes=F, ylim = ylimu)
axis(2);box()
mtext(text = "Smooth (Cov1)", side = 3, line = .75)
text(x = 0, y=1.5, labels = "Sample Paths", cex = .95, pos = 4)
par(mar=mar_l1)
matplot( y = 0, x = 0, type="n", ylab="", xlab = "", main="", ylim = yliml, xlim = c(0,1))
matlines(x=grid, y=FFSCB_t_band_s - hat_mu_s, col=1, lty=1, lwd=.75)
matlines(x=grid, y=Bs_band_s      - hat_mu_s, col=1, lty=2, lwd=.75)
matlines(x=grid, y=naive_t_band_s - hat_mu_s, col=1, lty=3, lwd=.75)
legend("left", legend = c("FFSCB (t distr.)", "Bootstrap", "naive (t distr.)"), 
       lty=c(1,2,3), bty="n", lwd = c(1), col=c("black"), cex = .87, seg.len=1)
mtext(text = "t", side = 1, line = 1.75)
text(x = 0, y=0.15, labels = "Centered Confidence Bands", cex = .95, pos = 4)
dev.off()
##
pdf(file = paste0(path_plot, "Fig_SIM_DGPr.pdf"), width = width, height = height)
par(mfrow=c(2,1), family = "serif", ps=13, cex.main=.99, font.main = 1, mar=mar_u2)
matplot(y = sim.dat_r,  x = grid, lwd=.5, col=gray(.25), type="l", lty=1, 
        ylab="", xlab = "t", main="", axes=F, ylim = ylimu)
box()
mtext(text = "Rough (Cov2)", side = 3, line = .75)
par(mar=mar_l2)
matplot( y = 0, x = 0, type="n", ylab="", xlab = "", main="", ylim = yliml, xlim = c(0,1), axes=F)
axis(1);box()
matlines(x=grid, y=FFSCB_t_band_r - hat_mu_r, col=1, lty=1, lwd=.75)
matlines(x=grid, y=Bs_band_r      - hat_mu_r,      col=1, lty=2, lwd=.75)
matlines(x=grid, y=naive_t_band_r - hat_mu_r, col=1, lty=3, lwd=.75)
legend("left", legend = c("FFSCB (t distr.)", "Bootstrap", "naive (t distr.)"), 
       lty=c(1,2,3), bty="n", lwd = c(1), col=c("black"), cex = .87, seg.len=2)
mtext(text = "t", side = 1, line = 1.75)
dev.off()
##
pdf(file = paste0(path_plot, "Fig_SIM_DGPs2r.pdf"), width = width, height = height)
par(mfrow=c(2,1), family = "serif", ps=13, cex.main=.99, font.main = 1, mar=mar_u3)
matplot(y = sim.dat_s2r,  x = grid, lwd=.5, col=gray(.25), type="l", lty=1, 
        ylab="", xlab = "t", main="", axes=F, ylim = ylimu)
box()
mtext(text = "Smooth to Rough (Cov3)", side = 3, line = .75)
par(mar=mar_l3)
matplot( y = 0, x = 0, type="n", ylab="", xlab = "", main="", ylim = yliml, xlim = c(0,1), axes=F)
axis(1);box()
matlines(x=grid, y=FFSCB_t_band_s2r- hat_mu_s2r, col=1, lty=1, lwd=.75)
matlines(x=grid, y=Bs_band_s2r     - hat_mu_s2r, col=1, lty=2, lwd=.75)
matlines(x=grid, y=naive_t_band_s2r- hat_mu_s2r, col=1, lty=3, lwd=.75)
legend("left", legend = c("FFSCB (t distr.)", "Bootstrap", "naive (t distr.)"), 
       lty=c(1,2,3), bty="n", lwd = c(1), col=c("black"), cex = .87, seg.len=2)
mtext(text = "t", side = 1, line = 1.75)
dev.off()





