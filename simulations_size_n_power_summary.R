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
N_seq         <- c(15,100)
p             <- 101
grid          <- make_grid(p, rangevals=c(0,1))
## 
IWT_bool      <- c(FALSE, TRUE)
##
t0            <- 0


## Wrangling
for(IWT in IWT_bool){
  if(IWT){ IWT_SimResults_df <- NULL }else{ SimResults_df <- NULL}
  for(DGP in DGP_seq){# IWT <- TRUE
    for(N in N_seq) {
      if ( N==min(N_seq) ) delta_seq <- delta_Nsmall else delta_seq <- delta_Nlarge
      for(delta in delta_seq) {# DGP <- "DGP3_local"; N <- 100; delta <- 0.08
        ## Load sim_df
        if(IWT){
          load(file = paste0(my_path, "Simulation_Results/IWT_", DGP, "_N=", N, "_Delta=", delta, ".RData"))
        }else{
          load(file = paste0(my_path, "Simulation_Results/",     DGP, "_N=", N, "_alpha=", alpha.level, "_t0=", t0, "_Delta=", delta, ".RData"))
        }
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
            #dplyr::select(band, DGP, N, delta, n_rep, alpha, avg_width, rfrq_excd, rfrq_excd_i1, rfrq_excd_i2) %>% 
            dplyr::select(band, DGP, N, delta, n_rep, alpha, avg_width, rfrq_excd) %>% 
            dplyr::bind_rows(IWT_SimResults_df, .)
        }else{
          SimResults_df <- SimResults_tmp %>% 
            #dplyr::select(band, DGP, N, delta, n_rep, alpha, avg_width, rfrq_excd, rfrq_excd_i1, rfrq_excd_i2) %>% 
            dplyr::select(band, DGP, N, delta, n_rep, alpha, avg_width, rfrq_excd) %>% 
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
SRes_df     <- SimResults_df     %>% mutate(band = as.character(band))
IWT_SRes_df <- IWT_SimResults_df %>% mutate(band = as.character(band))

Size_and_Power_df <- dplyr::bind_rows(SRes_df, IWT_SRes_df) %>% 
  ## back to 'band'-factor
  mutate(band = factor(band, level = c("FFSCB.t", "FFSCB.z", "KR.t", "KR.z", 
                                       "BEc", "Bs", "IWT")),
         DGP  = factor(DGP, levels = c("DGP1_shift", "DGP1_scale", "DGP1_local", 
                                       "DGP2_shift", "DGP2_scale", "DGP2_local", 
                                       "DGP3_shift", "DGP3_scale", "DGP3_local"))) %>% 
  dplyr::arrange(DGP, N, delta)

## ###############################################
## Size 'n' Power df
## ###############################################

## n=15
Size_and_Power_n15_df <- Size_and_Power_df %>% 
  filter(N==15) %>% 
  select(-alpha, -n_rep, -avg_width) %>% 
  spread(delta, rfrq_excd) %>% 
  select(DGP, N, band, `0`:`0.45`) %>%
  mutate_at(4:9, round, 3) %>% 
  arrange(DGP, N, band) 
#Size_and_Power_n15_df %>% print(n=Inf)

## n=100
Size_and_Power_n100_df <- Size_and_Power_df %>% 
  filter(N==100) %>% 
  select(-alpha, -n_rep, -avg_width) %>% 
  spread(delta, rfrq_excd) %>% 
  select(DGP, N, band, `0`:`0.1`) %>% 
  mutate_at(4:9, round, 3) %>% 
  arrange(DGP, N, band) 
#Size_and_Power_n100_df %>% print(n=Inf)

Size_n15_df <- Size_and_Power_n15_df %>% 
  filter(grepl("_shift", Size_and_Power_n15_df$DGP)) %>% 
  select(DGP, band, N, `0`) %>% 
  spread(band, `0`) %>% 
  mutate(DGP     = c("DGP1","DGP2","DGP3")) %>% 
  select(DGP, N, FFSCB.t, KR.t, IWT, BEc, Bs)

Size_n100_df <- Size_and_Power_n100_df %>% 
  filter(grepl("_shift", Size_and_Power_n100_df$DGP)) %>% 
  select(DGP, band, N, `0`) %>% 
  spread(band, `0`) %>% 
  mutate(DGP   = c("DGP1","DGP2","DGP3")) %>% 
  select(DGP, N, FFSCB.t, KR.t, IWT, BEc, Bs)
  
Size_df     <- dplyr::bind_rows(Size_n15_df, Size_n100_df)

AvgWidth_df <- Size_and_Power_df %>% 
  filter(grepl("_shift", Size_and_Power_df$DGP),
         delta==0) %>% 
  select(band, DGP, N, avg_width) %>% 
  spread(band, avg_width) %>% 
  select(-IWT) %>% 
  arrange(N, DGP) %>% 
  mutate(DGP     = c("DGP1","DGP2","DGP3","DGP1","DGP2","DGP3")) %>% 
  mutate_at(3:8, round, 3) %>% 
  rename(FFSCB_t_w = FFSCB.t, 
         FFSCB_z_w = FFSCB.z, 
         KR_t_w    = KR.t,
         KR_z_w    = KR.z,
         BEc_w     = BEc,
         Bs_w      = Bs)

## Table for size and average width:
full_join( Size_df, AvgWidth_df %>% select(-FFSCB_z_w, -KR_z_w))


## ###################
## Power plots
## ###################

mu_shift   <- cbind(meanf_shift(grid, 0), meanf_shift(grid, .1))
mu_scale   <- cbind(meanf_scale(grid, 0), meanf_scale(grid, .1))
mu_local   <- cbind(meanf_rect(grid,  0), meanf_rect(grid,  .1))
##
DGP3_shift <- Size_and_Power_n100_df %>% filter(DGP=="DGP3_shift") %>% select(-DGP,-N) %>% filter(band!="KR.z",band!="FFSCB.z")
DGP3_scale <- Size_and_Power_n100_df %>% filter(DGP=="DGP3_scale") %>% select(-DGP,-N) %>% filter(band!="KR.z",band!="FFSCB.z")
DGP3_local <- Size_and_Power_n100_df %>% filter(DGP=="DGP3_local") %>% select(-DGP,-N) %>% filter(band!="KR.z",band!="FFSCB.z")
##
AvgPWR_Boots_DGP3_shift <- mean(as.numeric(DGP3_shift[grepl("Bs",     DGP3_shift$band),-1])) %>% round(.,d=2)
AvgPWR_BEc_DGP3_shift   <- mean(as.numeric(DGP3_shift[grepl("BEc",    DGP3_shift$band),-1])) %>% round(.,d=2)
AvgPWR_KRt_DGP3_shift   <- mean(as.numeric(DGP3_shift[grepl("KR.t",   DGP3_shift$band),-1])) %>% round(.,d=2)
AvgPWR_FFt_DGP3_shift   <- mean(as.numeric(DGP3_shift[grepl("FFSCB.t",DGP3_shift$band),-1])) %>% round(.,d=2)
AvgPWR_IWT_DGP3_shift   <- mean(as.numeric(DGP3_shift[grepl("IWT",    DGP3_shift$band),-1])) %>% round(.,d=2)
##
AvgPWR_Boots_DGP3_scale <- mean(as.numeric(DGP3_scale[grepl("Bs",     DGP3_scale$band),-1])) %>% round(.,d=2)
AvgPWR_BEc_DGP3_scale   <- mean(as.numeric(DGP3_scale[grepl("BEc",    DGP3_scale$band),-1])) %>% round(.,d=2)
AvgPWR_KRt_DGP3_scale   <- mean(as.numeric(DGP3_scale[grepl("KR.t",   DGP3_scale$band),-1])) %>% round(.,d=2)
AvgPWR_FFt_DGP3_scale   <- mean(as.numeric(DGP3_scale[grepl("FFSCB.t",DGP3_scale$band),-1])) %>% round(.,d=2)
AvgPWR_IWT_DGP3_scale   <- mean(as.numeric(DGP3_scale[grepl("IWT",    DGP3_scale$band),-1])) %>% round(.,d=2)
##
AvgPWR_Boots_DGP3_local <- mean(as.numeric(DGP3_local[grepl("Bs",     DGP3_local$band),-1])) %>% round(.,d=2)
AvgPWR_BEc_DGP3_local   <- mean(as.numeric(DGP3_local[grepl("BEc",    DGP3_local$band),-1])) %>% round(.,d=2)
AvgPWR_KRt_DGP3_local   <- mean(as.numeric(DGP3_local[grepl("KR.t",   DGP3_local$band),-1])) %>% round(.,d=2)
AvgPWR_FFt_DGP3_local   <- mean(as.numeric(DGP3_local[grepl("FFSCB.t",DGP3_local$band),-1])) %>% round(.,d=2)
AvgPWR_IWT_DGP3_local   <- mean(as.numeric(DGP3_local[grepl("IWT",    DGP3_local$band),-1])) %>% round(.,d=2)

## Plots
width     <- 8.5
height    <- 6
mar_u1    <- c(2.5, 2, 2.5, 0.05)
mar_l1    <- c(3.1, 2, 0.5, 0.05)
mar_u2    <- c(2.5, 1.1, 2.5, 0.10)
mar_l2    <- c(3.1, 1.1, 0.5, 0.10)
mar_u3    <- c(2.5, 1.0, 2.5, 0.12)
mar_l3    <- c(3.1, 1.0, 0.5, 0.12)
##
ylimu     <- range(mu_shift, mu_scale, mu_local)
yliml     <- range(DGP3_shift[,-1], DGP3_scale[,-1], DGP3_local[,-1])
##
path_plot <- "/home/dom/Dropbox/Forschung/PRJ_OPEN/PRJ_Inference4_FDA_using_RFT/Manuscript/"
## 
pdf(file = paste0(path_plot, "Fig_SIM_PWR.pdf"), width = width, height = height)
layout(mat = matrix(c(1:6), nrow=2, ncol=3),
       heights = c(1, 2),   # Heights of the two rows
       widths = c(1, 1, 1)) # Widths of the two columns

#layout.show(6)
par(family = "serif", ps=13, cex.main=1, font.main = 1, cex.axis=1.3, mar=mar_u1)
matplot(y = mu_shift,  x = grid, col=1, lty=c(1,2), axes=FALSE,
        ylab="", xlab = "", main="", ylim = ylimu, type="n")
axis(1, at=seq(0,1,len=6))
axis(2, at=seq(0,1,len=6)); box()
matlines(y = mu_shift,  x = grid, col=1, lty=c(1,2))
legend("topleft", title="Mean functions",legend=c(expression(paste(theta," (",Delta==0.1,")")),expression(paste(theta[0]))), lty=c(2,1), bty="n",cex =1.5)
mtext(text = "Shift (Mean1)", side = 3, line = .75)
mtext(text = "t", side = 1, line = 1.75, cex = 1)
##
par(mar=mar_l1)
# FFSCB.t ==  4
# KR.t    ==  3
# BEc     ==  0
# Bs      ==  5
# IWT     ==  6
matplot(x=delta_Nlarge, y=t(DGP3_shift[,-1]), type="b", lty = 1, col=1, pch=c(4,3,0,5,6), axes=FALSE,cex=1.5)
axis(1,at=delta_Nlarge, labels = c("0","0.02","0.04","0.06","0.08","0.1"))
axis(2, at=c(0,0.05,0.2,0.4,0.6,0.8,1),labels = c("0",expression(alpha),"0.2","0.4","0.6","0.8","1"));box()
my_order <- order(c(AvgPWR_FFt_DGP3_shift, 
                    AvgPWR_KRt_DGP3_shift, 
                    AvgPWR_BEc_DGP3_shift,
                    AvgPWR_Boots_DGP3_shift,
                    AvgPWR_IWT_DGP3_shift),decreasing=T)
legend("topleft", title="Power (avg.)", 
       legend = c(paste0("FF-t (",  AvgPWR_FFt_DGP3_shift,")"),
                  paste0("KR-t (",  AvgPWR_KRt_DGP3_shift,")"), 
                  expression(paste(B[Ec]," (0.36)")), # AvgPWR_BEc_DGP3_shift
                  paste0("Boots* (", AvgPWR_Boots_DGP3_shift,")"), 
                  paste0("IWT (",   AvgPWR_IWT_DGP3_shift,")"))[my_order], 
       pch=c(4,3,0,5,6)[my_order], bty="n", col=c("black"), cex = 1.5, pt.cex = 1.5, title.adj = 0.25)
abline(h=0.05); text(x = 0.085, y = 0.02, labels = expression(paste(alpha==0.05)), cex=1.4)
mtext(text = expression(Delta), side = 1, line = 2.05)
##
par(mar=mar_u1)
matplot(y = mu_scale,  x = grid, col=1, lty=c(1,2), 
        ylab="", xlab = "t", main="", axes=F, ylim = ylimu, type="n")
matlines(y = mu_scale,  x = grid, col=1, lty=c(1,2))
legend("topleft", title="Mean functions",legend=c(expression(paste(theta," (",Delta==0.1,")")),expression(paste(theta[0]))), lty=c(2,1), bty="n",cex =1.5)
axis(1, at=seq(0,1,len=6)); box()
mtext(text = "Scale (Mean2)", side = 3, line = .75)
mtext(text = "t", side = 1, line = 1.75)
##
par(mar=mar_l1)
matplot(x=delta_Nlarge, y=t(DGP3_scale[,-1]), type="b", lty = 1, col=1, pch=c(4,3,0,5,6), axes=FALSE,cex=1.5)
axis(1,at=delta_Nlarge, labels = c("0","0.02","0.04","0.06","0.08","0.1"))
box()
my_order <- order(c(AvgPWR_FFt_DGP3_scale, 
                    AvgPWR_KRt_DGP3_scale, 
                    AvgPWR_BEc_DGP3_scale,
                    AvgPWR_Boots_DGP3_scale,
                    AvgPWR_IWT_DGP3_scale),decreasing=T)
legend("topleft", title="Power (avg.)", 
       legend = c(paste0("FF-t (",  AvgPWR_FFt_DGP3_scale,")"), 
                  paste0("KR-t (",  AvgPWR_KRt_DGP3_scale,")"), 
                  expression(paste(B[Ec]," (0.15)")), # AvgPWR_BEc_DGP3_scale
                  paste0("Boots* (", AvgPWR_Boots_DGP3_scale,")"),
                  paste0("IWT (",   AvgPWR_IWT_DGP3_scale,")"))[my_order], 
       pch=c(4,3,0,5,6)[my_order], bty="n", col=c("black"), cex = 1.5, pt.cex = 1.5, title.adj = 0.25)
abline(h=0.05)
mtext(text = expression(Delta), side = 1, line = 2.05)
##
par(mar=mar_u1)
matplot(y = mu_local,  x = grid, col=1, lty=c(1,2), 
        ylab="", xlab = "t", main="", axes=F, ylim = ylimu, type="n")
matlines(y = mu_local,  x = grid, col=1, lty=c(1,2))
legend("topleft", title="Mean functions", legend=c(expression(paste(theta," (",Delta==0.1,")")),expression(paste(theta[0]))), lty=c(2,1), bty="n",cex =1.5)
axis(1, at=seq(0,1,len=6)); box()
mtext(text = "Local (Mean3)", side = 3, line = .75)
mtext(text = "t", side = 1, line = 1.75)
##
par(mar=mar_l1)
matplot(x=delta_Nlarge, y=t(DGP3_local[,-1]), type="b", lty = 1, col=1, pch=c(4,3,0,5,6), axes=FALSE,cex=1.5)
axis(1,at=delta_Nlarge, labels = c("0","0.02","0.04","0.06","0.08","0.1"))
box()
my_order <- order(c(AvgPWR_FFt_DGP3_local, 
                    AvgPWR_KRt_DGP3_local, 
                    AvgPWR_BEc_DGP3_local,
                    AvgPWR_Boots_DGP3_local,
                    AvgPWR_IWT_DGP3_local),decreasing=T)
legend("topleft", title="Power (avg.)", 
       legend = c(paste0("FF-t (",  AvgPWR_FFt_DGP3_local,")"),
                  paste0("KR-t (",  AvgPWR_KRt_DGP3_local,")"), 
                  expression(paste(B[Ec]," (0.30)")), # AvgPWR_BEc_DGP3_local
                  paste0("Boots* (", AvgPWR_Boots_DGP3_local,")"), 
                  paste0("IWT (",   AvgPWR_IWT_DGP3_local,")"))[my_order], 
       pch=c(4,3,0,5,6)[my_order], bty="n", col=c("black"), cex = 1.5, pt.cex = 1.5, title.adj = 0.25)
abline(h=0.05)
mtext(text = expression(Delta), side = 1, line = 2.05)
dev.off()




## ####################
## Fairness
## ####################
p            <- 101
grid         <- make_grid(p, rangevals=c(0,1))
delta        <- 0
mu           <- meanf_shift(grid, delta)
t0           <- 0
alpha.level  <- 0.05
n_int        <- 3
tol          <- .Machine$double.eps^0.5
##
cov.m_s      <- make_cov_m(cov.f = covf.st.matern,    grid=grid, cov.f.params=c(2,        1/4))
cov.m_r      <- make_cov_m(cov.f = covf.st.matern,    grid=grid, cov.f.params=c(1/4,      1/4))
cov.m_s2r    <- make_cov_m(cov.f = covf.nonst.matern, grid=grid, cov.f.params=c(2,   1/4, 1/4))
##
cov.mu_s   <- cov.m_s/N
cov.mu_r   <- cov.m_r/N
cov.mu_s2r <- cov.m_s2r/N
##
tau_s      <- cov2tau_fun(cov.m_s)
tau_r      <- cov2tau_fun(cov.m_r)
tau_s2r    <- cov2tau_fun(cov.m_s2r)
##
b_n15_s    <- make_band_FFSCB_t(x = mu, diag.cov.x = diag(cov.mu_s),  tau = tau_s,  t0 = t0, df =15-1, conf.level = 1-alpha.level, n_int = n_int, tol = tol)
b_n15_r    <- make_band_FFSCB_t(x = mu, diag.cov.x = diag(cov.mu_r),  tau = tau_r,  t0 = t0, df =15-1, conf.level = 1-alpha.level, n_int = n_int, tol = tol)
b_n15_s2r  <- make_band_FFSCB_t(x = mu, diag.cov.x = diag(cov.mu_s2r),tau = tau_s2r,t0 = t0, df =15-1, conf.level = 1-alpha.level, n_int = n_int, tol = tol)
b_n100_s   <- make_band_FFSCB_t(x = mu, diag.cov.x = diag(cov.mu_s),  tau = tau_s,  t0 = t0, df =100-1, conf.level = 1-alpha.level, n_int = n_int, tol = tol)
b_n100_r   <- make_band_FFSCB_t(x = mu, diag.cov.x = diag(cov.mu_r),  tau = tau_r,  t0 = t0, df =100-1, conf.level = 1-alpha.level, n_int = n_int, tol = tol)
b_n100_s2r <- make_band_FFSCB_t(x = mu, diag.cov.x = diag(cov.mu_s2r),tau = tau_s2r,t0 = t0, df =100-1, conf.level = 1-alpha.level, n_int = n_int, tol = tol)
##
sgnf_levels_n15_s    <- b_n15_s$prob_t0    + c(b_n15_s$a_star*(1/3),    b_n15_s$a_star*(2/3),    b_n15_s$a_star)
sgnf_levels_n15_r    <- b_n15_r$prob_t0    + c(b_n15_r$a_star*(1/3),    b_n15_r$a_star*(2/3),    b_n15_r$a_star)
sgnf_levels_n15_s2r  <- b_n15_s2r$prob_t0  + c(b_n15_s2r$a_star*(1/3),  b_n15_s2r$a_star*(2/3),  b_n15_s2r$a_star)
sgnf_levels_n100_s   <- b_n100_s$prob_t0   + c(b_n100_s$a_star*(1/3),   b_n100_s$a_star*(2/3),   b_n100_s$a_star)
sgnf_levels_n100_r   <- b_n100_r$prob_t0   + c(b_n100_r$a_star*(1/3),   b_n100_r$a_star*(2/3),   b_n100_r$a_star)
sgnf_levels_n100_s2r <- b_n100_s2r$prob_t0 + c(b_n100_s2r$a_star*(1/3), b_n100_s2r$a_star*(2/3), b_n100_s2r$a_star)
##
sl_mat <- rbind(sgnf_levels_n15_s,
                sgnf_levels_n15_r,
                sgnf_levels_n15_s2r,
                sgnf_levels_n100_s,
                sgnf_levels_n100_r,
                sgnf_levels_n100_s2r)
##
Size_and_Power_df %>% 
  filter(band=="FFSCB.t", delta==0,
         grepl("_shift",Size_and_Power_df$DGP)) %>% 
  arrange(N, DGP) %>% 
  mutate(esl_1 = round(rfrq_excd_i1,3),
         esl_2 = round(rfrq_excd_i2,3),
         esl_3 = round(rfrq_excd,3),
         sl_1  = round(sl_mat[,1],2),
         sl_2  = round(sl_mat[,2],2),
         sl_3  = round(sl_mat[,3],2)) %>% 
  select(DGP, N, band, esl_1, esl_2, esl_3, sl_1, sl_2, sl_3) 
# select(DGP, N, band, esl_1, sl_1, esl_2, sl_2, esl_3, sl_3) 


Size_and_Power_df %>% 
  filter(band=="KR.t", delta==0,
         grepl("_shift",Size_and_Power_df$DGP)) %>% 
  arrange(N, DGP) %>% 
  mutate(esl_1 = round(rfrq_excd_i1,3),
         esl_2 = round(rfrq_excd_i2,3),
         esl_3 = round(rfrq_excd,3),
         sl_1  = round(sl_mat[,1],2),
         sl_2  = round(sl_mat[,2],2),
         sl_3  = round(sl_mat[,3],2)) %>% 
  select(DGP, N, band, esl_1, esl_2, esl_3, sl_1, sl_2, sl_3) 
# select(DGP, N, band, esl_1, sl_1, esl_2, sl_2, esl_3, sl_3) 

