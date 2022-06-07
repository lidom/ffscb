library("ffscb")
library("tidyverse")

fragm_df <- spnbmd        %>% 
  filter(sex == "fem")    %>% 
  group_by(idnum)         %>% 
  filter(length(idnum)>1) %>% 
  ungroup()
source(file = "Applications/Fragments_Data_Preparation.R")
Y_Wf_mat <- Y_W_mat[,apply(rbind(Y_W_mat,X_W_mat),2,function(x){!all(is.na(x))})]
X_Wf_mat <- X_W_mat[,apply(rbind(Y_W_mat,X_W_mat),2,function(x){!all(is.na(x))})]
N_Wf     <- N_W
Y_Af_mat <- Y_A_mat[,apply(rbind(Y_A_mat,X_A_mat),2,function(x){!all(is.na(x))})]
X_Af_mat <- X_A_mat[,apply(rbind(Y_A_mat,X_A_mat),2,function(x){!all(is.na(x))})]
N_Af     <- N_A
Y_Hf_mat <- Y_H_mat[,apply(rbind(Y_H_mat,X_H_mat),2,function(x){!all(is.na(x))})]
X_Hf_mat <- X_H_mat[,apply(rbind(Y_H_mat,X_H_mat),2,function(x){!all(is.na(x))})]
N_Hf     <- N_H
Y_Bf_mat <- Y_B_mat[,apply(rbind(Y_B_mat,X_B_mat),2,function(x){!all(is.na(x))})]
X_Bf_mat <- X_B_mat[,apply(rbind(Y_B_mat,X_B_mat),2,function(x){!all(is.na(x))})]
N_Bf     <- N_B

fragm_df <- spnbmd        %>% 
  filter(sex == "mal")    %>% 
  group_by(idnum)         %>% 
  filter(length(idnum)>1) %>% 
  ungroup()
source(file = "Applications/Fragments_Data_Preparation.R")
Y_Wm_mat <- Y_W_mat[,apply(rbind(Y_W_mat,X_W_mat),2,function(x){!all(is.na(x))})]
X_Wm_mat <- X_W_mat[,apply(rbind(Y_W_mat,X_W_mat),2,function(x){!all(is.na(x))})]
N_Wm     <- N_W
Y_Am_mat <- Y_A_mat[,apply(rbind(Y_A_mat,X_A_mat),2,function(x){!all(is.na(x))})]
X_Am_mat <- X_A_mat[,apply(rbind(Y_A_mat,X_A_mat),2,function(x){!all(is.na(x))})] 
N_Am     <- N_A
Y_Hm_mat <- Y_H_mat[,apply(rbind(Y_H_mat,X_H_mat),2,function(x){!all(is.na(x))})]
X_Hm_mat <- X_H_mat[,apply(rbind(Y_H_mat,X_H_mat),2,function(x){!all(is.na(x))})]
N_Hm     <- N_H
Y_Bm_mat <- Y_B_mat[,apply(rbind(Y_B_mat,X_B_mat),2,function(x){!all(is.na(x))})]
X_Bm_mat <- X_B_mat[,apply(rbind(Y_B_mat,X_B_mat),2,function(x){!all(is.na(x))})]
N_Bm     <- N_B


(N_f <- sum(c(N_Wf, N_Af, N_Hf, N_Bf)))
(N_m <- sum(c(N_Wm, N_Am, N_Hm, N_Bm)))

round(
  cbind(c(N_Wf, N_Af, N_Hf, N_Bf)/N_f,
        c(N_Wm, N_Am, N_Hm, N_Bm)/N_m)
  , d=2)


## Function for computing the local numbers of observations
n_ts <- function(X_mat) 
{
  p           <- nrow(X_mat)
  n           <- ncol(X_mat)
  n_ts_mat    <- matrix(NA, ncol = p, nrow = p)
  for (s in seq(1, p)) {
    for (t in seq(s, p)) {
      X_s <- X_mat[s, ]
      X_t <- X_mat[t, ]
      n_na <- sum(is.na(c(X_s * X_t)))
      if (n - n_na == 0) {
        n_ts_mat[s, t] <- 0
      }
      else {
        n_ts_mat[s, t] <- length(c(na.omit(X_s * X_t)))
      }
      n_ts_mat[t, s] <- n_ts_mat[s, t]
    }
  }
  return(n_ts_mat)
}


## ########################
Y_f_mat  <- cbind(Y_Wf_mat, Y_Af_mat, Y_Hf_mat, Y_Bf_mat)
X_f_mat  <- cbind(X_Wf_mat, X_Af_mat, X_Hf_mat, X_Bf_mat)
Y_m_mat  <- cbind(Y_Wm_mat, Y_Am_mat, Y_Hm_mat, Y_Bm_mat)
X_m_mat  <- cbind(X_Wm_mat, X_Am_mat, X_Hm_mat, X_Bm_mat)
##
Y_f_mat  <- Y_f_mat[,apply(X_f_mat,2,function(x){length(c(na.omit(x)))>=2})]
X_f_mat  <- X_f_mat[,apply(X_f_mat,2,function(x){length(c(na.omit(x)))>=2})]
Y_m_mat  <- Y_m_mat[,apply(X_m_mat,2,function(x){length(c(na.omit(x)))>=2})]
X_m_mat  <- X_m_mat[,apply(X_m_mat,2,function(x){length(c(na.omit(x)))>=2})]
##
(N_f      <- ncol(X_f_mat))
(N_m      <- ncol(X_m_mat))
##
apply(cbind(X_f_mat,X_m_mat),2,function(x){ max(x, na.rm = T) - min(x,na.rm = T)}) %>% summary()
##
hat_mu_f   <- rowMeans(Y_f_mat, na.rm = TRUE)
hat_mu_m   <- rowMeans(Y_m_mat, na.rm = TRUE)
##
n_f_ts     <- n_ts(Y_f_mat);  n_f <- diag(n_f_ts)
n_m_ts     <- n_ts(Y_m_mat);  n_m <- diag(n_m_ts)
##
cov_f      <- ffscb:::cov_partial_fd(Y_f_mat)
cov_m      <- ffscb:::cov_partial_fd(Y_m_mat)
cov_mat    <- ( (n_f_ts - 1) * cov_f + (n_m_ts - 1) * cov_m ) / (n_f_ts + n_m_ts -2)
diag_cov   <- diag(cov_mat)
diag_cov_x <- diag_cov / (n_f + n_m)
##
tau        <- ffscb:::cov2tau_fun(cov_mat)
##
b    <- ffscb:::confidence_band_fragm(x          = c(hat_mu_f - hat_mu_m), 
                                      diag.cov.x = diag_cov_x, 
                                      tau        = tau, 
                                      df         = min(c(n_f,n_m)) - 1, 
                                      type       = c("FFSCB.t"), 
                                      conf.level = (1-0.05), 
                                      n_int      = 2)

slct_FFSCB_t <- grepl(pattern = "FFSCB.t", colnames(b))

hat_mu_diff  <- b[,1]
FFSCB_t_band <- b[,slct_FFSCB_t] #FFSCB_res$band[,-1]


grid[FFSCB_t_band[,2]>0] %>% range


width     <- 7
height    <- 3.8
cex       <- 1
cexs      <- 0.95
mar       <- c(3.25, 2.1, 2.1, 2.1)
##
par(mfrow=c(1,2), family = "serif", ps=13, cex.main=.99, font.main = 1, mar=mar)
matplot(y=cbind(Y_f_mat, Y_m_mat),  
        x=cbind(X_f_mat, X_m_mat), type="n", xlab="", ylab="", ylim=c(.45,max(cbind(Y_f_mat, Y_m_mat),na.rm=T)))
matlines(y=Y_f_mat,  x=X_f_mat, lty=1, col=gray(.5))
matlines(y=hat_mu_f, x=grid,    lty=1, col="black", lwd = 2)
matlines(y=Y_m_mat,  x=X_m_mat, lty=2, col=gray(.5))
matlines(y=hat_mu_m, x=grid,    lty=2, col="black", lwd = 2)
legend("topleft", legend = c("Female", "Male"), 
       lty = c(1,2), col = c(gray(.5), gray(.5)), lwd=1.5, bty="n", cex = cexs)
legend("bottomright", legend = c("Estimated mean (female)",
                                 "Estimated mean (male)"), 
       lty = c(1,2), col = c("black", "black"), lwd=1.5, bty="n", cex = cexs)
mtext(text = "Spinal BMD", 3, line = 0.4, adj = 0, cex=cex)
mtext(text = "Age", 1, line = 2.25, cex=cexs)
##
matplot(y=FFSCB_t_band,  x=grid, type="n", xlab="", ylab="")
polygon(x=c(grid,rev(grid)), 
        y=c(FFSCB_t_band[,2],rev(FFSCB_t_band[,1])), col = gray(.75), border = gray(.75))
abline(  h = 0, lwd=0.7)
lines(   y = hat_mu_diff,  x = grid, col=1, lty=1)
axis(4, at = 0, labels = expression(H[0]:~theta[f]-theta[m]==0))
legend(x=9.5, y=-0.075, legend = c(expression(paste("Estimated mean diff.")), expression(FF["frag,t"]^"2,2S")), 
       lty=c(1,1), bty="n", lwd = c(1.5,10), col=c("black", gray(.75)), cex =cexs, seg.len=2)
mtext(text = "Differences in mean functions", 3, line = 0.4, adj = 0, cex=cex)
mtext(text = "Age", 1, line = 2.25, cex=cexs)
dev.off()


width     <- 7.1
height    <- 3.7
par(mfrow=c(1,2), family = "serif", ps=13, cex.main=.99, font.main = 1, mar=c(4.1, 4.1, 2.1, 2.1))
image(cov_f, col=gray(seq(0,1,len=32)), xlab = "Age", ylab = "Age"); box()
mtext(text = "Female", 3, line = 0.4, adj = 0, cex=cex)
image(cov_m, col=gray(seq(0,1,len=32)), xlab = "Age", ylab = "Age"); box()
mtext(text = "Male", 3, line = 0.4, adj = 0, cex=cex)
dev.off()