library("ffscb")

## Data (selecting 'extra' and 'normal cushioned')
slct_EC      <- grepl(pattern = "Extra_Cush_",  x = colnames(Biomechanics))
slct_NC      <- grepl(pattern = "Normal_Cush_", x = colnames(Biomechanics))
##
grid         <- Biomechanics[,1]/100   # Grid in [0,1]  
EC_mat       <- Biomechanics[,slct_EC] # Torque curves when running with extra  cushioned running shoes
NC_mat       <- Biomechanics[,slct_NC] # Torque curves when running with normal cushioned running shoes
Dff_mat      <- EC_mat - NC_mat        # Difference curves
N            <- ncol(Dff_mat)          # Sample size

## Computing the estimates
hat_mu       <- rowMeans(Dff_mat)
hat.cov      <- crossprod(t(Dff_mat - hat_mu)) / (N-1)
hat.cov.mu   <- hat.cov / N
hat.tau      <- tau_fun(Dff_mat) 
plot(y=hat.tau,x=grid,type="l") # Roughness parameter function tau(t)

## Selecting bands, t0, n_int, and significance level
type         <- c("KR.t","FFSCB.t")
t0           <- 0.2
n_int        <- 5
alpha.level  <- 0.05

## Computing the confidence bands: FF_t and KR_t
b            <- confidence_band(x=hat_mu, cov=hat.cov.mu, tau=hat.tau, t0=0, df=N-1,
                                    type=type, conf.level=(1-alpha.level), n_int=n_int)
plot(b); abline(h=0)

slct_KR_t    <- grepl(pattern = "KR.t",    colnames(b))
slct_FFSCB_t <- grepl(pattern = "FFSCB.t", colnames(b))

KR_t_band    <- b[,slct_KR_t] 
FFSCB_t_band <- b[,slct_FFSCB_t] #FFSCB_res$band[,-1]

KR_threshold_u    <- (KR_t_band[,1]    - hat_mu)/sqrt(diag(hat.cov.mu))
FFSCB_threshold_u <- (FFSCB_t_band[,1] - hat_mu)/sqrt(diag(hat.cov.mu))

## 1st Significant Area:
## FF-t
signif_int1_FF <- (FFSCB_t_band[,1] - 0) < 0
range(grid[signif_int1_FF])*100
## KR-t
signif_int1_KR <- (KR_t_band[,1] - 0) < 0
range(grid[signif_int1_KR])*100

## 2nd Significant Area:
## FF-t
signif_int2_FF <- (FFSCB_t_band[,2] - 0) > 0
range(grid[signif_int2_FF])*100
## KR-t (no signif. area)
signif_int2_KR <- (KR_t_band[,2] - 0) > 0
range(grid[signif_int2_KR])*100

## Evaluations points for computing point-wise p-values:
ep1 <- grid[which.min(FFSCB_t_band[,1] - 0)] 
ep2 <- grid[which.max(FFSCB_t_band[,2] - 0)]


## Evaluation point at which to compute the pvalue
ep_vec       <- c(ep1,ep2)
pvalues      <- get_pvalue_FFSCB_t(x=hat_mu, tau=hat.tau, t0=t0, diag.cov=diag(hat.cov.mu), df=N-1,
                                   eval.points = ep_vec, n_int = n_int)
round(pvalues, 3)

## #############
## Plots      ##
## #############
width     <- 7
height    <- 2.5
mar       <- c(4.1, 4.1, 4.1, 1.1)
##
layout(mat = matrix(c(1:3), nrow=1, ncol=3),
       widths = c(1, 1, 1))  # Widths of the two columns
##
par(family = "serif", ps=13, cex.main=1.3, cex.lab=1.25, cex.axis=1.3, font.main = 1, mar=mar)
matplot(y = EC_mat,  x = grid*100, lwd=.5, col=gray(.25), type="l", lty=1, 
        ylab="N m / kg", xlab = "% of Stance Phase", main="Extra Cushioned\nRunning Shoe")
matplot(y = NC_mat,  x = grid*100, lwd=.5, col=gray(.25), type="l", lty=1, 
        ylab="N m / kg", xlab = "% of Stance Phase", main="Normal Cushioned\nRunning Shoe")
matplot(y = Dff_mat, x = grid*100, lwd=.5, col=gray(.5), type="l", lty=1, ylim = c(-0.4, max(Dff_mat)), 
        ylab="N m / kg", xlab = "% of Stance Phase", main="Difference Curves\n(Extra - Normal)")
lines(y = hat_mu, x = grid*100)
legend("bottomright", legend = expression(paste("Estimated mean")), col="black", lty=1, bty="n", lwd = 1.5, cex=1.3)
par(mfrow=c(1,1))
dev.off()

width     <- 7
height    <- 3.125
cex       <- 1
cexs      <- 0.95
##
layout(mat = matrix(c(1,2,3,3), nrow=2, ncol=2),
       heights = c(1, 1),  # Heights of the two rows
       widths =  c(1, 1.25))  # Widths of the two columns

par(family = "serif", ps=13, cex.axis=1.05, font.main = 1)
par(mar=c(2.5, 2.1, 2, 1.1))
plot(y=hat.tau,x=grid*100, type="l", main = "", xlab = "", ylab="")
mtext(text = expression(paste("Roughness ",tau)), 3, line = 0.4, adj = 0, cex=cex)
##
par(mar=c(3.1, 2.1, 1.5, 1.1))
matplot(y=cbind(KR_threshold_u, FFSCB_threshold_u), x=grid*100, col="black", lty=c(2,1), type = "l", 
        main="", xlab = "", ylab="", ylim=c(2.7, 4.2), lwd=c(1.5,1))
legend("topright", 
       legend = c(expression(KR[t]), expression(FF[t])), 
       lty=c(2,1), bty="n", lwd = c(1.5,1.5), col=c("black", "black"), cex = cex, seg.len=2)
mtext(text = expression(paste("Thresholds ",u)), 3, line = 0.4, adj = 0, cex=cex)
mtext(text = "% of Stance Phase", 1, line = 2.25, cex=cexs)
## 
par(mar=c(3.1, 4.1, 2, 2.1))
matplot( y = b, x = grid*100, type="n", ylim=c(min(FFSCB_t_band),0.28), 
         ylab="", xlab = "", main="")
polygon(x=c(grid*100,rev(grid*100)), y=c(FFSCB_t_band[,2],rev(FFSCB_t_band[,1])), col = gray(.75), border = gray(.75))
matlines(x=grid*100, y=KR_t_band, lty=2, col=1)
abline(  h = 0, lwd=0.7)
lines(   y = hat_mu,  x = grid*100, col=1, lty=1)
axis(4, at = 0, labels = expression(H[0]:~theta==0))
legend("topright", legend = c(expression(paste("Estimated mean")), 
                              expression(KR[t]), expression(FF[t])), 
      lty=c(1,2,1), bty="n", lwd = c(1.5,1,10), col=c("black", "black",gray(.75)), cex =cex, seg.len=2)
mtext(text = "% of Stance Phase", 1, line = 2.25, cex=cexs)
mtext(text = "N m / kg", 2, line = 2.5, cex=cex)
mtext(text = "Simultaneous 95% Confidence Bands", 3, line = 0.4, cex=cex)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()


## ####################################
## Region of Interest (ROI) analysis ##
## ####################################

## Computing the FF_t confidence band
FFSCB_res <- make_band_FFSCB_t(x=hat_mu, diag.cov.x = diag(hat.cov.mu), tau=hat.tau, t0=t0, df=N-1, 
                               conf.level=(1-alpha.level), n_int=n_int)

## Results
FFSCB_res$t0
FFSCB_res$prob_t0
FFSCB_res$a_star 

## Confidence level for the region of interest ROI=[0,20%] :
(1 - ( round(FFSCB_res$prob_t0,3) + round(FFSCB_res$a_star,3)*(FFSCB_res$t0 - 0))  ) * 100 