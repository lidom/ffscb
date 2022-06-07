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

## Computing the estimates
hat_mu       <- rowMeans(Dff_mat)
hat.cov      <- crossprod(t(Dff_mat - hat_mu)) / (N-1)
hat.cov.mu   <- hat.cov / N
hat.tau      <- tau_fun(Dff_mat) 

## 
alpha.level  <- 0.05
n_int        <- 8

## Over each subinterval [0, 12.5%], [12.5%, 25%], ... , [87.5%, 100%]  
## we compute a 99.375% simultanouse confidence interval:
(1 - alpha.level*(1/n_int)) * 100 

## Computing the confidence bands
b                 <- confidence_band(x          = hat_mu, 
                                     cov        = hat.cov.mu, 
                                     tau        = hat.tau, 
                                     df         = N-1, 
                                     type       ="FFSCB.t", 
                                     conf.level = 1-alpha.level, 
                                     n_int      = n_int)
FF8_t_band        <- b[,-1] 
FF8_crit_value_u  <- (FF8_t_band[,1] - hat_mu)/sqrt(diag(hat.cov.mu))

## Significant Areas:
range(grid[(FF8_t_band[,1] - 0) < 0])*100
range(grid[(FF8_t_band[,2] - 0) > 0])*100


## Plots
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
lines(y = hat_mu, x = grid*100)
legend("bottomright", legend = expression(paste("Estimated mean")), col="black", lty=1, bty="n", lwd = 1.5, cex=1.3)
par(mfrow=c(1,1))
dev.off()

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
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()