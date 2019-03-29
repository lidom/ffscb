## ################################
## Make/Check package
remove.packages("ffscb")
devtools::document("ffscb")   
devtools::install_local("ffscb")
## #################################

# Load packages 
library("ffscb")

# library("here")
# here("Manuscript")

p         <- 201 
N         <- c(10,50,100)[2]
rangeval  <- c(0,1)
grid      <- make_grid(p, rangevals=rangeval)#, type == "close")
DGP       <- c("DGP1","DGP2","DGP3","DGP4","DGP5")[5]
##
## DGPs under the null-hypothesis
if(DGP=="DGP1"){
  mu        <- rep(0,p)
  mu0       <- mu
  cov.m     <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2, 1, 1))
  t0        <- grid[1]
}
if(DGP=="DGP2"){
  mu        <- rep(0,p)
  mu0       <- mu
  cov.m     <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(.75, 1, 1))
  t0        <- grid[1]
}
if(DGP=="DGP3"){
  mu        <- rep(0,p)
  mu0       <- mu
  cov.m     <- make_cov_m(cov.f = covf.st.matern.warp.power, grid=grid, cov.f.params=c(1.25, 1, 1, 2))
  t0        <- grid[p]
}
if(DGP=="DGP4"){
  mu        <- rep(0,p)
  mu0       <- mu
  cov.m     <- make_cov_m(cov.f = covf.st.matern.warp.sigmoid, grid=grid, cov.f.params=c(1.25, 1, 1))
  t0        <- grid[which(0.5==grid)]
}
## DGPs under the alternative-hypothesis
if(DGP=="DGP5"){
  mu        <- meanf_bump(x=grid, quarter = 4, height = .1) # plot(x=grid,y=mu)
  mu0       <- rep(0,p)
  cov.m     <- make_cov_m(cov.f = covf.st.matern.warp.power, grid=grid, cov.f.params=c(1.25, 1, 1, 2.5))
  t0        <- grid[p]
}
names(mu)  <- grid
names(mu0) <- grid

## check plot:
x  <-  make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
matplot(grid, x, type="l", lty=1); lines(grid, mu, lwd=2); confint(lm(x[1,]~1)) # naive.t
## 


set.seed(1110)
reps            <- 5000
type            <- c("naive.t", "Bs", "BEc", "KR.t", "FFSCB.t")
alpha.level     <- 0.10
n_int           <- 8
##
count_exceed    <- numeric(length(type)) 
count_exceed_t0 <- numeric(length(type)) 
crossings_loc   <- array(NA, dim = c(reps, p, length(type)))
# max_loc         <- matrix(NA, nrow=reps, ncol=length(type))
widths          <- numeric(length(type))
widths_sqr      <- numeric(length(type))
##
for(i in 1:reps){# 
  ## Generate data
  dat         <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
  ## Estimate mean, covariance, and tau
  hat_mu      <- rowMeans(dat)
  hat.cov.m   <- crossprod(t(dat - hat_mu)) / (N-1)
  hat.tau.v   <- tau_fun(dat)# plot(y=hat.tau.v,x=seq(0,1,len=p),type="l")
  ## Confidence bands
  b           <- confidence_band(x=hat_mu, cov=hat.cov.m, tau=hat.tau.v, t0=t0, N=N, type=type, 
                              conf.level=(1-alpha.level), n_int=n_int)# 
  # plot(b); abline(h=0)
  ##
  upper_Bands     <- b[,  2*(1:length(type))]
  lower_Bands     <- b[,1+2*(1:length(type))]
  ##
  tmp_up          <- upper_Bands < mu0    
  tmp_lo          <- lower_Bands > mu0
  # ## save locations of significant max_t(X(t)) locations:
  # tmp_min_up      <- apply(tmp_up,2,function(x){ifelse(length(hat_mu[x])>0,c(1:p)[x][which.min(hat_mu[x])],NA)})
  # tmp_max_lo      <- apply(tmp_lo,2,function(x){ifelse(length(hat_mu[x])>0,c(1:p)[x][which.max(hat_mu[x])],NA)})
  # max_loc[i,]     <- apply(cbind(tmp_min_up, tmp_max_lo), 1, 
  #                          function(x){ifelse(any(!is.na(x)), as.numeric(names(which.max(abs(hat_mu[x])))), NA)})
  ##
  ## counting of events: 'at least one crossing occured'
  tmp             <- tmp_up | tmp_lo
  count_exceed    <- count_exceed + as.numeric(apply(tmp, 2, function(x){any(x==TRUE)}))
  ## counting exceedances at t0:
  tmp_t0_up       <- upper_Bands[which(t0==grid),] < mu0[which(t0==grid)]      
  tmp_t0_lo       <- lower_Bands[which(t0==grid),] > mu0[which(t0==grid)]
  count_exceed_t0 <- count_exceed_t0 + as.numeric(tmp_t0_up | tmp_t0_lo)
  ## saving band-crossing locations:
  for(j in 1:length(type)){
    ## delta: number of gridpoints before t0:
    if(which(grid==t0) != p){delta <- (which(grid==t0)-1)}else{delta <- 0}
    ## crossing locations with respect to the upper Bands:
    tmp_cr_up  <- c(      locate_crossings(mu0[1:which(grid==t0)],upper_Bands[1:which(grid==t0),j],type="down"), # down-crossings (left of t0)
                    delta+locate_crossings(mu0[which(grid==t0):p],upper_Bands[which(grid==t0):p,j],type="up"  )) #   up-crossings (right of t0)
    ## crossing locations with respect to the lower Bands:
    tmp_cr_lo  <- c(      locate_crossings(mu0[1:which(grid==t0)],lower_Bands[1:which(grid==t0),j],type="up"  ), # down-crossings (left of t0)
                    delta+locate_crossings(mu0[which(grid==t0):p],lower_Bands[which(grid==t0):p,j],type="down")) #   up-crossings (right of t0)
    ## all (sorted) crossing locations together:
    tmp_cr_loc <- sort(c(tmp_cr_up, tmp_cr_lo))
    ## save crossing locations (if any):
    if(length(tmp_cr_loc)>0){crossings_loc[i,tmp_cr_loc,j]  <- grid[tmp_cr_loc]}
  }
  ## saving band-widths:
  widths      <- widths     + colMeans( upper_Bands - lower_Bands)
  widths_sqr  <- widths_sqr + colMeans((upper_Bands - lower_Bands)^2)
  ##
  if(i %% 500 ==0){cat("i /",reps,"=",i,"/",reps,"\n")}
  ##
} 
##
save(type, 
     reps,
     count_exceed,
     count_exceed_t0,
     crossings_loc,
     widths,
     widths_sqr, file = paste0(here("R/Simulation_Results/"),DGP,"_N=",N))


exceed_frq           <- count_exceed   /reps
names(exceed_frq)    <- names(widths)
exceed_t0_frq        <- count_exceed_t0/reps
names(exceed_t0_frq) <- names(widths)

exceed_frq




library("dplyr")
type_order <- pmatch(type, names(widths))
type       <- type[type_order]
crossings_df <- dplyr::tibble(
  !!type[1] :=  c(crossings_loc[,,1]),
  !!type[2] :=  c(crossings_loc[,,2]),
  !!type[3] :=  c(crossings_loc[,,3]),
  !!type[4] :=  c(crossings_loc[,,4]),
  !!type[5] :=  c(crossings_loc[,,5])
) %>% 
  tidyr::gather(band, cr_loc, factor_key=TRUE) %>% 
  dplyr::mutate(intervals = as.factor(findInterval(x   = cr_loc, 
                                                   vec = seq(0,1,len=5), 
                                                   rightmost.closed = T)))


crossings_df %>% 
  tidyr::drop_na() %>% 
  dplyr::group_by(band) %>% 
  dplyr::mutate(n_cr = n()) %>% 
  dplyr::group_by(band, intervals) %>% 
  dplyr::summarise(rfrq = n()/n_cr[1]) %>% 
  mutate_if(is.numeric, round, 3) %>% 
  print(n=Inf)





# crossings_df %>% 
#   tidyr::drop_na() %>% 
#   dplyr::group_by(band, intervals) %>% 
#   dplyr::summarise(rfrq = n()/reps) %>% 
#   mutate_if(is.numeric, round, 3) %>% 
#   print(n=Inf)




plot(density(c(na.omit(c(crossings_loc[,,4]))), from = 0, to=1), ylim=c(0,1))
points(y=rep(0,len=length(c(na.omit(c(crossings_loc[,,4]))))), 
       x=c(na.omit(c(crossings_loc[,,4]))))



stem(c(na.omit(c(crossings_loc[,,4]))), scale = .2, width = 20)








## -------------------------------------------



















##
widths
##
plot(b)
##
points(grid,mu,typ="l",lty=2)
##
# data.frame(type          = type, 
#            coverage      = coverage, 
#            me            = round(2*sqrt(coverage*(1-coverage)/reps),4), 
#            ave_width     = widths/reps,  
#            ave_sqr_width = widths_sqr/reps)
##
plot_data <- data.frame(grid=grid, counts_grid/reps)
names(plot_data)[-1] <- type
plot_data <- melt(plot_data, id.var="grid", value.name="coverage")
##
ggplot(plot_data,aes(x=grid,y=coverage,group=variable,color=variable))+
  geom_line(aes(lty=variable)) +
  ylab("Pointwise Coverage") + xlab("") + theme_bw() 


dat         <-  make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
hat_mu      <- rowMeans(dat)
hat.cov.m   <- crossprod(t(dat - hat_mu)) / (N-1)
hat.tau.v   <- tau_fun(dat)
b           <- confidence_band(x = hat_mu, cov = hat.cov.m, tau=hat.tau.v, N=N, type=type, conf.level = 0.95, 
                            n_int = 10, tol=.Machine$double.eps^0.25)
par(mfrow=c(2,1), mar=c(2.1, 4.1, 2.1, 2.1))
matplot(grid, dat, type="l", lty=1, ylab="X(t)",xlab = "")
plot(b, xlab="", ylab="Bands")
par(mfrow=c(2,1), mar=c(5.1, 4.1, 4.1, 2.1))

# matplot(x=grid,y=cbind(hat_mu, b[,6:7]),type="l",lty=c(1,2,2), col=1)
# matlines(x=grid,y=b[,4:5], lty=3, col=2)


