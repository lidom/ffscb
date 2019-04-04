## ################################
## Make/Check package
remove.packages("ffscb")
devtools::document("ffscb")   
devtools::install_local("ffscb")
## #################################

# Load packages 
library("ffscb")
help("ffscb")

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
  # mu        <- meanf_bump(x=grid, quarter = 4, height = .1) # plot(x=grid,y=mu)
  # mu0       <- rep(0,p)
  mu        <- rep(0.1,p)
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
count_exceed     <- numeric(length(type)) 
count_exceed_t0  <- numeric(length(type)) 
crossings_loc    <- array(NA, dim = c(reps, p, length(type)))
exceed_up_loc    <- array(NA, dim = c(reps, p, length(type)))
exceed_lo_loc    <- array(NA, dim = c(reps, p, length(type)))
intgr_widths     <- numeric(length(type))
intgr_widths_sqr <- numeric(length(type))
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
  upper_Bands <- b[,  2*(1:length(type))]
  lower_Bands <- b[,1+2*(1:length(type))]
  ##
  tmp_up      <- upper_Bands < mu0    
  tmp_lo      <- lower_Bands > mu0
  ##
  exceed_up_loc[i,,] <- tmp_up
  exceed_lo_loc[i,,] <- tmp_lo
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
  intgr_widths      <- intgr_widths     + colSums( upper_Bands - lower_Bands   )*diff(grid)[1]
  intgr_widths_sqr  <- intgr_widths_sqr + colSums((upper_Bands - lower_Bands)^2)*diff(grid)[1]
  ##
  if(i %% 500 ==0){cat("i /",reps,"=",i,"/",reps,"\n")}
  ##
} 
##


CI_type      <- stringr::str_remove(names(widths), paste0(".u.",(1-alpha.level)))
crossings_df <- dplyr::tibble(
  !!CI_type[1] :=  c(crossings_loc[,,1]),
  !!CI_type[2] :=  c(crossings_loc[,,2]),
  !!CI_type[3] :=  c(crossings_loc[,,3]),
  !!CI_type[4] :=  c(crossings_loc[,,4]),
  !!CI_type[5] :=  c(crossings_loc[,,5])
) %>% 
  tidyr::gather(band, cr_loc, factor_key=TRUE) %>% 
  dplyr::mutate(intervals = as.factor(findInterval(x   = cr_loc, 
                                                   vec = seq(0,1,len=5), 
                                                   rightmost.closed = T)))


rfrq_interv_df <- crossings_df %>% 
  tidyr::drop_na() %>% 
  dplyr::group_by(band) %>% 
  dplyr::mutate(n_cr = n()) %>% 
  dplyr::group_by(band, intervals) %>% 
  dplyr::summarise(rfrq = n()/n_cr[1]) %>% 
  mutate_if(is.numeric, round, 3) 


KL_df <- rfrq_interv_df %>% 
  dplyr::group_by(band) %>% 
  dplyr::summarise(KL_rfrq_interv=sum(rfrq*log(rfrq/rep(.25,4)))) 


sim_df <- dplyr::tibble("DGP"              = DGP,
                        "N"                = N,
                        "Band"             = CI_type,
                        "alpha.level"      = alpha.level,
                        "exceed_frq"       = count_exceed   /reps,
                        "exceed_t0_frq"    = count_exceed_t0/reps, 
                        "intgr_widths"     = intgr_widths,
                        "intgr_widths_sqr" = intgr_widths_sqr, 
                        "KL"               = round(KL_df$KL_rfrq_interv,2), 
                        "rfrq_1"           = rfrq_interv_df %>% dplyr::filter(intervals==1),
                        "rfrq_2"           = rfrq_interv_df %>% dplyr::filter(intervals==2),
                        "rfrq_3"           = rfrq_interv_df %>% dplyr::filter(intervals==3),
                        "rfrq_4"           = rfrq_interv_df %>% dplyr::filter(intervals==4))

save(sim_df, file = paste0(here("R/Simulation_Results/"), DGP,"_N=", N))



# ffscb_m <- apply(exceed_lo_loc[,,5], 1, function(x){tmp <- rep(NA,p);tmp[x] <- grid[x]; tmp})
# boots_m <- apply(exceed_lo_loc[,,1], 1, function(x){tmp <- rep(NA,p);tmp[x] <- grid[x]; tmp})
# 
# ffscb_v <- apply(ffscb_m, 2, function(x){ifelse(any(!is.na(x)), max(x,na.rm=T)-min(x,na.rm=T),NA)})
# boot_v  <- apply(boots_m, 2, function(x){ifelse(any(!is.na(x)), max(x,na.rm=T)-min(x,na.rm=T),NA)})
# 
# ffscb_v <- c(na.omit(ffscb_v))
# boot_v  <- c(na.omit(boot_v))
# 
# summary(ffscb_v)
# summary(boot_v)



