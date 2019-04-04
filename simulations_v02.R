## ################################
## Make/Check package
remove.packages("ffscb")
devtools::document("ffscb")   
devtools::install_local("ffscb")
## #################################

# Load packages 
# devtools::install_github("alessiapini/fdatest")
library("tidyverse")
library("fdatest")
library("ffscb")
help("ffscb")


##
p            <- 101
rangeval     <- c(0,1)
grid         <- make_grid(p, rangevals=rangeval)
type         <- c("naive.t", "Bs", "BEc", "KR.t", "FFSCB.t")
alpha.level  <- 0.10
n_int        <- 8
##
reps         <- 500
##
DGP_seq   <- c("DGP1_H0","DGP1_H1", 
               "DGP2_H0","DGP2_H1", 
               "DGP3_H0","DGP3_H1",
               "DGP4_H0","DGP4_H1")
N_seq     <- c(10,50,100)

## DGP <- DGP_seq[2]; N <- 100

##
for(DGP in DGP_seq) {
  ##
  set.seed(1110)
  ##
  for(N in N_seq) {
    ## 
    if( any(DGP==c("DGP1_H0","DGP1_H1")) ) {# smooth
      mu        <- meanf_poly(grid, c(1,2))
      if(grepl("H0", DGP)) { mu0 <- mu } else { mu0 <- mu - 0.1}
      cov.m     <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2, 1, 1))
      t0        <- grid[1]
    }
    if( any(DGP==c("DGP2_H0","DGP2_H1")) ) {# rough
      mu        <- meanf_poly(grid, c(1,2))
      if(grepl("H0", DGP)) { mu0 <- mu } else { mu0 <- mu - 0.1}
      cov.m     <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(.75, 1, 1))
      t0        <- grid[1]
    }
    if( any(DGP==c("DGP3_H0","DGP3_H1")) ) {# from smooth to rough
      mu        <- meanf_poly(grid, c(1,2))
      if(grepl("H0", DGP)) { mu0 <- mu } else { mu0 <- mu - 0.1}
      cov.m     <- make_cov_m(cov.f = covf.st.matern.warp.power, grid=grid, cov.f.params=c(1.25, 1, 1, 2.5))
      t0        <- grid[p]
    }
    if( any(DGP==c("DGP4_H0","DGP4_H1")) ) {# from smooth to rough to smooth
      mu        <- meanf_poly(grid, c(1,-1))
      if(grepl("H0", DGP)) { mu0 <- mu } else { mu0 <- mu - 0.1}
      cov.m     <- make_cov_m(cov.f = covf.st.matern.warp.sigmoid, grid=grid, cov.f.params=c(1.25, 1, 1))
      t0        <- grid[which(0.5==grid)]
    }
    names(mu)  <- grid
    names(mu0) <- grid
    
    ## check plot:
    #x  <-  make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
    #matplot(grid, x, type="l", lty=1); lines(grid, mu, lwd=2); confint(lm(x[1,]~1)) # naive.t
    ## 
    
    count_exceed      <- numeric(length(type)) 
    count_exceed_IWT  <- 0
    count_exceed_t0   <- numeric(length(type)) 
    crossings_loc     <- array(NA, dim = c(reps, p, length(type)))
    crossings_loc_IWT <- matrix(NA, nrow = reps, ncol = p)
    exceed_loc        <- array(NA, dim = c(reps, p, length(type)))
    exceed_loc_IWT    <- matrix(NA, nrow = reps, ncol = p)
    intgr_widths      <- numeric(length(type))
    intgr_widths_sqr  <- numeric(length(type))
    ##
    for(i in 1:reps){# 
      ## Generate data
      dat         <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
      ## Estimate mean, covariance, and tau
      hat_mu      <- rowMeans(dat)
      hat.cov.m   <- crossprod(t(dat - hat_mu)) / (N-1)
      hat.tau.v   <- tau_fun(dat)# plot(y=hat.tau.v,x=seq(0,1,len=p),type="l")
      ## Confidence bands
      b           <- confidence_band(x=hat_mu, cov=hat.cov.m, tau=hat.tau.v, t0=t0, N=N, type=type, conf.level=(1-alpha.level), n_int=n_int)
      # plot(b); abline(h=0)
      ## ==============================================================================================
      ## IWT1 function from the fdatest package
      IWT_messages       <- capture.output(IWT_res <- fdatest::IWT1(data = t(dat), mu = mu0))
      count_exceed_IWT   <- count_exceed_IWT + as.numeric(any(IWT_res$adjusted_pval < alpha.level))
      exceed_loc_IWT[i,] <- IWT_res$adjusted_pval < alpha.level
      tmp_cr_loc_IWT     <- locate_crossings(IWT_res$adjusted_pval, rep(alpha.level, p), type="down")# down-crossing of pval means upcrossing of empirical mean
      crossings_loc_IWT[i,tmp_cr_loc_IWT]  <- grid[tmp_cr_loc_IWT]
      
      get_pvalue_FFSCB_t(x=hat_mu, x0=NULL, tau=hat.tau.v, diag.cov = diag(hat.cov.m), t0 = t0, N = N)
      
      ## ==============================================================================================
      ##
      upper_Bands  <- b[,  2*(1:length(type))]
      lower_Bands  <- b[,1+2*(1:length(type))]
      ##
      exceed_loc[i,,] <- upper_Bands < matrix(mu0, nrow=p, ncol=length(type)) | lower_Bands > matrix(mu0, nrow=p, ncol=length(type))
      ##
      ## counting of events: 'at least one crossing occured'
      tmp             <- exceed_loc[i,,]
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
        tmp_cr_up  <- c(locate_crossings(mu0[1:which(grid==t0)],upper_Bands[1:which(grid==t0),j],type="down"), # down-crossings (left of t0)
                        delta+locate_crossings(mu0[which(grid==t0):p],upper_Bands[which(grid==t0):p,j],type="up"  )) #   up-crossings (right of t0)
        ## crossing locations with respect to the lower Bands:
        tmp_cr_lo  <- c(locate_crossings(mu0[1:which(grid==t0)],lower_Bands[1:which(grid==t0),j],type="up"  ), # down-crossings (left of t0)
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
    
    
    CI_type      <- c(stringr::str_replace(names(intgr_widths), paste0(".u.",(1-alpha.level)),""), "IWT")
    crossings_df <- dplyr::tibble(
      !!CI_type[1] :=  c(crossings_loc[,,1]),
      !!CI_type[2] :=  c(crossings_loc[,,2]),
      !!CI_type[3] :=  c(crossings_loc[,,3]),
      !!CI_type[4] :=  c(crossings_loc[,,4]),
      !!CI_type[5] :=  c(crossings_loc[,,5]),
      "IWT"         =  c(crossings_loc_IWT)
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
      dplyr::summarise(rfrq = n()/n_cr[1]) 
    
    
    KL_df <- rfrq_interv_df %>% 
      dplyr::group_by(band) %>% 
      dplyr::summarise(KL_rfrq_interv=sum(rfrq*log(rfrq/rep(.25,4)))) %>% 
      dplyr::ungroup() 
    
    
    # slct_comms <- c(apply(exceed_loc[,,1], 1, function(x){any(x==TRUE)})&
    #   apply(exceed_loc[,,2], 1, function(x){any(x==TRUE)})&
    #   apply(exceed_loc[,,3], 1, function(x){any(x==TRUE)})&
    #   apply(exceed_loc[,,4], 1, function(x){any(x==TRUE)})&
    #   apply(exceed_loc[,,5], 1, function(x){any(x==TRUE)})&
    #   apply(exceed_loc_IWT, 1, function(x){any(x==TRUE)}))
    # 
    # 
    # ## For all DGPs under the alternative: Compute the lengths of the exceedances
    # if( grepl("H1", DGP) ) {
    #   detected_H1_len_med <- NULL
    #   ##
    #   target_loc <- mu!=mu0
    #   ##
    #   for( i in 1:(length(CI_type)-1) ) {
    #     detected_H1_len_med <- c(detected_H1_len_med,
    #                              apply(exceed_loc[slct_comms,,i], 1, function(x){
    #                                if(any(x==TRUE)) {
    #                                  tmp <- x == target_loc
    #                                  length(tmp[tmp==TRUE])/length(target_loc[target_loc==TRUE])
    #                                } else {
    #                                  NA  
    #                                }}) %>% median(.,na.rm = T)) 
    #   }
    #   detected_H1_len_med <- c(detected_H1_len_med,
    #                            apply(exceed_loc_IWT[slct_comms,], 1, function(x){
    #                              if(any(x==TRUE)) {
    #                                tmp <- x == target_loc
    #                                length(tmp[tmp==TRUE])/length(target_loc[target_loc==TRUE])
    #                              } else {
    #                                NA  
    #                              }}) %>% median(.,na.rm = T)) 
    # } else {
    #   detected_H1_len_med <- rep(NA,length(CI_type))
    # }
    
    sim_df <- dplyr::tibble("DGP"              = DGP,
                            "N"                = N,
                            "Band"             = CI_type,
                            "alpha.level"      = alpha.level,
                            "exceed_frq"       = c(count_exceed, count_exceed_IWT) /reps,
                            "KL"               = KL_df$KL_rfrq_interv, 
                            "rfrq_I1"          = rfrq_interv_df %>% dplyr::filter(intervals==1) %>% pull(rfrq),
                            "rfrq_I2"          = rfrq_interv_df %>% dplyr::filter(intervals==2) %>% pull(rfrq),
                            "rfrq_I3"          = rfrq_interv_df %>% dplyr::filter(intervals==3) %>% pull(rfrq),
                            "rfrq_I4"          = rfrq_interv_df %>% dplyr::filter(intervals==4) %>% pull(rfrq),
                            ##
                            #"detected_H1_len_med"  = detected_H1_len_med,
                            ##
                            "exceed_t0_frq"    = c(count_exceed_t0,NA)/reps, 
                            "intgr_widths"     = c(intgr_widths,NA),
                            "intgr_widths_sqr" = c(intgr_widths_sqr,NA))
    
    save(sim_df, file = paste0("Simulation_Results/", DGP,"_N=", N))
    
  }
}




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



