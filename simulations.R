# Load packages 
# devtools::install_github("alessiapini/fdatest")
library("tidyverse")
library("fda")
library("fdatest")
library("parallel")
library("ffscb")

##
detectCores()
nworkers <- 6

## Setup ####################################################
p            <- 101
grid         <- make_grid(p, rangevals=c(0,1))
type         <- c("naive.t", "Bs", "BEc", "KR.t", "FFSCB.t")
alpha.level  <- 0.05
n_int        <- 10
##
n_reps_H0    <- 50000
n_reps_H1    <- 10000
##
DGP_seq      <- c("DGP1_shift","DGP1_scale","DGP1_local",
                  "DGP2_shift","DGP2_scale","DGP2_local", 
                  "DGP3_shift","DGP3_scale","DGP3_local",
                  "DGP4_shift","DGP4_scale","DGP4_local")
##
delta_Nsmall  <- c(0, seq(from = 0.04, to = 0.2, len = 5))
delta_Nlarge  <- c(0, seq(from = 0.02, to = 0.1, len = 5))
##
N_seq         <- c(10,100)
## #########################################################

##
for(DGP in DGP_seq) {
  ##
  set.seed(1110)
  ##
  for(N in N_seq) {
    ## Take the correct delta_seq corresponding to N
    if ( N==min(N_seq) ) delta_seq <- delta_Nsmall else delta_seq <- delta_Nlarge
    ##
    for(delta in delta_seq) {# DGP <- "DGP1_shift"; N <- 100; delta <- 0.1
      ## 
      if(grepl("shift", DGP)) { mu0 <- meanf_shift(grid, 0);      mu <- meanf_shift(grid, delta) }
      if(grepl("scale", DGP)) { mu0 <- meanf_scale(grid, 0);      mu <- meanf_scale(grid, delta) }
      if(grepl("local", DGP)) { mu0 <- meanf_localshift(grid, 0); mu <- meanf_localshift(grid, delta) }
      names(mu)  <- grid
      names(mu0) <- grid
      ##
      if(grepl("DGP1", DGP)) {# stationary: smooth 
        cov.m     <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(3/2, 1, 1/4))
        t0        <- grid[1]
      }
      if(grepl("DGP2", DGP)) {# stationary: rough
        cov.m     <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(1/2, 1, 1/4))
        t0        <- grid[1]
      }
      if(grepl("DGP3", DGP)) {# non-stationary: from smooth to rough
        cov.m     <- make_cov_m(cov.f = covf.st.matern.warp.power, grid=grid, cov.f.params=c(1.25, 1, 1/4, 2.5))
        t0        <- grid[p]
      }
      if(grepl("DGP4", DGP)) {# non-stationary: from smooth to rough to smooth
        cov.m     <- make_cov_m(cov.f = covf.st.matern.warp.sigmoid, grid=grid, cov.f.params=c(1.25, 1, 1/4))
        t0        <- grid[which(0.5==grid)]
      }
      ## check plot:
      # sim.dat  <-  make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
      # matplot(grid, sim.dat, type="l", lty=1); lines(grid, mu, lwd=2)
      ## 
      ## Number of Monte-Carlo repetitions
      n_reps     <- ifelse(delta==0, n_reps_H0, n_reps_H1)
      ##
      start_time <- Sys.time()
      ##
      res_mclapply <- mclapply(1:n_reps, function(reps) {
        ## Generate data
        dat         <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
        ## Estimate mean, covariance, and tau
        hat_mu      <- rowMeans(dat)
        hat.cov.m   <- crossprod(t(dat - hat_mu)) / (N-1)
        hat.tau.v   <- tau_fun(dat) # plot(y=hat.tau.v,x=seq(0,1,len=p),type="l")
        ##
        ## Confidence bands
        b               <- confidence_band(x=hat_mu, cov=hat.cov.m, tau=hat.tau.v, t0=t0, N=N, 
                                           type=type, conf.level=(1-alpha.level), n_int=n_int)# 
        # plot(b); lines(x=grid, y=mu, col="blue"); lines(x=grid, y=mu0, lty=2, col="blue")
        ##
        upper_Bands     <- b[,  2*(1:length(type))]
        lower_Bands     <- b[,1+2*(1:length(type))]
        ##
        exceed_loc      <- upper_Bands < matrix(mu0, nrow=p, ncol=length(type)) | lower_Bands > matrix(mu0, nrow=p, ncol=length(type))
        ##
        ## saving exceedances events ('at least one crossing occured?')
        exceedances     <- as.numeric(apply(exceed_loc, 2, function(x){any(x==TRUE)}))
        ##
        ## saving exceedances at t0:
        tmp_t0_up       <- upper_Bands[which(t0==grid),] < mu0[which(t0==grid)]      
        tmp_t0_lo       <- lower_Bands[which(t0==grid),] > mu0[which(t0==grid)]
        exceedances_t0  <- as.numeric(tmp_t0_up | tmp_t0_lo)
        ##
        ## saving band x mu0 crossing locations:
        crossings_loc <- matrix(NaN, nrow=p, ncol=length(type))
        for(j in 1:length(type)){
          if(all(exceed_loc[,j]==TRUE )){ crossings_loc[,j] <- rep( 100,p) }# no crossing since mu0 is completely outside of the band
          if(all(exceed_loc[,j]==FALSE)){ crossings_loc[,j] <- rep(-100,p) }# no crossing since mu0 is completely inside  of the band
          if(any(exceed_loc[,j]==TRUE) & any(exceed_loc[,j]==FALSE)){   
            ## 'ngbt0': number of gridpoints before t0:
            if(which(grid==t0) != p){ ngbt0 <- (which(grid==t0)-1) }else{ ngbt0 <- 0 }
            ## crossing locations with respect to the upper Bands:
            tmp_cr_up  <- c(locate_crossings(mu0[1:which(grid==t0)],upper_Bands[1:which(grid==t0),j],type="down"),         # down-crossings (left of t0)
                            ngbt0 + locate_crossings(mu0[which(grid==t0):p],upper_Bands[which(grid==t0):p,j],type="up"  )) #   up-crossings (right of t0)
            ## crossing locations with respect to the lower Bands:
            tmp_cr_lo  <- c(locate_crossings(mu0[1:which(grid==t0)],lower_Bands[1:which(grid==t0),j],type="up"  ),         # down-crossings (left of t0)
                            ngbt0 + locate_crossings(mu0[which(grid==t0):p],lower_Bands[which(grid==t0):p,j],type="down")) #   up-crossings (right of t0)
            ## all (sorted) crossing locations together:
            tmp_cr_loc <- sort(c(tmp_cr_up, tmp_cr_lo))
            ## save crossing locations (if any):
            if( length(tmp_cr_loc)>0 ){ crossings_loc[tmp_cr_loc,j]  <- grid[tmp_cr_loc] }
          }
        }
        #apply(crossings_loc,2,function(x)any(!is.na(x)))
        ## saving widths of the bands:
        intgr_widths_sqr  <- colSums((upper_Bands - lower_Bands)^2)*diff(grid)[1]
        ## Names of methods in correct order:
        Band_type         <- stringr::str_replace(string  = names(intgr_widths_sqr), 
                                                  pattern = paste0(".u.", (1-alpha.level)),"")
        ## simulation data
        sim_df <- dplyr::tibble(band     = as_factor(rep(Band_type, each = p)),
                                excd     = rep(exceedances,         each = p),
                                excd_t0  = rep(exceedances_t0,      each = p),
                                excd_loc = c(exceed_loc),
                                cros_loc = c(crossings_loc),
                                wdth     = rep(intgr_widths_sqr,    each = p))# glimpse(sim_df)
        ##
        return(sim_df)
      }, mc.cores = nworkers)
      ##
      end_time <- Sys.time()
      run_time <- end_time - start_time
      ##
      ## Combine all simulation results (rowbinding the n_reps many tibbles in 'res_mclapply')
      sim_df <- dplyr::bind_rows(res_mclapply) %>% 
        dplyr::mutate(
          intervals = cut(x=cros_loc, breaks=c(-101,seq(0,1,len=5),101), labels = c('compl_in',1,2,3,4,'compl_out')),# Assign crossing-locations to one of four equidistant intervals in [0,1]
          cros_loc  = replace(cros_loc, cros_loc == -100 | cros_loc == 100, NA),# Setting the 'compl_in/out'  indicators to NA
          run       = rep(c(1:n_reps), each=length(type)*p),   # Numbering the single simulation runs
          n_rep     = rep(n_reps, times=length(type)*p*n_reps),# Total number of simulation runs  
          delta     = rep(delta,  times=length(type)*p*n_reps),# Delta (from H0 to H1)  
          N         = rep(N,      times=length(type)*p*n_reps),# Sample size
          DGP       = rep(DGP,    times=length(type)*p*n_reps))# Name of DGP
      ##
      ## Feedback
      cat(DGP, ", N=", N, ", Delta=", delta, ", Run-Time=", run_time, " (", attr(run_time, "units"),")\n", sep="")
      ##
      save(sim_df, file = paste0("Simulation_Results/", DGP, "_N=", N, "_Delta=", delta,".RData"))
    }# delta-loop
  }# N-loop
}# DGP-loop

