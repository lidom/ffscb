## Load packages 
library("tidyverse")
library("parallel")
library("ffscb")
##
## path to simulation results:
my_path <- "/home/dom/Dropbox/Forschung/PRJ_OPEN/PRJ_Inference4_FDA_using_RFT/"
##
detectCores()
nworkers <- 7


Error_Checker <- function(x){
  bool.result <- inherits(x, "try-error")
  bool.result <- unname(bool.result)
  return(bool.result)
}

## Setup ####################################################
p            <- 101
grid         <- make_grid(p, rangevals=c(0,1))
type         <- c("Bs", "BEc", "KR.z", "KR.t", "FFSCB.z", "FFSCB.t")
alpha.level  <- 0.05
n_int        <- 4
##
n_reps_H0    <- 20000
n_reps_H1    <- 10000
##
DGP_seq      <- c("DGP1_shift","DGP1_scale","DGP1_local",
                  "DGP2_shift","DGP2_scale","DGP2_local", 
                  "DGP3_shift","DGP3_scale","DGP3_local")[7:9]
##
delta_Nsmall  <- c(0, seq(from = 0.05, to = 0.45, len = 5))
delta_Nlarge  <- c(0, seq(from = 0.02, to = 0.1,  len = 5))
##
N_seq         <- c(10, 100)
## #########################################################

##
for(DGP in DGP_seq) {
  ##
  set.seed(123)
  ##
  for(N in N_seq) {
    ##
    if(any(DGP==c("DGP1_shift","DGP2_shift","DGP3_shift"))){
      ## For H0 (i.e., delta_seq == 0 <=> mu0==mu ) only DGP*i*_shift is needed since 
      ## DGP*i*_scale and DGP*i*_local equivalent for all i=1,2,3 if delta==0. 
      if( N==min(N_seq) ) delta_seq <- delta_Nsmall     else delta_seq <- delta_Nlarge     # Take the correct delta_seq corresponding to N
    }else{
      if( N==min(N_seq) ) delta_seq <- delta_Nsmall[-1] else delta_seq <- delta_Nlarge[-1] # Take the correct delta_seq corresponding to N
    }
    ##
    for(delta in delta_seq) {# DGP <- "DGP3_shift"; N <- 100; delta <- 0.1
      ## 
      if(grepl("shift", DGP)) { mu0 <- meanf_shift(grid, 0);  mu <- meanf_shift(grid, delta) }
      if(grepl("scale", DGP)) { mu0 <- meanf_scale(grid, 0);  mu <- meanf_scale(grid, delta) }
      if(grepl("local", DGP)) { mu0 <- meanf_rect( grid, 0);  mu <- meanf_rect( grid, delta) }
      names(mu)  <- grid
      names(mu0) <- grid
      ##
      if(grepl("DGP1", DGP)) {# stationary: smooth 
        cov.m     <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(3/2, 1/4))
        t0        <- 0.5#grid[1]
      }
      if(grepl("DGP2", DGP)) {# stationary: rough
        cov.m     <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(1/2, 1/4))
        t0        <- 0.5#grid[1]
      }
      if(grepl("DGP3", DGP)) {# non-stationary: from smooth to rough
        cov.m     <- make_cov_m(cov.f = covf.nonst.matern, grid=grid, cov.f.params=c(4/2, 1/2, 1/4))
        t0        <- 0.5#grid[50]
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
        check <- TRUE
        while(check){
          ## Generate data
          dat         <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
          ## Estimate mean, covariance, and tau
          hat_mu      <- rowMeans(dat)
          hat.cov     <- crossprod(t(dat - hat_mu)) / (N-1)
          hat.cov.mu  <- hat.cov / N
          hat.tau     <- tau_fun(dat) # plot(y=hat.tau,x=seq(0,1,len=p),type="l")
          ##
          ## Confidence bands
          b <- try(confidence_band(x=hat_mu, cov=hat.cov.mu, tau=hat.tau, t0=t0, df=N-1, 
                                   type=type, conf.level=(1-alpha.level), n_int=n_int), 
                   silent = TRUE)
          if(Error_Checker(b)){ check <- TRUE; cat("Error") } else { check <- FALSE }
        }
        # plot(b); lines(x=grid, y=mu, col="blue"); lines(x=grid, y=mu0, lty=2, col="blue")
        ##
        upper_Bands       <- b[,  2*(1:length(type))]
        lower_Bands       <- b[,1+2*(1:length(type))]
        ##
        exceed_loc        <- upper_Bands < matrix(mu0, nrow=p, ncol=length(type)) | lower_Bands > matrix(mu0, nrow=p, ncol=length(type))
        ##
        ## saving exceedances events ('at least one crossing occured?')
        exceedances       <- as.numeric(apply(exceed_loc, 2, function(x){any(x==TRUE)}))
        ## saving exceedances events per interval [0,1/2] and [1/2,1]
        exceedances_int1 <- numeric(length(type))
        exceedances_int2 <- numeric(length(type))
        ##
        for(j in 1:length(type)){
          exceed_tmp <- sapply(X   = grid[exceed_loc[,j]], 
                               FUN = function(x){unique(
                                 ## The following assures that an exceedance at x=0.5 gets counted for both intervals [0,1/2] and [1/2,1]
                                 c(cut(x=x, breaks=seq(0,1,len=3), labels = c(1,2), include.lowest = TRUE),
                                   cut(x=x, breaks=seq(0,1,len=3), labels = c(1,2), include.lowest = TRUE, right = FALSE)))})
          exceed_tmp <- unlist(exceed_tmp)
          exceedances_int1[j] <- as.numeric(any(exceed_tmp == '1'))
          exceedances_int2[j] <- as.numeric(any(exceed_tmp == '2'))
        }      
        ## saving exceedances at t0:
        tmp_t0_up       <- upper_Bands[which(t0==grid),] < mu0[which(t0==grid)]      
        tmp_t0_lo       <- lower_Bands[which(t0==grid),] > mu0[which(t0==grid)]
        exceedances_t0  <- as.numeric(tmp_t0_up | tmp_t0_lo)
        ##
        ## saving crossing (band vs. mu0) locations:
        ## (If mu0 is completely in or out of the band, there is no crossing location.)
        n_crossing_int1 <- numeric(length(type))
        n_crossing_int2 <- numeric(length(type))
        for(j in 1:length(type)){
          crossings_loc <- rep(FALSE, p)
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
          if( length(tmp_cr_loc)>0 ){ crossings_loc[tmp_cr_loc]  <- TRUE }
          ##
          crossing_intervals <- cut(x=grid[crossings_loc], breaks=seq(0,1,len=3), labels = c(1,2), include.lowest = TRUE)
          n_crossing_int1[j] <- length(crossing_intervals[crossing_intervals == '1'])
          n_crossing_int2[j] <- length(crossing_intervals[crossing_intervals == '2'])
        }
        ## widths of the bands:
        intgr_widths_sqr  <- colSums((upper_Bands - lower_Bands)^2)*diff(grid)[1]
        ## Names of methods in correct order:
        Band_type         <- stringr::str_replace(string  = names(intgr_widths_sqr), pattern = paste0(".u.", (1-alpha.level)),"")
        ## simulation data
        sim_df <- dplyr::tibble(band      = as_factor(Band_type),# band types
                                excd      = exceedances,         # was there an exceedance event at all?
                                excd_i1   = exceedances_int1,    # was there an exceedance event in interval 1?
                                excd_i2   = exceedances_int2,    # was there an exceedance event in interval 2?
                                excd_t0   = exceedances_t0,      # was there an exceedance event at t0?
                                n_cros_i1 = n_crossing_int1,     # number of band-crossing in interval 1?
                                n_cros_i2 = n_crossing_int2,     # number of band-crossing in interval 2?
                                wdth      = intgr_widths_sqr)    # width of the bands 
        ## glimpse(sim_df)
        return(sim_df)
      }, mc.cores = nworkers)
      ##
      end_time <- Sys.time()
      run_time <- end_time - start_time
      ##
      sim_df <- dplyr::bind_rows(res_mclapply) %>%           # row-binding the n_reps-many tibbles contained in 'res_mclapply'
        dplyr::mutate(                                       # Adding variables:
          run       = rep(c(1:n_reps), each=length(type)),   # Numbering the single simulation runs
          n_rep     = rep(n_reps, times=length(type)*n_reps),# Number of Monte-Carlo simulation
          delta     = rep(delta,  times=length(type)*n_reps),# Delta   
          N         = rep(N,      times=length(type)*n_reps),# Sample size
          t0        = rep(t0,     times=length(type)*n_reps),# Value of t0
          DGP       = rep(DGP,    times=length(type)*n_reps))# Name of DGP
      ##
      ## Feedback
      cat(DGP, ", N=", N, ", Delta=", delta, ", Run-Time=", run_time, " (", attr(run_time, "units"),")\n", sep="")
      ##
      save(sim_df, file = paste0(my_path, "Simulation_Results/", DGP, "_N=", N, "_alpha=", alpha.level, "_t0=", t0, "_Delta=", delta,".RData"))
    }# delta-loop
  }# N-loop
}# DGP-loop


## Under H0 (mu0 == mu <=> delta == 0) are DGP*i*_scale and DGP*i*_local equivalent for all i=1,2,3. 
## Therefore, we used only one MC-Simulation (DGP*i*_shift).
## The following code addes the results for DGP*i*_scale (delta==0) and DGP*i*_local (delta==0) by 
## copying the result of DGP*i*_shift (delta==0).
for(i in 1:3){
  for(dgp in c(paste0("DGP",i,"_scale"), paste0("DGP",i,"_local"))) {
    for(N in N_seq) {
      ##
      load(file = paste0(my_path, "Simulation_Results/DGP",i,"_shift", "_N=", N, "_Delta=0.RData"))
      sim_df <- sim_df %>% mutate(DGP = dgp) # Replace name of DGP
      save(sim_df, file = paste0(my_path, "Simulation_Results/", dgp, "_N=", N, "_Delta=0.RData"))
      rm(sim_df)
    }
  }
}


