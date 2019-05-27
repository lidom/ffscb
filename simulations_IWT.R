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
alpha.level  <- 0.05
n_int        <- 8
##
n_reps_H0    <- 5000
n_reps_H1    <- 5000
##
DGP_seq      <- c("DGP1_shift","DGP1_scale","DGP1_local",
                  "DGP2_shift","DGP2_scale","DGP2_local", 
                  "DGP3_shift","DGP3_scale","DGP3_local",
                  "DGP4_shift","DGP4_scale","DGP4_local")
##
delta_Nsmall  <- c(0, seq(from = 0.05, to = 0.45, len = 5))
delta_Nlarge  <- c(0, seq(from = 0.02, to = 0.1,  len = 5))
##
N_seq         <- c(10,100)
## #########################################################

##
for(DGP in DGP_seq) {
  ##
  set.seed(1110)
  ##
  for(N in N_seq) {
    ##
    if(DGP=="DGP1_shift"){
      ## H0 (i.e., delta_seq == 0 <=> mu0==mu) only one time, since equal for all other DGPs. 
      if( N==min(N_seq) ) delta_seq <- delta_Nsmall     else delta_seq <- delta_Nlarge     # Take the correct delta_seq corresponding to N
    }else{
      if( N==min(N_seq) ) delta_seq <- delta_Nsmall[-1] else delta_seq <- delta_Nlarge[-1] # Take the correct delta_seq corresponding to N
    }
    ##
    for(delta in delta_seq) {# DGP <- DGP_seq[1]; N <- N_seq[1]; delta <- max(delta_Nsmall)
      ## 
      if(grepl("shift", DGP)) { mu0 <- meanf_shift(grid, 0);  mu <- meanf_shift(grid, delta) }
      if(grepl("scale", DGP)) { mu0 <- meanf_scale(grid, 0);  mu <- meanf_scale(grid, delta) }
      if(grepl("local", DGP)) { mu0 <- meanf_rect( grid, 0);  mu <- meanf_rect( grid, delta) }
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
      n_reps       <- ifelse(delta==0, n_reps_H0, n_reps_H1)
      ##
      start_time   <- Sys.time()
      ##
      res_mclapply <- mclapply(1:n_reps, function(reps) {
        ## Generate data
        dat         <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
        ##
        ## ==============================================================================================
        ## IWT1 function from the fdatest package
        IWT_messages                       <- capture.output(IWT_res <- fdatest::IWT1(data = t(dat), mu = mu0))
        exceedances_IWT                    <- as.numeric(any(IWT_res$adjusted_pval < alpha.level))
        exceed_loc_IWT                     <- IWT_res$adjusted_pval < alpha.level
        ##
        if(all(exceed_loc_IWT==TRUE )){ crossings_loc_IWT <- rep( 100,p) }# '100': no crossing since mu0 is completely outside of the band
        if(all(exceed_loc_IWT==FALSE)){ crossings_loc_IWT <- rep(-100,p) }#'-100': no crossing since mu0 is completely inside  of the band
        if(any(exceed_loc_IWT==TRUE) & any(exceed_loc_IWT==FALSE)){   
          tmp_cr_loc_IWT                     <- locate_crossings(IWT_res$adjusted_pval, rep(alpha.level, p), type="down")# down-crossing of pval means upcrossing of empirical mean
          crossings_loc_IWT                  <- rep(NA, times=p)
          crossings_loc_IWT[tmp_cr_loc_IWT]  <- grid[tmp_cr_loc_IWT]
        }
        ##
        ## ==============================================================================================
        ## simulation data
        sim_df <- dplyr::tibble(band     = as_factor(rep("IWT", p)),
                                excd     = rep(exceedances_IWT, p),
                                excd_t0  = rep(NA,              p),
                                excd_loc = exceed_loc_IWT,
                                cros_loc = crossings_loc_IWT,
                                wdth     = rep(NA,              p))# glimpse(sim_df)
        ##
        return(sim_df)
      }, mc.cores = nworkers)
      ##
      end_time <- Sys.time()
      run_time <- end_time - start_time 
      ##
      ## Combine all simulation results
      sim_df <- dplyr::bind_rows(res_mclapply) %>% 
        dplyr::mutate(
          intervals = cut(x=cros_loc, breaks=c(-101,seq(0,1,len=5),101), labels = c('compl_in',1,2,3,4,'compl_out')),# Assign crossing-locations to one of four equidistant intervals in [0,1]
          cros_loc  = replace(cros_loc, cros_loc == -100 | cros_loc == 100, NA),# Setting the 'compl_in/out'  indicators to NA
          run       = rep(c(1:n_reps),    each=p),# Numbering the single simulation runs
          n_rep     = rep(n_reps, times=p*n_reps),# Total number of simulation runs  
          delta     = rep(delta,  times=p*n_reps),# delta values   
          N         = rep(N,      times=p*n_reps),# Sample size
          DGP       = rep(DGP,    times=p*n_reps))# Name of DGP
      ##
      ## Feedback
      cat("IWT_",DGP, ", N=", N, ", Delta=", delta, ", Run-Time=", run_time, " (",attr(run_time, "units"),")\n", sep="")
      ##
      save(sim_df, file = paste0("Simulation_Results/IWT_", DGP, "_N=", N, "_Delta=", delta,".RData"))
    }# delta-loop
  }# N-loop
}# DGP-loop


## Under H0 (mu0 == mu) are all DGPs equivalent, therefore, we use only one MC-Simulation.
for(DGP in DGP_seq[-1]) {
  for(N in N_seq) {
    ##
    load(file = paste0("Simulation_Results/IWT_", DGP_seq[1], "_N=", N, "_Delta=0.RData"))
    sim_df %>% mutate(DGP = rep(DGP, times=length(type)*p*n_reps)) # Replace name of DGP
    save(sim_df, file = paste0("Simulation_Results/IWT_", DGP, "_N=", N, "_Delta=0.RData"))
    rm(sim_df)
  }
}

