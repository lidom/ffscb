## Load packages 
library("tidyverse")
library("parallel") 
library("ffscb")
##
detectCores()
nworkers <- 7

## Setup ####################################################
p             <- 101
grid          <- make_grid(p, rangevals=c(0,1))
alpha.level   <- 0.05
n_int         <- c(3,4,8,9)[1]
tol           <- .Machine$double.eps^0.5
##
n_reps_H0     <- 50000
n_reps_H1     <- 10000
##
DGP_seq       <- c("DGP1_shift","DGP1_scale","DGP1_local",
                   "DGP2_shift","DGP2_scale","DGP2_local", 
                   "DGP3_shift","DGP3_scale","DGP3_local")
##
type          <- c("Bs", "BEc", "KR.z", "KR.t", "FFSCB.z", "FFSCB.t")
##
delta_Nsmall  <- c(0, seq(from = 0.05, to = 0.45, len = 5))
delta_Nlarge  <- c(0, seq(from = 0.02, to = 0.1,  len = 5))
##
N_seq         <- c(15, 100)
##
t0_seq        <- c(0, 0.25, 0.50, 0.75, 1)
## #########################################################

##
for(t0 in t0_seq){
  for(DGP in DGP_seq) {
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
      for(delta in delta_seq) {# 
        ## 
        if(grepl("shift", DGP)) { mu0 <- meanf_shift(grid, 0);  mu <- meanf_shift(grid, delta) }
        if(grepl("scale", DGP)) { mu0 <- meanf_scale(grid, 0);  mu <- meanf_scale(grid, delta) }
        if(grepl("local", DGP)) { mu0 <- meanf_rect( grid, 0);  mu <- meanf_rect( grid, delta) }
        names(mu)  <- grid
        names(mu0) <- grid
        ##
        if(grepl("DGP1", DGP)) {# stationary: smooth 
          cov.m     <- make_cov_m(cov.f = covf.st.matern,    grid=grid, cov.f.params=c(3/2, 1/4))
        }
        if(grepl("DGP2", DGP)) {# stationary: rough
          cov.m     <- make_cov_m(cov.f = covf.st.matern,    grid=grid, cov.f.params=c(1/2, 1/4))
        }
        if(grepl("DGP3", DGP)) {# non-stationary: from smooth to rough
          cov.m     <- make_cov_m(cov.f = covf.nonst.matern, grid=grid, cov.f.params=c(2, 1/4, 1/4))
        }
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
          hat.cov     <- crossprod(t(dat - hat_mu)) / (N-1)
          hat.cov.mu  <- hat.cov / N
          hat.tau     <- tau_fun(dat)  
          ##
          ## Confidence bands
          b <- confidence_band(x=hat_mu, cov=hat.cov.mu, tau=hat.tau, t0=t0, df=N-1, 
                               type=type, conf.level=(1-alpha.level), n_int=n_int, tol=tol)
          ##
          upper_Bands       <- b[,  2*(1:length(type)), drop=FALSE]
          lower_Bands       <- b[,1+2*(1:length(type)), drop=FALSE]
          ##
          exceed_loc        <- upper_Bands < matrix(mu0, nrow=p, ncol=length(type)) | lower_Bands > matrix(mu0, nrow=p, ncol=length(type))
          ##
          ## saving exceedances events ('at least one crossing occured?')
          exceedances       <- as.numeric(apply(exceed_loc, 2, function(x){any(x==TRUE)}))
          ##
          if(t0 == 0){
            ## saving exceedances events per interval int1=[0,1/3] and int1=[1,2/3]
            exceedances_int1  <- as.numeric(apply(exceed_loc, 2, function(x){any(x[grid < 1/3]==TRUE)}))
            exceedances_int2  <- as.numeric(apply(exceed_loc, 2, function(x){any(x[grid < 2/3]==TRUE)}))
          }
          if(t0 == 0.25){
            ## saving exceedances events per interval int1=[0,1/3] and int1=[1,2/3]
            exceedances_int1  <- as.numeric(apply(exceed_loc, 2, function(x){any(x[grid <= 1/4]==TRUE)}))
            exceedances_int2  <- as.numeric(apply(exceed_loc, 2, function(x){any(x[grid >= 1/4]==TRUE)}))
          }
          if(t0 == 0.5){
            ## saving exceedances events per interval int1=[0,1/3] and int1=[1,2/3]
            exceedances_int1  <- as.numeric(apply(exceed_loc, 2, function(x){any(x[grid <= 1/2]==TRUE)}))
            exceedances_int2  <- as.numeric(apply(exceed_loc, 2, function(x){any(x[grid >= 1/2]==TRUE)}))
          }
          if(t0 == 0.75){
            ## saving exceedances events per interval int1=[0,1/3] and int1=[1,2/3]
            exceedances_int1  <- as.numeric(apply(exceed_loc, 2, function(x){any(x[grid <= 3/4]==TRUE)}))
            exceedances_int2  <- as.numeric(apply(exceed_loc, 2, function(x){any(x[grid >= 3/4]==TRUE)}))
          }
          ##
          if(t0 == 1){
            ## saving exceedances events per interval int1=[2/3,1] and int1=[1/3,1]
            exceedances_int1  <- as.numeric(apply(exceed_loc, 2, function(x){any(x[grid > 2/3]==TRUE)}))
            exceedances_int2  <- as.numeric(apply(exceed_loc, 2, function(x){any(x[grid > 1/3]==TRUE)}))
          }
          ##
          ## saving exceedances at t0:
          tmp_t0_up       <- upper_Bands[which(t0==grid),, drop=FALSE] < mu0[which(t0==grid)]      
          tmp_t0_lo       <- lower_Bands[which(t0==grid),, drop=FALSE] > mu0[which(t0==grid)]
          exceedances_t0  <- as.numeric(tmp_t0_up | tmp_t0_lo)
          ##
          ## widths of the bands:
          avg_width  <- colMeans((upper_Bands - lower_Bands))
          ## Names of methods in correct order:
          Band_type  <- stringr::str_replace(string  = names(avg_width), pattern = paste0(".u.", (1-alpha.level)),"")
          ## simulation data
          sim_df <- dplyr::tibble(band      = as_factor(Band_type),# band types
                                  excd      = exceedances,         # was there an exceedance event at all?
                                  excd_i1   = exceedances_int1,    # was there an exceedance event in int1?
                                  excd_i2   = exceedances_int2,    # was there an exceedance event in int2?
                                  excd_t0   = exceedances_t0,      # was there an exceedance event at t0?
                                  wdth      = avg_width)           # width of the bands 
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
        cat(DGP, "_N=", N, "_n_int=", n_int, "_t0=", t0, "_Delta=", delta, 
            ", Run-Time=", run_time, " (", attr(run_time, "units"),")\n", sep="")
        ##
        save(sim_df, file = paste0("Simulation_Results/", DGP, "_N=", N, "_n_int=", n_int, "_t0=", t0, "_Delta=", delta,".RData"))
      }# delta-loop
    }# N-loop
  }# DGP-loop
}# t0 loop

## Under H0 (mu0 == mu <=> delta == 0) are DGP*i*_scale and DGP*i*_local equivalent for all i=1,2,3. 
## Therefore, we used only one MC-Simulation (DGP*i*_shift).
## The following code addes the results for DGP*i*_scale (delta==0) and DGP*i*_local (delta==0) by 
## copying the result of DGP*i*_shift (delta==0).
for(i in 1:3){
  for(dgp in c(paste0("DGP",i,"_scale"), paste0("DGP",i,"_local"))) {
    for(N in N_seq) {
      for(t0 in t0_seq) {
        ##
        load(file = paste0("Simulation_Results/DGP",i,"_shift", "_N=", N, "_n_int=", n_int, "_t0=", t0, "_Delta=0.RData"))
        sim_df <- sim_df %>% mutate(DGP = dgp) # Replace name of DGP
        save(sim_df, file =  paste0("Simulation_Results/", dgp, "_N=", N, "_n_int=", n_int, "_t0=", t0, "_Delta=0.RData"))
        rm(sim_df)
      }
    }
  }
}