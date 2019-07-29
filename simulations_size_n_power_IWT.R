# Load packages 
# devtools::install_github("alessiapini/fdatest")
library("tidyverse")
library("parallel") 
library("ffscb")    
library("fdatest")
##
## path to simulation results:
my_path <- "/home/dom/Dropbox/Forschung/PRJ_OPEN/PRJ_Inference4_FDA_using_RFT/"
##
detectCores()
nworkers <- 6


Error_Checker <- function(x){
  bool.result <- inherits(x, "try-error")
  bool.result <- unname(bool.result)
  return(bool.result)
}


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
                  "DGP3_shift","DGP3_scale","DGP3_local")
##
delta_Nsmall  <- c(0, seq(from = 0.05, to = 0.45, len = 5))
delta_Nlarge  <- c(0, seq(from = 0.02, to = 0.1,  len = 5))
##
N_seq         <- c(15,100)
## #########################################################

##
for(DGP in DGP_seq) {
  ##
  set.seed(1110)
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
    for(delta in delta_seq) {# DGP <- DGP_seq[1]; N <- N_seq[1]; delta <- max(delta_Nsmall)
      ## 
      if(grepl("shift", DGP)) { mu0 <- meanf_shift(grid, 0);  mu <- meanf_shift(grid, delta) }
      if(grepl("scale", DGP)) { mu0 <- meanf_scale(grid, 0);  mu <- meanf_scale(grid, delta) }
      if(grepl("local", DGP)) { mu0 <- meanf_rect( grid, 0);  mu <- meanf_rect( grid, delta) }
      names(mu)  <- grid
      names(mu0) <- grid
      ##
      if(grepl("DGP1", DGP)) {# stationary: smooth 
        cov.m     <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(2, 1/4))
        t0        <- 0
      }
      if(grepl("DGP2", DGP)) {# stationary: rough
        cov.m     <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(1/4, 1/4))
        t0        <- 0
      }
      if(grepl("DGP3", DGP)) {# non-stationary: from smooth to rough
        cov.m     <- make_cov_m(cov.f = covf.nonst.matern, grid=grid, cov.f.params=c(2, 1/4, 1/4))
        t0        <- 0
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
        IWT_messages      <- capture.output(IWT_res <- fdatest::IWT1(data = t(dat), mu = mu0))
        exceedances       <- as.numeric(any(IWT_res$adjusted_pval < alpha.level))
        ##
        ## ==============================================================================================
        ## simulation data
        sim_df <- dplyr::tibble(band     = as_factor("IWT"),# band types
                                excd     = exceedances,         # was there an exceedance event at all?
                                wdth     = NA)    # width of the bands 
        ## glimpse(sim_df)
        return(sim_df)
      }, mc.cores = nworkers)
      ##
      end_time <- Sys.time()
      run_time <- end_time - start_time 
      ##
      ## Combine all simulation results
      sim_df <- dplyr::bind_rows(res_mclapply) %>%           # row-binding the n_reps-many tibbles contained in 'res_mclapply'
        dplyr::mutate(                                       # Adding variables:
          run       = 1:n_reps,                 # Numbering the single simulation runs
          n_rep     = rep(n_reps, times=n_reps),# Number of Monte-Carlo simulation
          delta     = rep(delta,  times=n_reps),# Delta   
          N         = rep(N,      times=n_reps),# Sample size
          DGP       = rep(DGP,    times=n_reps))# Name of DGP
      ##
      ## Feedback
      cat("IWT_",DGP, ", N=", N, ", Delta=", delta, ", Run-Time=", run_time, " (",attr(run_time, "units"),")\n", sep="")
      ##
      save(sim_df, file = paste0(my_path, "Simulation_Results/IWT_", DGP, "_N=", N, "_Delta=", delta,".RData"))
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
      load(file = paste0(my_path, "Simulation_Results/IWT_DGP",i,"_shift", "_N=", N, "_Delta=0.RData"))
      sim_df <- sim_df %>% mutate(DGP = dgp) # Replace name of DGP
      save(sim_df, file = paste0(my_path, "Simulation_Results/IWT_", dgp, "_N=", N, "_Delta=0.RData"))
      rm(sim_df)
    }
  }
}
