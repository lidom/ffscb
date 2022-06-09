## Load packages 
library("tidyverse")
library("parallel") 
library("ffscb")
## -----------------------
## R packages of Fabian Telschow
## Need to be installed and loaded 
## Here commented out since they lead to errors at windows-OS
#library("SIRF")
#library("RFT")
#library("SampleFields")
## -----------------------

##
detectCores()
##
nworkers <- 1
##
## Setup ####################################################
p            <- 101
grid         <- make_grid(p, rangevals=c(0,1))
alpha.level  <- 0.05
##
n_reps_H0    <- 50000
n_reps_H1    <- 10000
##
DGP_seq      <- c("DGP1_shift","DGP1_scale","DGP1_local",
                  "DGP2_shift","DGP2_scale","DGP2_local", 
                  "DGP3_shift","DGP3_scale","DGP3_local")[1:3]
##
delta_Nsmall  <- c(0, seq(from = 0.05, to = 0.45, len = 5))
delta_Nlarge  <- c(0, seq(from = 0.02, to = 0.1,  len = 5))
##
N_seq         <- c(15, 100)
##
type          <- c("Bs", "BEc", "FFSCB.z", "FFSCB.t")
## #########################################################


## #########################################################
## Check computation speed:
library("microbenchmark")
library("fdatest")

set.seed(123)

N <- 100; delta <- 0; 
## Mean and Cov of DGP3
mu      <- meanf_shift(grid, 0)
cov.m   <- make_cov_m(cov.f = ffscb::covf.nonst.matern, grid=grid, cov.f.params=c(2, 1/4, 1/4)) 

run_bechmark <- TRUE
if(run_bechmark){
  benchmark_result <- microbenchmark(
    "FF1.z"={
      ## Generate data
      dat     <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
      ## Estimate mean, covariance, and tau
      hat_mu      <- rowMeans(dat)
      hat.cov     <- crossprod(t(dat - hat_mu)) / (N-1)
      hat.cov.mu  <- hat.cov / N
      hat.tau     <- tau_fun(dat)
      
      try(confidence_band(x=hat_mu, cov=hat.cov.mu, tau=hat.tau, df=N-1, 
                          type="FFSCB.z", conf.level=(1-alpha.level), n_int=1))
    },
    "FF2.z"={
      ## Generate data
      dat     <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
      ## Estimate mean, covariance, and tau
      hat_mu      <- rowMeans(dat)
      hat.cov     <- crossprod(t(dat - hat_mu)) / (N-1)
      hat.cov.mu  <- hat.cov / N
      hat.tau     <- tau_fun(dat)
      
      try(confidence_band(x=hat_mu, cov=hat.cov.mu, tau=hat.tau, df=N-1, 
                          type="FFSCB.z", conf.level=(1-alpha.level), n_int=2))
    },
    "FF4.z"={
      ## Generate data
      dat     <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
      ## Estimate mean, covariance, and tau
      hat_mu      <- rowMeans(dat)
      hat.cov     <- crossprod(t(dat - hat_mu)) / (N-1)
      hat.cov.mu  <- hat.cov / N
      hat.tau     <- tau_fun(dat)
      
      try(confidence_band(x=hat_mu, cov=hat.cov.mu, tau=hat.tau, df=N-1, 
                          type="FFSCB.z", conf.level=(1-alpha.level), n_int=4))
    },
    "FF1.t"={
      ## Generate data
      dat     <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
      ## Estimate mean, covariance, and tau
      hat_mu      <- rowMeans(dat)
      hat.cov     <- crossprod(t(dat - hat_mu)) / (N-1)
      hat.cov.mu  <- hat.cov / N
      hat.tau     <- tau_fun(dat) 
      
      try(confidence_band(x=hat_mu, cov=hat.cov.mu, tau=hat.tau, df=N-1, 
                          type="FFSCB.t", conf.level=(1-alpha.level), n_int=1))
    },
    "FF2.t"={
      ## Generate data
      dat     <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
      ## Estimate mean, covariance, and tau
      hat_mu      <- rowMeans(dat)
      hat.cov     <- crossprod(t(dat - hat_mu)) / (N-1)
      hat.cov.mu  <- hat.cov / N
      hat.tau     <- tau_fun(dat)
      
      try(confidence_band(x=hat_mu, cov=hat.cov.mu, tau=hat.tau, df=N-1, 
                          type="FFSCB.t", conf.level=(1-alpha.level), n_int=2))
    },
    "FF4.t"={
      ## Generate data
      dat     <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
      ## Estimate mean, covariance, and tau
      hat_mu      <- rowMeans(dat)
      hat.cov     <- crossprod(t(dat - hat_mu)) / (N-1)
      hat.cov.mu  <- hat.cov / N
      hat.tau     <- tau_fun(dat)
      
      try(confidence_band(x=hat_mu, cov=hat.cov.mu, tau=hat.tau, df=N-1, 
                          type="FFSCB.t", conf.level=(1-alpha.level), n_int=4))
    },
    "BEc"={
      ## Generate data
      dat     <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
      ## Estimate mean, covariance, and tau
      hat_mu      <- rowMeans(dat)
      hat.cov     <- crossprod(t(dat - hat_mu)) / (N-1)
      hat.cov.mu  <- hat.cov / N
      hat.tau     <- tau_fun(dat)
      try(confidence_band(x=hat_mu, cov=hat.cov.mu, tau=hat.tau, df=N-1, 
                          type="BEc", conf.level=(1-alpha.level), n_int=n_int))
    },
    "Bs"={
      ## Generate data
      dat     <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
      ## Estimate mean, covariance, and tau
      hat_mu      <- rowMeans(dat)
      hat.cov     <- crossprod(t(dat - hat_mu)) / (N-1)
      hat.cov.mu  <- hat.cov / N
      hat.tau     <- tau_fun(dat)
      try(confidence_band(x=hat_mu, cov=hat.cov.mu, tau=hat.tau, df=N-1, 
                          type="Bs", conf.level=(1-alpha.level), n_int=n_int))
    },
    "GKF.t"={
      ## Generate data
      dat       <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
      Y_obj     <- SampleFields::RandomField(field=dat, locations=grid)
      try(SIRF::scb_moments(Y = Y_obj, 
                            method = list(name = "GKF", field = "t"),
                            level=(1-alpha.level),
                            se.est    = "estimate"))
    },
    "MultBs"={
      # Generate data
      dat       <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
      Y_obj     <- SampleFields::RandomField(field=dat, locations=grid)
      try(SIRF::scb_moments(Y = Y_obj, 
                            method = list(name = "MultBoot", field = "t"),
                            level=(1-alpha.level),
                            se.est    = "estimate"))
    },
    "IWT"={
      dat       <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
      try((ITP1fourier(t(dat),maxfrequency=10,B=1000)))#fdatest::IWT1(data = t(dat), mu = rep(0,len=length(grid)))
    }#,
  )
}

boxplot(benchmark_result)

summary(benchmark_result, unit="s")


##
for(DGP in DGP_seq) {
  ##
  set.seed(1)
  ##
  for(N in N_seq) {
    ##
    if(any(DGP==c("DGP1_shift","DGP2_shift","DGP3_shift"))){
      ## For H0 (i.e., delta == 0 <=> mu0==mu ) only DGP*i*_shift is needed since
      ## DGP*i*_scale and DGP*i*_local are then equivalent for all i=1,2,3.
      ##
      ## Take the correct delta_seq corresponding to N
      if( N==min(N_seq) ) delta_seq <- delta_Nsmall     else delta_seq <- delta_Nlarge
    }else{
      ## Take the correct delta_seq corresponding to N
      if( N==min(N_seq) ) delta_seq <- delta_Nsmall[-1] else delta_seq <- delta_Nlarge[-1] 
    }
    ##
    for(delta in delta_seq) {# DGP <- "DGP3_shift"; N <- 100; delta <- 0; n_int <- 4; type  <- c("FFSCB.z", "FFSCB.t")
      ## 
      if(grepl("shift", DGP)) { mu0 <- meanf_shift(grid, 0);  mu <- meanf_shift(grid, delta) }
      if(grepl("scale", DGP)) { mu0 <- meanf_scale(grid, 0);  mu <- meanf_scale(grid, delta) }
      if(grepl("local", DGP)) { mu0 <- meanf_rect( grid, 0);  mu <- meanf_rect( grid, delta) }
      names(mu)  <- grid
      names(mu0) <- grid
      ##
      if(grepl("DGP1", DGP)) {# stationary: smooth 
        cov.m     <- make_cov_m(cov.f = ffscb::covf.st.matern,    grid=grid, cov.f.params=c(3/2, 1/4)) # cov2cor(cov.m)[1,p]; range(cov2cor(cov.m))
      }
      if(grepl("DGP2", DGP)) {# stationary: rough
        cov.m     <- make_cov_m(cov.f = ffscb::covf.st.matern,    grid=grid, cov.f.params=c(1/2, 1/4)) # cov2cor(cov.m)[1,p]; range(cov2cor(cov.m))
      }
      if(grepl("DGP3", DGP)) {# non-stationary: from smooth to rough
        cov.m     <- make_cov_m(cov.f = ffscb::covf.nonst.matern, grid=grid, cov.f.params=c(2, 1/4, 1/4)) # cov2cor(cov.m)[1,p]; range(cov2cor(cov.m))
      }
      ## check plot:
      ## sim.dat  <-  make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
      ## matplot(grid, sim.dat, type="l", lty=1); lines(grid, mu, lwd=2)
      ## 
      ## Number of Monte-Carlo repetitions
      n_reps     <- ifelse(delta==0, n_reps_H0, n_reps_H1)
      ##
      start_time <- Sys.time()
      ##
      res_mclapply <- mclapply(1:n_reps, function(reps) {
        ##
        ## Generate data
        dat         <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
        ## Estimate mean, covariance, and tau
        hat_mu      <- rowMeans(dat)
        hat.cov     <- crossprod(t(dat - hat_mu)) / (N-1)
        hat.cov.mu  <- hat.cov / N
        hat.tau     <- tau_fun(dat) 
        
        ## #################################################################################
        ## Confidence bands
        ## #################################################################################
        ##
        ## Multiplyer Boostrap and tGKF as proposed by 
        ## Simultaneous Confidence Bands for Functional Data Using the Gaussian Kinematic Formula
        ## by Fabian J.E. Telschow, Armin Schwartzman
        ##
        Y_obj     <- SampleFields::RandomField(field=dat, locations=grid)
        MultBoot  <- SIRF::scb_moments(Y = Y_obj, 
                                       method = list(name = "MultBoot", field = "t"),
                                       level=(1-alpha.level),
                                       se.est    = "estimate")
        GKF       <- SIRF::scb_moments(Y = Y_obj, 
                                       method = list(name = "GKF", field = "t"),
                                       level=(1-alpha.level),
                                       se.est    = "estimate")
        ## Preparing the bands for joining with other bands
        b_MultBoot           <- MultBoot$scb[,c(3,1)]
        b_GKF                <- GKF$scb[    ,c(3,1)]
        colnames(b_MultBoot) <- c("MultBoot.t.u.0.95", "MultBoot.t.l.0.95")
        colnames(b_GKF)      <- c("GKF.t.u.0.95",      "GKF.t.l.0.95")
        ## #################################################################################
        
        ## 4 intervals
        b <- confidence_band(x=hat_mu, cov=hat.cov.mu, tau=hat.tau, df=N-1, 
                             type=type, conf.level=(1-alpha.level), n_int=4)
        ## 2 intervals
        b2 <- confidence_band(x=hat_mu, cov=hat.cov.mu, tau=hat.tau, df=N-1, 
                              type=c("FFSCB.z", "FFSCB.t"), conf.level=(1-alpha.level), n_int=2)
        
        ## plot(b); lines(x=grid, y=mu, col="blue"); lines(x=grid, y=mu0, lty=2, col="blue")
        ##
        ## Join all bands
        b2           <- b2[,-1]
        colnames(b2) <- c("FF2int.z.u.0.95", "FF2int.z.l.0.95", "FF2int.t.u.0.95", "FF2int.t.l.0.95")
        b            <- cbind(b, b2, b_GKF, b_MultBoot)
        ##
        n_bands      <- (ncol(b) - 1)/2
        ##
        upper_Bands  <- b[,  2*(1:n_bands), drop=FALSE]
        lower_Bands  <- b[,1+2*(1:n_bands), drop=FALSE]
        ##
        exceed_loc   <- upper_Bands < matrix(mu0, nrow=p, ncol=n_bands) |
          lower_Bands > matrix(mu0, nrow=p, ncol=n_bands)
        ##
        ## saving exceedances events ('at least one crossing occurred in [0,1]?')
        excd    <- as.numeric(apply(exceed_loc, 2, function(x){any(x==TRUE)}))
        ##
        ## saving exceedances events in subintervals [0, 1/4], [1/4, 2/4], ... [3/4,1]
        excd4_i <- matrix(0, nrow = n_bands, ncol = 4)
        for(i in 1:4){
          lo <- (i-1)/4
          up <- ifelse(i/4 == 1, 1.01, i/4)
          ##
          excd4_i[,i] <- as.numeric(apply(exceed_loc, 2, function(x){any(x[lo <= grid & grid < up]==TRUE)}))
        }
        ##
        ## saving exceedances events in subintervals [0, 1/2], [1/2, 1]
        excd2_1 <- as.numeric(excd4_i[,1] ==1 | excd4_i[,2] ==1)
        excd2_2 <- as.numeric(excd4_i[,3] ==1 | excd4_i[,4] ==1)
        excd2_i <- cbind(excd2_1, excd2_2)
        ##
        ## widths of the bands:
        avg_width <- colMeans((upper_Bands - lower_Bands))
        ## Names of methods in correct order:
        Band_type <- stringr::str_replace(string  = names(avg_width),
                                          pattern = paste0(".u.", (1-alpha.level)),"")
        Band_type[Band_type=="FFSCB.z"] <- "FF4int.z"
        Band_type[Band_type=="FFSCB.t"] <- "FF4int.t"
        ## simulation data
        sim_df <- dplyr::tibble(band    = as_factor(Band_type),# band types
                                excd    = excd,         # was there an exceedance event at all?
                                ##
                                excd4_i1 = excd4_i[,1], # was there an exceedance event in [  0,1/4]?
                                excd4_i2 = excd4_i[,2], # was there an exceedance event in [1/4,2/4]?
                                excd4_i3 = excd4_i[,3], # was there an exceedance event in [2/4,3/4]?
                                excd4_i4 = excd4_i[,4], # was there an exceedance event in [3/4,  1]?
                                ##
                                excd2_i1 = excd2_i[,1], # was there an exceedance event in [0,  1/2]?
                                excd2_i2 = excd2_i[,2], # was there an exceedance event in [1/2,  1]?
                                wdth    = avg_width)    # width of the bands
        ## glimpse(sim_df)
        return(sim_df)
      }, mc.cores = nworkers)
      ##
      end_time <- Sys.time()
      run_time <- end_time - start_time
      ##
      n_bands <- length(res_mclapply[[1]]$band)
      ##
      sim_df <- dplyr::bind_rows(res_mclapply) %>% ## row-binding n_rep tibbles in 'res_mclapply'
        dplyr::mutate( # Adding variables:
          run       = rep(c(1:n_reps), each=n_bands),   # Numbering simulation runs
          n_rep     = rep(n_reps, times=n_bands*n_reps),# Number of Monte-Carlo simulation
          delta     = rep(delta,  times=n_bands*n_reps),# Delta
          N         = rep(N,      times=n_bands*n_reps),# Sample size
          DGP       = rep(DGP,    times=n_bands*n_reps))# Name of DGP
      ##
      ## Feedback
      cat(DGP, "_N=", N, "_tau_", tau, "_Delta=", delta, 
          ", Run-Time=", run_time, " (", attr(run_time, "units"),")\n", sep="")
      ##
      save(sim_df, file = paste0("Simulation_Results_Revision_1/", DGP, "_N=", N, "_tau_", tau, "_Delta=", delta,".RData"))
    }# delta-loop
  }# N-loop
}# DGP-loop


## Under H0 (mu0 == mu <=> delta == 0) are DGP*i*_scale and DGP*i*_local equivalent for all i=1,2,3. 
## Therefore, we used only one MC-Simulation (DGP*i*_shift).
## The following code adds the results for DGP*i*_scale (delta==0) and DGP*i*_local (delta==0) by 
## filling in the result of DGP*i*_shift (delta==0).
for(i in 1:3){
  for(dgp in c(paste0("DGP",i,"_scale"), paste0("DGP",i,"_local"))) {
    for(N in N_seq) {
      ##
      load(file = paste0("Simulation_Results_Revision_1/DGP",i,"_shift", "_N=", N, "_tau_", tau, "_Delta=0.RData"))
      sim_df <- sim_df %>% mutate(DGP = dgp) # Replace name of DGP
      save(sim_df, file =  paste0("Simulation_Results_Revision_1/", dgp, "_N=", N, "_tau_", tau, "_Delta=0.RData"))
      rm(sim_df)
    }
  }
}