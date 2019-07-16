## #######################################################
##
## Summarizing the simulation results
##
## #######################################################

library("ffscb")
library("tidyverse")
## 

## path to simulation results:
my_plot_path       <- "/home/dom/Dropbox/Forschung/PRJ_OPEN/PRJ_Inference4_FDA_using_RFT/Manuscript/"

## Looping-Variables
DGP_seq       <- c("DGP1_shift","DGP1_scale","DGP1_local",
                   "DGP2_shift","DGP2_scale","DGP2_local", 
                   "DGP3_shift","DGP3_scale","DGP3_local")
##
delta_Nsmall  <- c(0, seq(from = 0.05, to = 0.45, len = 5))
delta_Nlarge  <- c(0, seq(from = 0.02, to = 0.1,  len = 5))
alpha.level   <- 0.05
N_seq         <- c(10,100)
p             <- 101
grid          <- make_grid(p, rangevals=c(0,1))
## 
IWT_bool      <- c(FALSE, TRUE)


## Wrangling
for(IWT in IWT_bool){
  if(IWT){ IWT_SimResults_df <- NULL }else{ SimResults_df <- NULL}
  for(DGP in DGP_seq){# IWT <- TRUE
    for(N in N_seq) {
      if ( N==min(N_seq) ) delta_seq <- delta_Nsmall else delta_seq <- delta_Nlarge
      for(delta in delta_seq) {# DGP <- "DGP3_shift"; N <- 10; delta <- 0
        ## Load sim_df
        if(IWT){
          load(file = paste0(my_path, "Simulation_Results/IWT_", DGP, "_N=", N, "_Delta=", delta, ".RData"))
        }else{
          load(file = paste0(my_path, "Simulation_Results/",     DGP, "_N=", N, "_Delta=", delta, ".RData"))
        }
        ## Compute which share of the difference between mu and mu0 was correctly found
        if(grepl("shift", DGP)) { mu0 <- meanf_shift(grid, 0);  mu <- meanf_shift(grid, delta) }
        if(grepl("scale", DGP)) { mu0 <- meanf_scale(grid, 0);  mu <- meanf_scale(grid, delta) }
        if(grepl("local", DGP)) { mu0 <- meanf_rect( grid, 0);  mu <- meanf_rect( grid, delta) }
        ##
        SimResults_tmp <- sim_df %>% 
          dplyr::group_by(band) %>% 
          dplyr::summarise(rfrq_excd    = mean(excd),
                           avg_width    = mean(wdth),
                           n_rep        = unique(sim_df$n_rep),
                           DGP          = unique(sim_df$DGP),
                           delta        = unique(sim_df$delta),
                           N            = unique(sim_df$N),
                           alpha        = alpha.level) 
        ##
        ## Row-Binding all 'SimResults_tmp' data frames:
        if(IWT){
          IWT_SimResults_df <- SimResults_tmp %>% 
            dplyr::select(band, DGP, N, delta, n_rep, alpha, rfrq_excd) %>% 
            dplyr::bind_rows(IWT_SimResults_df, .)
        }else{
          SimResults_df <- SimResults_tmp %>% 
            dplyr::select(band, DGP, N, delta, n_rep, alpha, rfrq_excd) %>% 
            dplyr::bind_rows(SimResults_df, .)
        }
      }
    }
  }
}

## ###############################################
## Joining the aggregated simulation results
## ###############################################

## 'band'-factor to 'band'-character (needed for the rowbinding)
SimResults_df     <- SimResults_df     %>% mutate(band = as.character(band))
IWT_SimResults_df <- IWT_SimResults_df %>% mutate(band = as.character(band))

Size_and_Power_df <- dplyr::bind_rows(SimResults_df, IWT_SimResults_df) %>% 
  ## back to 'band'-factor
  mutate(band = as.factor(band)) %>% 
  dplyr::arrange(DGP, N, delta)

## ###############################################
## Building data frames  
## ###############################################

## n=10
Size_and_Power_n10_df <- Size_and_Power_df %>% 
  filter(band!="naive.t" & N==10) %>% 
  select(-alpha, -n_rep, -N) %>% 
  spread(delta, rfrq_excd) %>% 
  select(DGP, band, `0`:`0.45`) %>%
  arrange(DGP, band) 

Size_and_Power_n10_df %>% print(n=Inf)

Size_n10_df <- Size_and_Power_n10_df %>% 
  filter(grepl("_local", Size_and_Power_n10_df$DGP)) %>% 
  select(DGP, band, `0`) %>% 
  spread(band, `0`) %>% 
  mutate(DGP=c("DGP1","DGP2","DGP3")) %>% 
  select(DGP, FFSCB.t, FFSCB.z, KR.t, KR.z, IWT, BEc, Bs)

## n=100
Size_and_Power_n100_df <- Size_and_Power_df %>% 
  filter(band!="naive.t" & N==100) %>% 
  select(-alpha, -n_rep, -N) %>% 
  spread(delta, rfrq_excd) %>% 
  select(DGP, band, `0`:`0.1`) %>% 
  arrange(DGP, band) 

Size_and_Power_n100_df %>% print(n=Inf)

Size_n100_df <- Size_and_Power_n100_df %>% 
  filter(grepl("_local", Size_and_Power_n100_df$DGP)) %>% 
  select(DGP, band, `0`) %>% 
  spread(band, `0`) %>% 
  mutate(DGP=c("DGP1","DGP2","DGP3")) %>% 
  select(DGP, FFSCB.t, FFSCB.z, KR.t, KR.z, IWT, BEc, Bs)

## Plots

Size_and_Power_n100_df




## #########################################
## Plots of the Data Generating Processes
## #########################################

## Load packages 
library("tidyverse")
library("parallel")
library("ffscb")
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
type         <- c("Bs", "BEc", "KR.z", "KR.t", "FFSCB.z", "FFSCB.t")
alpha.level  <- 0.05
n_int        <- 10
##
n_reps_H0    <- 50000
n_reps_H1    <- 10000
##
DGP_seq      <- c("DGP1_shift","DGP1_scale","DGP1_local",
                  "DGP2_shift","DGP2_scale","DGP2_local", 
                  "DGP3_shift","DGP3_scale","DGP3_local")[7:9]
##
delta_Nsmall  <- c(0, seq(from = 0.05, to = 0.45, len = 5))
delta_Nlarge  <- c(0, seq(from = 0.02, to = 0.1,  len = 5))
##
N_seq         <- c(10,100)
## #########################################################

##
for(DGP in DGP_seq) {
  ##
  set.seed(1111)
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
        cov.m     <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(3/2, 1, 1/4))
        t0        <- grid[1]
      }
      if(grepl("DGP2", DGP)) {# stationary: rough
        cov.m     <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=c(1/2, 1, 1/4))
        t0        <- grid[1]
      }
      if(grepl("DGP3", DGP)) {# non-stationary: from smooth to rough
        cov.m     <- make_cov_m(cov.f = covf.st.matern.warp.power, grid=grid, cov.f.params=c(4/5, 1, 1/4, 4))
        t0        <- grid[p]
      }
      ## check plot:
      # sim.dat  <-  make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
      # matplot(grid, sim.dat, type="l", lty=1); lines(grid, mu, lwd=2)
      ## 
