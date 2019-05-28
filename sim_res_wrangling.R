## #######################################################
##
## Data Wrangling (single sim_df's to one SimResults_df)
##
## #######################################################

library("ffscb")
library("tidyverse")
## 
# source("simulations.R")
# source("simulations_IWT.R")
##

## Looping-Variables
DGP_seq       <- c("DGP1_shift","DGP1_scale","DGP1_local",
                   "DGP2_shift","DGP2_scale","DGP2_local", 
                   "DGP3_shift","DGP3_scale","DGP3_local",
                   "DGP4_shift","DGP4_scale","DGP4_local")
delta_Nsmall  <- c(0, seq(from = 0.05, to = 0.45, len = 5))
delta_Nlarge  <- c(0, seq(from = 0.02, to = 0.1,  len = 5))
alpha.level   <- 0.05
N_seq         <- c(10,100)
p             <- 101
grid          <- make_grid(p, rangevals=c(0,1))
## 
## Wrangle IWT-simulation results or all other simulation results?
IWT <- c(FALSE, TRUE)[1]


if(IWT){
  IWT_SimResults_df <- NULL
}else{
  SimResults_df     <- NULL
}


## Wrangling
for(DGP in DGP_seq){
  for(N in N_seq) {
    if ( N==min(N_seq) ) delta_seq <- delta_Nsmall else delta_seq <- delta_Nlarge
    for(delta in delta_seq) {# DGP <- "DGP1_shift"; N <- 10; delta <- 0
      
      ## Load sim_df
      if(IWT){
        load(file = paste0("IWT_Simulation_Results/", DGP, "_N=", N, "_Delta=", delta,".RData"))
      }else{
        load(file = paste0("Simulation_Results/",     DGP, "_N=", N, "_Delta=", delta,".RData"))
      }
      
      ## Compute which share of the difference between mu and mu0 was correctly found
      if(grepl("shift", DGP)) { mu0 <- meanf_shift(grid, 0);  mu <- meanf_shift(grid, delta) }
      if(grepl("scale", DGP)) { mu0 <- meanf_scale(grid, 0);  mu <- meanf_scale(grid, delta) }
      if(grepl("local", DGP)) { mu0 <- meanf_rect( grid, 0);  mu <- meanf_rect( grid, delta) }
      ##
      ## Compute the relative frequency of crossings per interval [0,1/4],[1/4,2/4],[2/4,3/4],[3/4,1]
      SimResults_tmp <- sim_df %>% 
        dplyr::group_by(band) %>% 
        dplyr::summarise(rfrq_excd    = mean(excd),
                         rfrq_excd_t0 = mean(excd_t0),
                         rfrq_excd_i1 = mean(excd_i1),
                         rfrq_excd_i2 = mean(excd_i2),
                         rfrq_excd_i3 = mean(excd_i3),
                         rfrq_excd_i4 = mean(excd_i4),
                         rfrq_cros_i1 = mean(cros_i1),
                         rfrq_cros_i2 = mean(cros_i2),
                         rfrq_cros_i3 = mean(cros_i3),
                         rfrq_cros_i4 = mean(cros_i4),
                         avg_width    = mean(wdth),
                         n_rep        = unique(sim_df$n_rep),
                         DGP          = unique(sim_df$DGP),
                         delta        = unique(sim_df$delta),
                         N            = unique(sim_df$N)) %>% 
        # Note: relative frequencies of zero result in: 0*log(0/.25) => 0*-Inf*0 => NaN
        dplyr::mutate(KL = rfrq_excd_i1*log(rfrq_excd_i1/.25) + 
                        rfrq_excd_i2*log(rfrq_excd_i2/.25) + 
                        rfrq_excd_i3*log(rfrq_excd_i3/.25) + 
                        rfrq_excd_i4*log(rfrq_excd_i4/.25),
                      alpha = alpha.level)  
      ##
      ## Row-Binding all 'SimResults_tmp' data frames:
      if(IWT){
        IWT_SimResults_df <- SimResults_tmp %>% 
          dplyr::select(band, DGP, N, delta, n_rep, alpha, rfrq_excd, rfrq_excd_t0, 
                        rfrq_excd_i1, rfrq_excd_i2, rfrq_excd_i3, rfrq_excd_i4, 
                        rfrq_cros_i1, rfrq_cros_i2, rfrq_cros_i3, rfrq_cros_i4,
                        KL) %>% 
          dplyr::bind_rows(SimResults_df, .)
      }else{
        SimResults_df <- SimResults_tmp %>% 
          dplyr::select(band, DGP, N, delta, n_rep, alpha, rfrq_excd, rfrq_excd_t0, 
                        rfrq_excd_i1, rfrq_excd_i2, rfrq_excd_i3, rfrq_excd_i4, 
                        rfrq_cros_i1, rfrq_cros_i2, rfrq_cros_i3, rfrq_cros_i4,
                        KL) %>% 
          dplyr::bind_rows(SimResults_df, .)
      }
    }
  }
}

## save(SimResults_df,     file = "Simulation_Results/Aggregated_SimResults.RData")
## save(IWT_SimResults_df, file = "Simulation_Results/IWT_Aggregated_SimResults.RData")

load(file = "Simulation_Results/Aggregated_SimResults.RData")

SimResults_df %>% print(n=50)

SimResults_df %>% dplyr::filter(DGP=="DGP1_shift") %>% print(n=Inf)


SimResults_df %>% dplyr::filter(band=="FFSCB.t" & DGP=="DGP1_shift") %>% 
  pull(KL) %>% round(.,digits=2)









# ## Selecting common rejection cases (for having a comparable basis):
# slct_comms <- c(apply(exceed_loc[,,1], 1, function(x){any(x==TRUE)})&
#                   apply(exceed_loc[,,2], 1, function(x){any(x==TRUE)})&
#                   apply(exceed_loc[,,3], 1, function(x){any(x==TRUE)})&
#                   apply(exceed_loc[,,4], 1, function(x){any(x==TRUE)})&
#                   apply(exceed_loc[,,5], 1, function(x){any(x==TRUE)})&
#                   apply(exceed_loc_IWT, 1, function(x){any(x==TRUE)}))
# 
# 
# ## Compute the median of the percentages of correct point-wise decisions within the domain [0,1]
# perc_correct_median <- NULL
# target_loc          <- mu!=mu0
# ##
# for( i in 1:(length(Band_type)-1) ) {
#   perc_correct_median <- c(perc_correct_median,
#                            apply(exceed_loc[slct_comms,,i], 1, function(x){
#                              if(any(x==TRUE)) { tmp <- x == target_loc; length(tmp[tmp==TRUE])/p
#                              } else { NA }}) %>% median(.,na.rm = T))
# }
# perc_correct_median <- c(perc_correct_median,
#                          apply(exceed_loc_IWT[slct_comms,], 1, function(x){
#                            if(any(x==TRUE)) { tmp <- x == target_loc; length(tmp[tmp==TRUE])/p
#                            } else { NA }}) %>% median(.,na.rm = T))
# 
# sim_df_tmp <- dplyr::tibble("DGP"              = DGP,
#                             "N"                = N,
#                             "delta"            = delta,
#                             "Band"             = Band_type,
#                             "alpha.level"      = alpha.level,
#                             "exceed_frq"       = c(count_exceed, count_exceed_IWT) /n_reps,
#                             "KL"               = KL_df$KL_rfrq_interv, 
#                             "rfrq_I1"          = rfrq_interv_df %>% dplyr::filter(intervals==1) %>% pull(rfrq),
#                             "rfrq_I2"          = rfrq_interv_df %>% dplyr::filter(intervals==2) %>% pull(rfrq),
#                             "rfrq_I3"          = rfrq_interv_df %>% dplyr::filter(intervals==3) %>% pull(rfrq),
#                             "rfrq_I4"          = rfrq_interv_df %>% dplyr::filter(intervals==4) %>% pull(rfrq),
#                             ##
#                             "perc_correct_median"  = perc_correct_median,
#                             ##
#                             "exceed_t0_frq"    = c(count_exceed_t0,NA)/n_reps, 
#                             "intgr_widths"     = c(intgr_widths,NA),
#                             "intgr_widths_sqr" = c(intgr_widths_sqr,NA))
# 
# sim_df <- dplyr::bind_rows(sim_df, sim_df_tmp)




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



