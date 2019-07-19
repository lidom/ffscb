## #######################################################
##
## Summarizing the simulation results
##
## #######################################################

library("ffscb")
library("tidyverse")
## 

## path to simulation results:
my_path       <- "/home/dom/Dropbox/Forschung/PRJ_OPEN/PRJ_Inference4_FDA_using_RFT/"

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
##
t0            <- 0


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
          load(file = paste0(my_path, "Simulation_Results/",     DGP, "_N=", N, "_alpha=", alpha.level, "_t0=", t0, "_Delta=", delta, ".RData"))
        }
        ##
        SimResults_tmp <- sim_df %>% 
          dplyr::group_by(band) %>% 
          dplyr::summarise(rfrq_excd    = mean(excd),
                           rfrq_excd_i1 = mean(excd_i1),
                           rfrq_excd_i2 = mean(excd_i2),
                           rfrq_excd_i1v2 = mean(excd_i1)/mean(excd_i2),
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




