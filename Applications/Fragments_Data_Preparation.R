idnum_W <- unique(fragm_df$idnum[fragm_df$ethnic == "White"])
idnum_A <- unique(fragm_df$idnum[fragm_df$ethnic == "Asian"])
idnum_H <- unique(fragm_df$idnum[fragm_df$ethnic == "Hispanic"])
idnum_B <- unique(fragm_df$idnum[fragm_df$ethnic == "Black"])
##
N_W     <- length(idnum_W)
N_A     <- length(idnum_A)
N_H     <- length(idnum_H)
N_B     <- length(idnum_B)
##
# grid_W  <- sort(unique(fragm_df$age[fragm_df$ethnic == "White"]))
# grid_A  <- sort(unique(fragm_df$age[fragm_df$ethnic == "Asian"]))
# grid_H  <- sort(unique(fragm_df$age[fragm_df$ethnic == "Hispanic"]))
# grid_B  <- sort(unique(fragm_df$age[fragm_df$ethnic == "Black"]))
# ##
# xx_min  <- max(min(grid_W), min(grid_A), min(grid_H), min(grid_B))
# xx_max  <- min(max(grid_W), max(grid_A), max(grid_H), max(grid_B))
##
range_m <- spnbmd %>% filter(sex=="mal") %>% .$age %>% range(.)
range_f <- spnbmd %>% filter(sex=="fem") %>% .$age %>% range(.)
##
#grid    <- seq(max(range_f[1],range_m[1]), min(range_f[2],range_m[2]), len=26)
grid    <- seq(9.576, 24.448, len=26)
##
m       <- length(grid)
##
Y_W_mat <- matrix(NA, nrow = m, ncol = N_W)
X_W_mat <- matrix(NA, nrow = m, ncol = N_W)
counter <- 1
for(i in idnum_W){# i <- idnum_W[1]
  xx   <- fragm_df$age[   fragm_df$idnum == i]
  yy   <- fragm_df$spnbmd[fragm_df$idnum == i]
  appr <- approx(x=xx, y=yy, xout = grid[grid >= min(xx) & grid <= max(xx)]) 
  ##
  Y_W_mat[grid >= min(xx) & grid <= max(xx), counter] <- appr$y
  X_W_mat[grid >= min(xx) & grid <= max(xx), counter] <- appr$x
  ##
  counter <- counter + 1
}
##
Y_A_mat <- matrix(NA, nrow = m, ncol = N_A)
X_A_mat <- matrix(NA, nrow = m, ncol = N_A)
counter <- 1
for(i in idnum_A){
  xx   <- fragm_df$age[   fragm_df$idnum == i]
  yy   <- fragm_df$spnbmd[fragm_df$idnum == i]
  appr <- approx(x=xx, y=yy, xout = grid[grid >= min(xx) & grid <= max(xx)]) 
  ##
  Y_A_mat[grid >= min(xx) & grid <= max(xx), counter] <- appr$y
  X_A_mat[grid >= min(xx) & grid <= max(xx), counter] <- appr$x
  ##
  counter <- counter + 1
}
##
Y_H_mat <- matrix(NA, nrow = m, ncol = N_H)
X_H_mat <- matrix(NA, nrow = m, ncol = N_H)
counter <- 1
for(i in idnum_H){
  xx   <- fragm_df$age[   fragm_df$idnum == i]
  yy   <- fragm_df$spnbmd[fragm_df$idnum == i]
  appr <- approx(x=xx, y=yy, xout = grid[grid >= min(xx) & grid <= max(xx)]) 
  ##
  Y_H_mat[grid >= min(xx) & grid <= max(xx), counter] <- appr$y
  X_H_mat[grid >= min(xx) & grid <= max(xx), counter] <- appr$x
  ##
  counter <- counter + 1
}
##
Y_B_mat <- matrix(NA, nrow = m, ncol = N_B)
X_B_mat <- matrix(NA, nrow = m, ncol = N_B)
counter <- 1
for(i in idnum_B){
  xx   <- fragm_df$age[   fragm_df$idnum == i]
  yy   <- fragm_df$spnbmd[fragm_df$idnum == i]
  appr <- approx(x=xx, y=yy, xout = grid[grid >= min(xx) & grid <= max(xx)]) 
  ##
  Y_B_mat[grid >= min(xx) & grid <= max(xx), counter] <- appr$y
  X_B_mat[grid >= min(xx) & grid <= max(xx), counter] <- appr$x
  ##
  counter <- counter + 1
}
