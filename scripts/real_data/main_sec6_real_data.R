# =============================================================================
# File: main_sec6_real_data.R
#
# Purpose:
#   Reproduce the empirical results reported in the paper.
#
# Required input:
#   data/processed/processed_data.rds
#
# Outputs: 
#   results/figures/RD_1_2.png
#   results/figures/RD_3.png
#   results/figures/RD_4.png
#
# Prerequisite:
#   Run scripts/real_data/prepare_real_data.R before executing this script.
# =============================================================================


############################################################################### 
###########################  1. Preparation    ################################ 
############################################################################### 

#-------------------------------------------------------------------------------
#                         1.1 Required packages                         
#------------------------------------------------------------------------------- 
library(here) 
library(RSpectra) 
library(lbfgs)
library(foreach)
library(doParallel)
library(doRNG)
library(doSNOW)
library(ggplot2) 
library(tidyr) 
library(ggtext)
library(patchwork)
library(dplyr)
library(cowplot)

#-------------------------------------------------------------------------------
#                        1.2 Load functions and data                       
#-------------------------------------------------------------------------------
here::i_am("scripts/real_data/main_sec6_real_data.R")
source(here::here("src", "main_functions.R"))  # Load functions
processed_data <- readRDS( here::here("data", "processed", "processed_data.rds") ) # Load processed data
RD_X <- processed_data$RD_X # feature
RD_Y <- processed_data$RD_Y # response




################################################################################## 
############################### 2. Numerical performance #########################  
################################################################################## 

# 2.1 using bic, vary $b/h$ and compare sparse private, sparse nonprivate, dense, and non-feature-based
# 2.2 vary sample size. stable performance for our proposed method.
# 2.3 compare the selection of $s$ only on private sparse and highlight our use of BIC.
# 2.4 vary privacy parameter ε.

## common setting 
n_test <- nrow(RD_X)
repetitions <- 500 

## Mean quantile loss
MQL <- function(res, tau) { mean(res * (tau - as.numeric(res < 0))) }

## Remove constant/degenerate columns; keep the remaining columns
col_keep_func <- function(X){
  finite_cnt <- colSums(is.finite(X))                      # Count of finite (non-NA, non-Inf) entries per column
  sds        <- apply(X, 2, sd, na.rm=TRUE)                # Column-wise standard deviation (ignoring NAs)
  keep       <- finite_cnt > 1 & is.finite(sds) & sds > 0  # Keep columns with ≥2 finite values and positive SD
  if (!any(keep)) stop("No columns left after filtering.") # Guard: avoid empty result
  cols_keep <- colnames(X)[keep]                           # Names of columns to retain
  return(cols_keep)                                        # Return surviving column names
} 


#-------------------------------------------------------------------------------
#                         2.1 Vary b (Figure 6, (a))                          
#------------------------------------------------------------------------------- 

comp_test_1 <- function(X, Y, train_indices, test_indices, b = b, h = h,
                         s_vec = s_vec, epsilon = epsilon, delta = delta) {
  
  ########################################## data ##########################################
  tau <- b / (b + h)
  
  # original train/test split
  X_train0 <- X[train_indices, , drop = FALSE]
  X_test   <- X[test_indices , , drop = FALSE]
  
  # delete constant columns using ORIGINAL training set (to keep test aligned)
  cols_keep1 <- col_keep_func(X_train0)
  X_train0   <- X_train0[, cols_keep1, drop = FALSE]
  X_test     <- X_test[,   cols_keep1, drop = FALSE]
  Y_train0 <- Y[train_indices]
  Y_test   <- Y[test_indices]
  X_test_new <- cbind(1, X_test)
  n_tr0 <- nrow(X_train0)
  
  ########################################## methods #######################################
  c_w <- 0.5
  lr_1 <- 2
  T_1  <- 200
  lr_2 <- 0.1
  c_gamma <- 0.05
  # -----------------------------
  # Split ORIGINAL training into init vs main
  # -----------------------------
  sub_n   <- 0.2 * n_tr0 # so that: n_init = 0.25 n_train
  subsamp <- round(sub_n, 0)
  idx_init <- 1:subsamp
  idx_main <- (subsamp + 1):n_tr0
  
  # init sample (ONLY for BIC / beta0)
  X_subsample <- X_train0[idx_init, , drop = FALSE]
  cols_keep2  <- col_keep_func(X_subsample)
  X_subsample <- X_subsample[, cols_keep2, drop = FALSE]
  Y_subsample <- Y_train0[idx_init]
  p_sub <- ncol(X_subsample) + 1
  
  # main training sample (USED by ALL methods in scheme A)
  X_train <- X_train0[idx_main, , drop = FALSE]
  Y_train <- Y_train0[idx_main]
  n_tr <- nrow(X_train)
  p_tr <- ncol(X_train) + 1
  
  # iterations based on MAIN training
  T_2 <- 4 * ceiling(log(n_tr / (log(p_tr)))) 
  
  ########################################## BIC for s and beta0 (init only) #######################################
  BIC_vec  <- c()
  beta_mat <- NULL
  
  for (s in s_vec) {
    w_pre <- c_w * ((s * log(p_sub) + log(sub_n)) / sub_n)^(0.4)
    QT_pre <- QT_highdim(X_subsample, Y_subsample,lr = lr_1, beta0 = NULL, tau = tau, w = w_pre,
                          s = s, T = T_1, kernel = "Gaussian", intercept = TRUE)
    Q_tau <- sum(tau * QT_pre$residuals * (QT_pre$residuals >= 0) +
                   (tau - 1) * QT_pre$residuals * (QT_pre$residuals < 0))
    BIC_vec  <- c(BIC_vec, log(Q_tau) + s * log(sub_n) * log(p_sub) / (2 * sub_n))
    beta_mat <- rbind(beta_mat, QT_pre$beta)
  }
  s_select   <- s_vec[which.min(BIC_vec)]
  beta0_init <- beta_mat[which.min(BIC_vec), ]
  w_select <- c_w * ((s_select * log(p_tr) + log(n_tr)) / n_tr)^(0.4)
  
  ########################################## Align beta0 to MAIN feature space #######################################
  feat_tr  <- colnames(X_train)
  feat_sub <- colnames(X_subsample)
  
  beta0_full <- numeric(p_tr)
  beta0_full[1] <- beta0_init[1]  # intercept
  
  if (!is.null(feat_tr) && !is.null(feat_sub)) {
    pos <- match(feat_sub, feat_tr)
    ok  <- !is.na(pos)
    beta0_full[1 + pos[ok]] <- beta0_init[-1][ok]
  }
  
  ########################################## 1) SQR #######################################
  QT_new <- QT_highdim( X_train, Y_train, lr = lr_2, beta0 = beta0_full, tau = tau, w = w_select, 
                        s = s_select, T = T_2, kernel = "Gaussian", intercept = TRUE )
  SQR_beta <- QT_new$beta
  SQR_pred <- X_test_new %*% SQR_beta
  cost_SQR <- (b + h) * MQL(Y_test - SQR_pred, tau)
  
  ########################################## 2) DP-SQR #######################################
  B_QT <- c_gamma * sqrt(log(p_tr) + log(n_tr))
  QT_priv <- noisyQT_highdim(X_train, Y_train,epsilon = epsilon, delta = delta,lr = lr_2, 
                             beta0 = beta0_full, tau = tau,w = w_select, s = s_select, T = T_2, B = B_QT,
                              kernel = "Gaussian", intercept = TRUE)
  
  DPSQR_beta <- QT_priv$beta
  DPSQR_pred <- X_test_new %*% DPSQR_beta
  cost_DPSQR <- (b + h) * MQL(Y_test - DPSQR_pred, tau)
  
  ########################################## 3) cs-ERM  #######################################
  cs_ERM <- lbfgs( loss_fn, grad_fn, vars = rep(0, p_tr), X = X_train, Y = Y_train,
                   tau = tau, w = w_select)
  cs_ERM_beta <- cs_ERM$par
  cs_ERM_pred <- X_test_new %*% cs_ERM_beta
  cost_cs_ERM <- (b + h) * MQL(Y_test - cs_ERM_pred, tau)
  
  ########################################## 4) SAA #######################################
  hat_q <- as.numeric(quantile(Y_train, probs = tau, type = 1))
  cost_SAA <- mean(b * pmax(Y_test - hat_q, 0) + h * pmax(hat_q - Y_test, 0))
  
  return(c(cost_SQR, cost_DPSQR, cost_cs_ERM, cost_SAA))
}

start <- Sys.time()
wk <- max(1, parallel::detectCores() - 3)
cl <- makeCluster(wk, type = "SOCK") 
registerDoSNOW(cl)
pb <- txtProgressBar(min = 0, max = repetitions, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress, preschedule = FALSE) 
set_seed <- 12345
test_1 <- foreach(r = 1:repetitions, .combine = 'cbind',
                  .options.snow = opts,
                  .options.RNG  = set_seed,
                  .packages = c("lbfgs","RSpectra")) %dorng% {    
                    h     <- 1
                    b_vec <- seq(1 ,19, by = 1) 
                    s_vec <- seq(2 ,20, by = 1)
                    epsilon <- 0.5
                    delta <- 0.0001 
                    train_indices <- sample.int(n_test, size = floor(0.7 * n_test))
                    test_indices  <- setdiff(seq_len(n_test), train_indices)
                    
                    sapply(b_vec, function(b)
                      comp_test_1(RD_X, RD_Y, train_indices, test_indices, b = b, h = h, s_vec = s_vec , epsilon=epsilon,delta = delta))
                  }
close(pb)
stopCluster(cl)
end <- Sys.time()
print(end-start)   

b_vec <- seq(1, 19, by = 1) 
arr <- array(test_1, dim = c( 4 , length(b_vec), repetitions))     
test_1_mean  <- apply(arr, c(1, 2), mean, na.rm = TRUE) 
test_1_mean

L0ERM <- log(test_1_mean[1,])
DPL0ERM <- log(test_1_mean[2,])
csERM <- log(test_1_mean[3,]) 
SAA <- log(test_1_mean[4,]) 

plot_df_1 <- data.frame( b = b_vec, DPL0ERM = DPL0ERM, L0ERM = L0ERM, csERM = csERM,  SAA = SAA)
plot_df_long <- plot_df_1 %>% pivot_longer(cols = c("DPL0ERM", "L0ERM","csERM", "SAA"), names_to = "Method",values_to = "Cost") #library(tidyr) 
method_colors <- c(  "DPL0ERM" = "deepskyblue","L0ERM" = "tomato","csERM" = "#4B3F72", "SAA" = "#2ca02c")

figure1 <- ggplot(plot_df_long, aes(x = b, y = Cost, color = Method)) +
  geom_line(size = 1, linetype = "solid") +
  geom_point(size = 2) +
  labs(x = "b (Underage cost)", y = "Logarithmic average out-of-sample cost") +
  scale_color_manual(
    values = method_colors,
    breaks = c("DPL0ERM", "L0ERM","csERM", "SAA"),
    labels = c( "DP \u2113<sub>0</sub>-ERM",  "\u2113<sub>0</sub>-ERM",   "cs-ERM", "SAA" )
     ) +
  scale_x_continuous(breaks = seq(min(plot_df_1$b), max(plot_df_1$b), by = 1)) +    
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(size = 1.2, color = "black"), 
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold",  colour = "black"), 
    legend.title = element_text(face = "bold", size = 13), 
    legend.text  = element_markdown(size = 12),  
    legend.position = "right"
  ) +
  guides(color = guide_legend(title = "Method"))

#figure1 




#-------------------------------------------------------------------------------
#                         2.2 Vary n (Figure 6, (b))                      
#------------------------------------------------------------------------------- 
comp_test_2 <- function(X, Y,  sub_n_test = sub_n_test, b = b, h = h,  s_vec = s_vec,
                          epsilon = epsilon, delta = delta) {
  
  # -----------------------------
  # 0) Subsample a testbed dataset
  # -----------------------------
  n_test <- length(Y)
  sub_indices <- sample.int(n_test, size = sub_n_test)
  X <- as.matrix(X[sub_indices, , drop = FALSE])
  Y <- as.numeric(Y)[sub_indices]
  
  # train/test split within the subsample
  train_indices <- sample.int(sub_n_test, size = floor(0.7 * sub_n_test))
  test_indices  <- setdiff(seq_len(sub_n_test), train_indices)
  
  ########################################## data ##########################################
  tau <- b / (b + h)
  X_train0 <- X[train_indices, , drop = FALSE]
  X_test   <- X[test_indices , , drop = FALSE]
  
  # delete constant columns using ORIGINAL training set
  cols_keep1 <- col_keep_func(X_train0)
  X_train0   <- X_train0[, cols_keep1, drop = FALSE]
  X_test     <- X_test[,   cols_keep1, drop = FALSE]
  Y_train0 <- Y[train_indices]
  Y_test   <- Y[test_indices]
  X_test_new <- cbind(1, X_test)
  n_tr0 <- nrow(X_train0)
  
  ########################################## methods #######################################
  c_w <- 0.5
  lr_1 <- 2
  T_1  <- 200
  lr_2 <- 0.1
  c_gamma <- 0.05
  # -----------------------------
  # Split ORIGINAL training into init vs main 
  # -----------------------------
  sub_n   <- 0.2 * n_tr0
  subsamp <- round(sub_n, 0)
  
  idx_init <- 1:subsamp
  idx_main <- (subsamp + 1):n_tr0
  
  # init sample (ONLY for BIC / beta0)
  X_subsample <- X_train0[idx_init, , drop = FALSE]
  cols_keep2  <- col_keep_func(X_subsample)
  X_subsample <- X_subsample[, cols_keep2, drop = FALSE]
  Y_subsample <- Y_train0[idx_init]
  p_sub <- ncol(X_subsample) + 1
  
  # main training sample  
  X_train <- X_train0[idx_main, , drop = FALSE]
  Y_train <- Y_train0[idx_main]
  n_tr <- nrow(X_train)
  p_tr <- ncol(X_train) + 1
  
  # iterations based on MAIN training
  T_2 <- 4 * ceiling(log(n_tr / (log(p_tr)))) 
  
  ########################################## BIC ##########################################
  BIC_vec  <- c()
  beta_mat <- NULL
  for (s in s_vec) {
    w_pre <- c_w * ((s * log(p_sub) + log(sub_n)) / sub_n)^(0.4)
    QT_pre <- QT_highdim(X_subsample, Y_subsample,lr = lr_1, beta0 = NULL, tau = tau, w = w_pre,
                          s = s, T = T_1, kernel = "Gaussian", intercept = TRUE)
    Q_tau <- sum(tau * QT_pre$residuals * (QT_pre$residuals >= 0) +
                   (tau - 1) * QT_pre$residuals * (QT_pre$residuals < 0))
    BIC_vec  <- c(BIC_vec, log(Q_tau) + s * log(sub_n) * log(p_sub) / (2 * sub_n))
    beta_mat <- rbind(beta_mat, QT_pre$beta)
  }
  s_select   <- s_vec[which.min(BIC_vec)]
  beta0_init <- beta_mat[which.min(BIC_vec), ]
  w_select <- c_w * ((s_select * log(p_tr) + log(n_tr)) / n_tr)^(0.4)
  
  ########################################## Align beta0 to   feature space #######################################
  feat_tr  <- colnames(X_train)
  feat_sub <- colnames(X_subsample)
  
  beta0_full <- numeric(p_tr)
  beta0_full[1] <- beta0_init[1]
  
  if (!is.null(feat_tr) && !is.null(feat_sub)) {
    pos <- match(feat_sub, feat_tr)
    ok  <- !is.na(pos)
    beta0_full[1 + pos[ok]] <- beta0_init[-1][ok]
  }
  
  ########################################## 1) SQR  #######################################
  QT_new <- QT_highdim(X_train, Y_train,lr = lr_2, beta0 = beta0_full, tau = tau,w = w_select,
                        s = s_select, T = T_2, kernel = "Gaussian", intercept = TRUE)
  SQR_beta <- QT_new$beta
  SQR_pred <- X_test_new %*% SQR_beta
  cost_SQR <- (b + h) * MQL(Y_test - SQR_pred, tau)
  
  ########################################## 2) DP-SQR  #######################################
  B_QT <- c_gamma * sqrt(log(p_tr) + log(n_tr))
  QT_priv <- noisyQT_highdim( X_train, Y_train,epsilon = epsilon, delta = delta,lr = lr_2, 
                              beta0 = beta0_full, tau = tau,w = w_select, s = s_select, 
                              T = T_2, B = B_QT, kernel = "Gaussian", intercept = TRUE )
  DPSQR_beta <- QT_priv$beta
  DPSQR_pred <- X_test_new %*% DPSQR_beta
  cost_DPSQR <- (b + h) * MQL(Y_test - DPSQR_pred, tau)
  
  ########################################## 3) cs-ERM  #######################################
  cs_ERM <- lbfgs( loss_fn, grad_fn, vars = rep(0, p_tr), X = X_train, Y = Y_train,
                     tau = tau, w = w_select )
  cs_ERM_beta <- cs_ERM$par
  cs_ERM_pred <- X_test_new %*% cs_ERM_beta
  cost_cs_ERM <- (b + h) * MQL(Y_test - cs_ERM_pred, tau)
  
  ########################################## 4) SAA  #######################################
  hat_q <- as.numeric(quantile(Y_train, probs = tau, type = 1))
  cost_SAA <- mean(b * pmax(Y_test - hat_q, 0) + h * pmax(hat_q - Y_test, 0))
  
  return(c(cost_SQR, cost_DPSQR, cost_cs_ERM, cost_SAA))
}

start <- Sys.time()
wk <- max(1, parallel::detectCores() - 3)
cl <- makeCluster(wk, type = "SOCK") 
registerDoSNOW(cl)
pb <- txtProgressBar(min = 0, max = repetitions, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress, preschedule = FALSE) 
set_seed <- 12345
test_2 <- foreach(r = 1:repetitions, .combine = 'cbind',
                  .options.snow = opts,
                  .options.RNG  = set_seed, 
                  .packages = c("lbfgs","RSpectra")) %dorng% {    
                    h     <- 1
                    b = 10 
                    epsilon <- 0.5
                    delta <- 0.0001
                    s_vec <- seq(2 ,20, by = 1)
                    sample_sizes <- seq(500,800,by=50)
                    sapply(sample_sizes, function(n_test)
                      comp_test_2(RD_X, RD_Y ,sub_n_test=n_test, b = b, h = h, s_vec = s_vec , epsilon=epsilon,delta=delta))
                  }
close(pb)
stopCluster(cl)
end <- Sys.time()
print(end-start)   

sample_sizes <- seq(500,800,by=50) 
arr <- array(test_2, dim = c( 4 , length(sample_sizes), repetitions))     
test_2_mean  <- apply(arr, c(1, 2), mean, na.rm = TRUE) 
test_2_mean

plot_df2 <- data.frame( n = sample_sizes, L0ERM    = log(test_2_mean[1, ]), DPL0ERM  = log(test_2_mean[2, ]), csERM    = log(test_2_mean[3, ]), SAA = log(test_2_mean[4, ]) )
long2 <- plot_df2 %>% pivot_longer(-n, names_to = "Method", values_to = "Cost") %>% mutate(Panel = "(b) Setting 2") 

method_colors <- c( "DPL0ERM" = "deepskyblue", "L0ERM"  = "tomato", "csERM"   = "#4B3F72", "SAA" = "#2ca02c" ) 
figure2 <- ggplot(long2, aes(x = n, y = Cost, color = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = method_colors,
    breaks = c("DPL0ERM", "L0ERM", "csERM", "SAA"),
    labels = c( "DP \u2113<sub>0</sub>-ERM", "\u2113<sub>0</sub>-ERM", "cs-ERM", "SAA" )
    ) +
  scale_x_continuous(
    breaks = seq(min(long2$n, na.rm = TRUE),
                 max(long2$n, na.rm = TRUE), by = 50)
  ) +
  labs(
    x = "n (Sample size)",
    y = "Logarithmic average out-of-sample cost",
    color = "Method",
    title = "b = 10"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_blank(),
    axis.line.x  = element_line(size = 1.2, color = "black"),
    axis.line.y  = element_line(size = 1.2, color = "black"),
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(face = "bold", color = "black"),
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_markdown(size = 12),
    legend.position = "right"
  ) 
# figure2  

### Combining figure1 & figure2

fig1_no_y <- figure1 + labs(y = NULL)
fig2_no_y <- figure2 + labs(y = NULL)
combined <- (fig1_no_y | fig2_no_y) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_markdown(size = 12),       
    plot.background = element_rect(fill = "white", colour = NA)
  )
final_1_2 <-  ggdraw() +
  draw_plot(combined, x = 0.08, y = 0.04, width = 0.90, height = 0.92) + 
  draw_label("(a)", x = 0.28, y = 0.97, size = 14,  hjust = 0.5) +
  draw_label("(b)", x = 0.68, y = 0.97, size = 14,  hjust = 0.5) +
  draw_label("Logarithmic average out-of-sample cost",
             x = 0.03, y = 0.52, angle = 90, vjust = 1, size = 14, fontface = "bold") +
  theme(plot.background = element_rect(fill = "white", colour = NA))

#final_1_2 

## Save the results  
figure_file_1_2 <- here::here( "results", "figures", "RD_1_2.png" )
ggsave(filename = figure_file_1_2, plot = final_1_2, width = 12, height = 6, units = "in", dpi = 300)

#-------------------------------------------------------------------------------
#                         2.3 Vary s (Figure 7)                         
#-------------------------------------------------------------------------------   
comp_test_3 <- function(X, Y, train_indices, test_indices, b = b, h = h, 
                        s_vec = s_vec, epsilon = epsilon, delta = delta) {
  ########################################## data ##########################################
  tau <- b / (b + h)
  X_train0 <- X[train_indices, , drop = FALSE]
  X_test   <- X[test_indices , , drop = FALSE]
  cols_keep1 <- col_keep_func(X_train0)
  X_train0   <- X_train0[, cols_keep1, drop = FALSE]  # delete constant columns
  X_test     <- X_test[,   cols_keep1, drop = FALSE]
  Y_train0 <- Y[train_indices]
  Y_test   <- Y[test_indices]
  X_test_new <- cbind(1, X_test)
  p_tr0 <- ncol(X_train0) + 1
  n_tr0 <- nrow(X_train0)
  
  ########################################## parameters #######################################
  c_w <- 0.5
  lr_1 <- 2
  T_1  <- 200
  lr_2 <- 0.1
  c_gamma <- 0.05
  # -----------------------------
  # Split training into init vs main
  # -----------------------------
  sub_n   <- 0.2 * n_tr0
  subsamp <- round(sub_n, 0)
  
  idx_init <- 1:subsamp
  idx_main <- (subsamp + 1):n_tr0
  
  # init sample (for BIC / warm start)
  X_subsample <- X_train0[idx_init, , drop = FALSE]
  cols_keep2  <- col_keep_func(X_subsample)
  X_subsample <- X_subsample[, cols_keep2, drop = FALSE]
  Y_subsample <- Y_train0[idx_init]
  p_sub <- ncol(X_subsample) + 1
  
  # main training sample (remove init)
  X_train <- X_train0[idx_main, , drop = FALSE]
  Y_train <- Y_train0[idx_main]
  p_tr <- ncol(X_train) + 1
  n_tr <- nrow(X_train)
  
  # DP iterations based on main training
  T_2 <- 4 * ceiling(log(n_tr / (log(p_tr))))
  
  ########################################## BIC over s (init data only) #######################################
  BIC_vec  <- c()
  beta_mat <- NULL
  for (s in s_vec) {
    w_pre <- c_w * ((s * log(p_sub) + log(sub_n)) / sub_n)^(0.4)
    QT_pre <- QT_highdim(X_subsample, Y_subsample,lr = lr_1, beta0 = NULL, tau = tau, w = w_pre,
                          s = s, T = T_1, kernel = "Gaussian", intercept = TRUE )
    Q_tau <- sum(tau * QT_pre$residuals * (QT_pre$residuals >= 0) +
                   (tau - 1) * QT_pre$residuals * (QT_pre$residuals < 0))
    BIC_vec  <- c(BIC_vec, log(Q_tau) + s * log(sub_n) * log(p_sub) / (2 * sub_n))
    beta_mat <- rbind(beta_mat, QT_pre$beta)
  }
  s_select   <- s_vec[which.min(BIC_vec)]
  beta0_init <- beta_mat[which.min(BIC_vec), ]
  w_select <- c_w * ((s_select * log(p_tr) + log(n_tr)) / n_tr)^(0.4)
  
  ########################################## Align beta0 to main feature space #######################################
  feat_tr  <- colnames(X_train)      # main training features
  feat_sub <- colnames(X_subsample)  # init subsample features
  
  beta0_full <- numeric(p_tr)
  beta0_full[1] <- beta0_init[1]
  
  if (!is.null(feat_tr) && !is.null(feat_sub)) {
    pos <- match(feat_sub, feat_tr)
    ok  <- !is.na(pos)
    beta0_full[1 + pos[ok]] <- beta0_init[-1][ok]
  }
  
  ########################################## DP-SQR with BIC-selected s (fit on main training) #######################################
  B_QT <- c_gamma * sqrt(log(p_tr) + log(n_tr))
  QT_cont_priv <- noisyQT_highdim( X_train, Y_train,epsilon = epsilon, delta = delta,lr = lr_2, 
                                   beta0 = beta0_full, tau = tau,w = w_select, s = s_select, 
                                   T = T_2, B = B_QT,kernel = "Gaussian", intercept = TRUE )
  DPSQR_beta_bic <- QT_cont_priv$beta
  DPSQR_pred_bic <- X_test_new %*% DPSQR_beta_bic
  cost_DPSQR_bic <- (b + h) * MQL(Y_test - DPSQR_pred_bic, tau)
  
  ########################################## DP-SQR for each s in s_vec (fit on main training) #######################################
  rrre <- c()
  
  for (s in s_vec) {
    w <- c_w * ((s * log(p_sub) + log(sub_n)) / sub_n)^(0.4)
    QT_init <- QT_highdim(X_subsample, Y_subsample,lr = lr_1, beta0 = NULL, tau = tau, w = w,
                          s = s, T = T_1, kernel = "Gaussian", intercept = TRUE )
    beta0_init_s <- QT_init$beta
    beta0_full_s <- numeric(p_tr)
    beta0_full_s[1] <- beta0_init_s[1]
    
    if (!is.null(feat_tr) && !is.null(feat_sub)) {
      pos <- match(feat_sub, feat_tr)
      ok  <- !is.na(pos)
      beta0_full_s[1 + pos[ok]] <- beta0_init_s[-1][ok]
    }
    
    B_QT <- c_gamma * sqrt(log(p_tr) + log(n_tr))
    QT_cont_priv_s <- noisyQT_highdim(X_train, Y_train,epsilon = epsilon, delta = delta,lr = lr_2,
                                      beta0 = beta0_full_s, tau = tau,w = w, s = s, T = T_2, 
                                      B = B_QT, kernel = "Gaussian", intercept = TRUE )
    DPSQR_beta <- QT_cont_priv_s$beta
    DPSQR_pred <- X_test_new %*% DPSQR_beta
    cost_DPSQR <- (b + h) * MQL(Y_test - DPSQR_pred, tau)
    rrre <- c(rrre, cost_DPSQR)
  }
  
  return(c(cost_DPSQR_bic, rrre))
}

start <- Sys.time()
wk <- max(1, parallel::detectCores() - 3)
cl <- makeCluster(wk, type = "SOCK") 
registerDoSNOW(cl)
pb <- txtProgressBar(min = 0, max = repetitions, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress, preschedule = FALSE) 
set_seed <- 12345
test_3 <- foreach(r = 1:repetitions, .combine = 'cbind',
                  .options.snow = opts,
                  .options.RNG  = set_seed, 
                  .packages = c("lbfgs","RSpectra")) %dorng% {    
                    h     <- 1  
                    s_vec <- seq(2 ,20, by = 2)
                    b_vec <- c(5,10,15)
                    epsilon <- 0.5
                    delta <- 0.0001 
                    train_indices <- sample.int(n_test, size = floor(0.7 * n_test))
                    test_indices  <- setdiff(seq_len(n_test), train_indices)
                    
                    sapply(b_vec, function(b)
                      comp_test_3(RD_X, RD_Y, train_indices, test_indices, b = b, h = h,  s_vec = s_vec , epsilon=epsilon,delta=delta))
                  }
close(pb)
stopCluster(cl)
end <- Sys.time()
print(end-start)   

s_vec  <-  seq(2,20,by=2)
b_vec  <- c(5,10,15) 
arr <- array(test_3, dim = c( length(s_vec)+1 ,length(b_vec ), repetitions))     
test_3_mean  <- apply(arr, c(1, 2), mean, na.rm = TRUE) 
test_3_mean

test_3_b5_bic <- log(test_3_mean[1,1])
test_3_b10_bic <- log(test_3_mean[1,2])
test_3_b15_bic <- log(test_3_mean[1,3])
test_3_b5_s <- log(test_3_mean[-1,1])
test_3_b10_s <- log(test_3_mean[-1,2])
test_3_b15_s <- log(test_3_mean[-1,3])


plot_df_3 <- data.frame( s   = s_vec , b5  = test_3_b5_s, b10 = test_3_b10_s, b15 = test_3_b15_s )
plot_df_long <- plot_df_3 %>% pivot_longer(cols = c("b5", "b10", "b15"), names_to = "Method", values_to = "Cost")
method_colors <- c("b5" = "deepskyblue", "b10" = "tomato", "b15" = "#2ca02c")
baselines <- data.frame( Method = c("b5","b10","b15"), y = c(test_3_b5_bic, test_3_b10_bic, test_3_b15_bic), Series = "BIC baseline" )

figure_3 <- ggplot() +
  geom_line( data = plot_df_long, aes(x = s, y = Cost, color = Method), linewidth = 1 ) +
  geom_point( data = plot_df_long, aes(x = s, y = Cost, color = Method), size = 2 ) +
   
  geom_hline(
    data = baselines,
    aes(yintercept = y, color = Method, linetype = Series),
    linewidth = 0.9
  ) +
  
  labs( x = "s (Sparsity level)", y = "Logarithmic average out-of-sample cost" ) +
  scale_color_manual(
    name   = "Underage cost",
    values = method_colors,
    breaks = c("b5","b10","b15"),
    labels = c("b = 5","b = 10","b = 15")
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = c("BIC baseline" = "dashed"),
    breaks = "BIC baseline",
    labels = "BIC baseline"
  ) +
  scale_x_continuous( breaks = s_vec  ) +
   
  guides( linetype = guide_legend(override.aes = list(color = "black")) ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.line   = element_line(size = 1.2, color = "black"),
    axis.title  = element_text(face = "bold", size = 14),
    axis.text   = element_text(face = "bold", color = "black"),
    legend.text = element_text(size = 12),
    legend.title= element_blank(),
    legend.position = "right",
    legend.key.width = grid::unit(0.8, "cm")   
  )

#figure_3
  
## Save the results  
figure_file_3 <- here::here( "results", "figures", "RD_3.png" )
ggsave(filename = figure_file_3, plot = figure_3, width = 6, height = 5, units = "in", dpi = 300)


#-------------------------------------------------------------------------------
#                         2.4 Vary ε  (Figure 8)                        
#-------------------------------------------------------------------------------  
comp_test_4 <- function(X, Y, train_indices, test_indices, b = b, h = h, 
                        s_vec = s_vec, epsilon_vec = epsilon_vec, delta = delta) {
  ########################################## data ##########################################
  tau <- b / (b + h)
  X_train0 <- X[train_indices, , drop = FALSE]
  X_test   <- X[test_indices , , drop = FALSE]
  cols_keep1 <- col_keep_func(X_train0)
  X_train0   <- X_train0[, cols_keep1, drop = FALSE]   # delete constant columns
  X_test     <- X_test[,   cols_keep1, drop = FALSE]
  Y_train0 <- Y[train_indices]
  Y_test   <- Y[test_indices]
  X_test_new <- cbind(1, X_test)
  n_tr0 <- nrow(X_train0)
  
  ########################################## methods #######################################
  ## SQR & DP-SQR parameters
  c_w <- 0.5
  lr_1 <- 2
  T_1  <- 200
  lr_2 <- 0.1 
  c_gamma <- 0.05
  # -----------------------------
  # Split training into init vs main
  # -----------------------------
  sub_n   <- 0.2 * n_tr0
  subsamp <- round(sub_n, 0)
  
  idx_init <- 1:subsamp
  idx_main <- (subsamp + 1):n_tr0
  
  # init sample (for BIC / beta0)
  X_subsample <- X_train0[idx_init, , drop = FALSE]
  cols_keep2  <- col_keep_func(X_subsample)
  X_subsample <- X_subsample[, cols_keep2, drop = FALSE]
  Y_subsample <- Y_train0[idx_init]
  p_sub <- ncol(X_subsample) + 1
  
  # main training sample (remove init)
  X_train <- X_train0[idx_main, , drop = FALSE]
  Y_train <- Y_train0[idx_main]
  n_tr <- nrow(X_train)
  p_tr <- ncol(X_train) + 1
  T_2 <- 4 * ceiling(log(n_tr / (log(p_tr))))
  
  ########################################## BIC for s and beta0 #######################################
  BIC_vec  <- c()
  beta_mat <- NULL
  for (s in s_vec) {
    w_pre <- c_w * ((s * log(p_sub) + log(sub_n)) / sub_n)^(0.4)
    QT_pre <- QT_highdim(X_subsample, Y_subsample,lr = lr_1, beta0 = NULL, tau = tau, w = w_pre,
                          s = s, T = T_1, kernel = "Gaussian", intercept = TRUE )
    Q_tau <- sum(tau * QT_pre$residuals * (QT_pre$residuals >= 0) +
                   (tau - 1) * QT_pre$residuals * (QT_pre$residuals < 0))
    BIC_vec  <- c(BIC_vec, log(Q_tau) + s * log(sub_n) * log(p_sub) / (2 * sub_n))
    beta_mat <- rbind(beta_mat, QT_pre$beta)
  }
  s_select   <- s_vec[which.min(BIC_vec)]
  beta0_init <- beta_mat[which.min(BIC_vec), ]
  w_select <- c_w * ((s_select * log(p_tr) + log(n_tr)) / n_tr)^(0.4)
  
  ########################################## Align beta0 to main feature space #######################################
  feat_tr  <- colnames(X_train)
  feat_sub <- colnames(X_subsample)
  
  beta0_full <- numeric(p_tr)
  beta0_full[1] <- beta0_init[1]  # intercept
  
  if (!is.null(feat_tr) && !is.null(feat_sub)) {
    pos <- match(feat_sub, feat_tr)
    ok  <- !is.na(pos)
    beta0_full[1 + pos[ok]] <- beta0_init[-1][ok]
  }
  
  ########################################## 1) SQR on main training #######################################
  QT_new <- QT_highdim( X_train, Y_train, lr = lr_2, beta0 = beta0_full, tau = tau, w = w_select,
                         s = s_select, T = T_2, kernel = "Gaussian", intercept = TRUE )
  SQR_beta <- QT_new$beta
  SQR_pred <- X_test_new %*% SQR_beta
  cost_SQR <- (b + h) * MQL(Y_test - SQR_pred, tau)
  
  ########################################## 2) DP-SQR on main training #######################################
  rre <- c()
  B_QT <- c_gamma * sqrt(log(p_tr) + log(n_tr))
  for (epsilon in epsilon_vec) {
    QT_priv <- noisyQT_highdim( X_train, Y_train, epsilon = epsilon, delta = delta, lr = lr_2,
                                beta0 = beta0_full, tau = tau, w = w_select, s = s_select, 
                                T = T_2, B = B_QT, kernel = "Gaussian", intercept = TRUE )
    DPSQR_beta <- QT_priv$beta
    DPSQR_pred <- X_test_new %*% DPSQR_beta
    cost_DPSQR <- (b + h) * MQL(Y_test - DPSQR_pred, tau)
    rre <- c(rre, cost_DPSQR)
  }
  
  return(c(cost_SQR, rre))
}


start <- Sys.time()
wk <- max(1, parallel::detectCores() - 3)
cl <- makeCluster(wk, type = "SOCK") 
registerDoSNOW(cl)
pb <- txtProgressBar(min = 0, max = repetitions, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress, preschedule = FALSE) 
set_seed <- 12345  #123
test_4 <- foreach(r = 1:repetitions, .combine = 'cbind',
                  .options.snow = opts,
                  .options.RNG  = set_seed, 
                  .packages = c("RSpectra")) %dorng% {    
                    h     <- 1
                    b_vec <- c(10) 
                    s_vec <- seq(2 ,20, by = 1)
                    epsilon_vec <- c( 0.1,0.3,0.5  )
                    delta <- 0.0001  
                    train_indices <- sample.int(n_test, size = floor(0.7 * n_test))
                    test_indices  <- setdiff(seq_len(n_test), train_indices)
                    sapply(b_vec, function(b)
                      comp_test_4(RD_X, RD_Y, train_indices, test_indices, b = b, h = h, s_vec = s_vec , epsilon_vec=epsilon_vec,delta=delta))
                  }
close(pb)
stopCluster(cl)
end <- Sys.time()
print(end-start)  ## 30s


## boxplot data
test_4_re <- as.matrix(test_4)
test_4_np_b10 <- test_4_re[1,] 
test_4_priv01_b10 <- test_4_re[2,]
test_4_priv03_b10 <- test_4_re[3,]
test_4_priv05_b10 <- test_4_re[4,]

df_b10 <- data.frame(  b10_Nonprivate = test_4_np_b10, b10_DP_01 = test_4_priv01_b10, b10_DP_03 = test_4_priv03_b10,  b10_DP_05 = test_4_priv05_b10 )
long_b10 <- df_b10 |> pivot_longer(everything(), names_to = "Method", values_to = "Value") 
method_colors <- c(  "b10_Nonprivate"="tomato", "b10_DP_01" ="deepskyblue", "b10_DP_03" ="#2ca02c", "b10_DP_05" ="#4B3F72" )
x_labels_b10 <- c(
  "b10_Nonprivate" = "<b>Nonprivate</b>",
  "b10_DP_01"      = "<b>*epsilon* = 0.1</b>",
  "b10_DP_03"      = "<b>*epsilon* = 0.3</b>",
  "b10_DP_05"      = "<b>*epsilon* = 0.5</b>"
)
 
figure_4 <- ggplot(long_b10, aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_x_discrete(labels = x_labels_b10) +
  labs(x = NULL, y = "Boxplots of out-of-sample cost", title = "b = 10") +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_blank(),
    panel.grid   = element_blank(),
    axis.line.x  = element_line(colour = "black", linewidth = 1.2),
    axis.line.y  = element_line(colour = "black", linewidth = 1.2),
    axis.text.x  = ggtext::element_markdown(size = 12, colour = "black"),
    axis.text.y  = element_text(face = "bold", size = 12, colour = "black"),
    axis.title.y = element_text(face = "bold", size = 14),
    plot.title   = element_blank()
  ) 
 
#figure_4  

## Save the results  
figure_file_4 <- here::here( "results", "figures", "RD_4.png" )
ggsave(filename = figure_file_4, plot = figure_4, width = 8, height = 4, units = "in", dpi = 300)
