############################################################################### 
###########################  1. Preparation    ################################ 
############################################################################### 

#-------------------------------------------------------------------------------
#                         1.1 Required packages                         
#------------------------------------------------------------------------------- 
library(here)
library(parallel)  
#-------------------------------------------------------------------------------
#                           1.2 Load functions                         
#------------------------------------------------------------------------------- 
here::i_am("scripts/simulations/supp_grid.R") 
source(here::here("src", "main_functions.R")) 
source(here::here("src", "supp_functions_grid.R"))

################################################################################# 
############################  2. Execution algorithm  ###########################  
################################################################################# 

#-------------------------------------------------------------------------------
#               DGP models :y = f0(x) + eps, with Q_tau(eps)=0
#-------------------------------------------------------------------------------

f0 <- function(X) { 
  
  x2 <- X[,2]; x3 <- X[,3] 
  x4 <- X[,4]; x5 <- X[,5]    
  
  re <- 0.5 +   x2 -   x3 +   x4 * x5 + 0.8* x2* x3 * x4 + 0.8 *exp(x5)
  
  return(re)
}

gen_data <- function(n, p, tau) {
  
  X_org <- matrix(rnorm(n * (p-1)), nrow = n, ncol = p-1)  
  X <- cbind(1,X_org)
  f_0 <- f0(X)    
  sigma <- 1
  Noise <-  rnorm(n,mean = -sigma*qnorm(tau), sd = sigma)   
  Y <- f_0 + Noise
  list(X = X_org, Y = Y)
  
}

test <-  function(X, Y, train_indices, test_indices, s_vec, b, h, n_bins = 3, k =30) {
  
  # train data and test data
  X_train <- X[train_indices, , drop = FALSE]
  X_test  <- X[test_indices,  , drop = FALSE]
  Y_train <- Y[train_indices]
  Y_test  <- Y[test_indices] 
  tau <- b/(b+h)
  
  q_grid <- seq(quantile(Y_train, 0.01, na.rm = TRUE), quantile(Y_train, 0.99, na.rm = TRUE), length.out = 100)
  
  
  # --------------------------------------------------------------------
  # ------------------------------ bin ---------------------------------
  # --------------------------------------------------------------------
  
  fit_bin <- fit_binning_policy(X_train, Y_train, n_bins )
  qhat_bin <- predict_binning_policy(fit_bin, X_test, q_grid, b, h)
  cost_bin <- mean(nv_cost_vec(qhat_bin, Y_test, b, h))
  
  
  # --------------------------------------------------------------------
  # ------------------------------ knn ---------------------------------
  # --------------------------------------------------------------------
  
  fit_knn <- fit_knn_policy(X_train, Y_train)
  qhat_knn <- predict_knn_policy(fit_knn, X_test, k , q_grid, b, h)
  cost_knn <- mean(nv_cost_vec(qhat_knn, Y_test, b, h))
  
  
  # --------------------------------------------------------------------
  # ------------------------------ Linear-SQR -------------------------- 
  # --------------------------------------------------------------------
 
  n_total <- nrow(X_train)
  size_len <- n_total / 5
  init_idx <- sample.int(n_total, size = size_len, replace = FALSE)
  train_idx_new <- setdiff(seq_len(n_total), init_idx)
  
  X_init <- X_train[init_idx, , drop = FALSE]
  Y_init <- Y_train[init_idx]
  
  X_tr <- X_train[train_idx_new, , drop = FALSE]
  Y_tr <- Y_train[train_idx_new]
  
  p_tr <- ncol(X_tr) + 1      # Number of predictors including intercept
  n_tr <- nrow(X_tr)          # Training sample size
  n_init <- nrow(X_init)      # Initial sample size  
  ###### methods ###### 
  c_w = 0.5                                    # Constant for smoothing bandwidth
  lr_1 = 2                                     # Learning rate for warm-up (subsample stage)
  T_1 = 200                                    # Iterations for warm-up
  lr_2 = 0.1     # Step size (learning rate) 
  T_2 =  4*ceiling(log(n_tr /(log (p_tr))))   # Iteration count based on sample/feature size
  BIC_vec <- c()                               # Store BIC values over s
  beta_mat <- NULL                             # Store initial betas over s
  s_vec_raw <- s_vec[s_vec <= p_tr]
  for(s in s_vec_raw){ 
    w_pre <- c_w * ((s * log(p_tr) + log(n_init)) / n_init)^{0.4}                                                                   # Smoothing bandwidth for subsample 
    QT_pre <- QT_highdim(X_init,Y_init,lr=lr_1,beta0=NULL,tau=tau,w=w_pre,s=s,T=T_1,kernel="Gaussian",intercept=TRUE)   # Warm-up sparse QR
    Q_tau <- sum(tau * QT_pre$residuals * (QT_pre$residuals >= 0) + (tau - 1) * QT_pre$residuals * (QT_pre$residuals < 0))        # Quantile check loss
    BIC_vec <- c(BIC_vec,log(Q_tau) + s * log(n_init) * log(p_tr) / (2 * n_init))                                                   # BIC-type criterion
    beta_mat <- rbind(beta_mat,QT_pre$beta)                                                                                       # Store candidate beta
  }
  s_select <- s_vec_raw[which.min(BIC_vec)]                                 # Selected sparsity
  w_select  <- c_w * ((s_select * log(p_tr) + log(n_tr)) / n_tr)^{0.4}  # Bandwidth on full data  
  beta0_init <- beta_mat[which.min(BIC_vec),]                           # Initial beta from best subsample fit
  QT_new <- QT_highdim(X_tr,Y_tr, lr = lr_2, beta0 = beta0_init, tau = tau, w = w_select, s = s_select, T =  T_2,  kernel = "Gaussian", intercept = TRUE)
  SQR_beta <- QT_new$beta  
  X_test_new <- cbind(1,X_test)
  SQR_pred <- X_test_new%*% SQR_beta
  cost_SQR <- mean(nv_cost_vec(SQR_pred, Y_test, b, h))  
  
  
  # --------------------------------------------------------------------
  # ------------------------------ Expand-SQR -------------------------- 
  # --------------------------------------------------------------------
  
  X_train_df <- as.data.frame(X_train)
  X_test_df  <- as.data.frame(X_test)
   
  if (is.null(colnames(X_train_df))) {
    colnames(X_train_df) <- paste0("x", seq_len(ncol(X_train_df)))
    colnames(X_test_df)  <- colnames(X_train_df)
  }
  colnames(X_test_df) <- colnames(X_train_df)
  
  rhs <- paste(colnames(X_train_df), collapse = " + ")
  form_12 <- as.formula(paste0("~ (", rhs, ")^2 - 1"))
 
  X_train_12 <- model.matrix(form_12, data = X_train_df)
  X_test_12  <- model.matrix(form_12, data = X_test_df)
 
  inter_cols <- grepl(":", colnames(X_train_12))
   
  mu_inter <- colMeans(X_train_12[, inter_cols, drop = FALSE])
  sd_inter <- apply(X_train_12[, inter_cols, drop = FALSE], 2, sd)
  sd_inter[sd_inter == 0] <- 1
  
  X_train_12_sc <- X_train_12
  X_test_12_sc  <- X_test_12
  
  X_train_12_sc[, inter_cols] <- scale(
    X_train_12[, inter_cols, drop = FALSE],
    center = mu_inter,
    scale = sd_inter
  )
  
  X_test_12_sc[, inter_cols] <- scale(
    X_test_12[, inter_cols, drop = FALSE],
    center = mu_inter,
    scale = sd_inter
  )
  
  
  X_init_12 <- X_train_12_sc[init_idx, , drop = FALSE] 
  X_tr_12 <- X_train_12_sc[train_idx_new, , drop = FALSE] 
  
  p_tr_12 <- ncol(X_tr_12) + 1      # Number of predictors including intercept 
  ###### methods ######  
  T_2_new =  4*ceiling(log(n_tr /(log (p_tr_12))))   # Iteration count based on sample/feature size
  BIC_vec <- c()                               # Store BIC values over s
  beta_mat <- NULL                             # Store initial betas over s
  s_vec_12  <- s_vec[s_vec <= p_tr_12]
  for(s in s_vec_12){ 
    w_pre <- c_w * ((s * log(p_tr_12) + log(n_init)) / n_init)^{0.4}                                                                   # Smoothing bandwidth for subsample 
    QT_pre <- QT_highdim(X_init_12,Y_init,lr=lr_1,beta0=NULL,tau=tau,w=w_pre,s=s,T=T_1,kernel="Gaussian",intercept=TRUE)   # Warm-up sparse QR
    Q_tau <- sum(tau * QT_pre$residuals * (QT_pre$residuals >= 0) + (tau - 1) * QT_pre$residuals * (QT_pre$residuals < 0))        # Quantile check loss
    BIC_vec <- c(BIC_vec,log(Q_tau) + s * log(n_init) * log(p_tr_12) / (2 * n_init))                                                   # BIC-type criterion
    beta_mat <- rbind(beta_mat,QT_pre$beta)                                                                                       # Store candidate beta
  }
  s_select <- s_vec_12[which.min(BIC_vec)]                                 # Selected sparsity
  w_select  <- c_w * ((s_select * log(p_tr_12) + log(n_tr)) / n_tr)^{0.4}  # Bandwidth on full data  
  beta0_init <- beta_mat[which.min(BIC_vec),]                           # Initial beta from best subsample fit
  QT_new_12 <- QT_highdim(X_tr_12,Y_tr, lr = lr_2, beta0 = beta0_init, tau = tau, w = w_select, s = s_select, T =  T_2_new,  kernel = "Gaussian", intercept = TRUE)
  SQR_beta_12 <- QT_new_12$beta  
  X_test_12_new <- cbind(1,X_test_12_sc)
  SQR_pred_12 <- X_test_12_new%*% SQR_beta_12
  cost_SQR_12 <- mean(nv_cost_vec(SQR_pred_12, Y_test, b, h)) 
  
  return(c( cost_SQR, cost_SQR_12,cost_bin, cost_knn))
  
}  

p_vec <- c(seq(5,9,by=1),seq(10,30,by=5))           
n <- 1000 
b <- 2                                                          
h <- 1 
tau <- b/(b+h) 
repetitions <- 500   
s_vec <- seq(2,30,2) 


output_dir <- here::here( "results", "figure_data")
dir.create( output_dir, recursive = TRUE, showWarnings = FALSE )
Ncores = 167 # default: 1
Re  <- NULL   
for(p in p_vec){  
  ##  parallel 
  re_vec <- rep(NA,repetitions) 
  re_vec <- mclapply(1:repetitions, function(r){            # require 'parallel'
    cat('---- Test Time = ', r, '----- \r')
    set.seed(2025*r)
    # True sparsity level
    n_all  <- n + 5000
    data <- gen_data(n_all, p, tau)
    X   <- data$X
    Y   <- data$Y 
    
    test_indices  <- sample.int(n_all, size = 5000, replace = FALSE)
    train_indices <- setdiff(seq_len(n_all), test_indices)
    
    re <- test(X,  Y ,train_indices,test_indices, s_vec=s_vec, b=b,h=h )    
    return(re) 
    
  }, mc.cores = Ncores)    
  
  Re_1_num <- colMeans(t(matrix(unlist(re_vec), nrow = 4)))  # Return four costs: ell_0 ERM, cs-ERM, Oracle, ERM
  
  Re_1 <- data.frame(
    p  = p,
    linear_SQR = Re_1_num[1],
    Expand_SQR  = Re_1_num[2],
    Bin = Re_1_num[3],
    Knn  = Re_1_num[4] 
  )
  Re <- rbind(Re, Re_1) 
   
} # p


output_file <- file.path( output_dir, "Result_grid.csv" )
write.csv(Re,file = output_file,row.names = FALSE)








