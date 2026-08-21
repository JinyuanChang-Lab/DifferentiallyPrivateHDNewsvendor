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
here::i_am("scripts/simulations/supp_LDP.R")
source(here::here("src", "main_functions.R"))
source(here::here("src", "supp_functions_LDP.R"))

 


################################################################################# 
############################  2. Execution algorithm  ###########################  
################################################################################# 

#-------------------------------------------------------------------------------
#                               2.1 Section 5.3                        
#-------------------------------------------------------------------------------

nv_loss <- function(q,d,b,h){ mean(h*pmax(q-d,0)+b*pmax(d-q,0)) } #  Newsvendor loss: average overage + underage cost

#--------------------- main function 3 ------------------------
CDP_LDP_HDP_comp <- function(X,Y,init_indices,train_indices,test_indices,   b=b,h=h,K_vec=K_vec, s_vec=s_vec, beta_true = beta_true,epsilon=epsilon,delta=delta){ 
  ###### data  #####
  tau = b/(b+h)                                # Target quantile determined by cost ratio
  X_init <- X[init_indices,]
  X_train <- X[train_indices,]                 # Training features (high-dimensional)
  X_test <- X[test_indices,]                   # Test features (high-dimensional)
  Y_init <- Y[init_indices]
  Y_train <- Y[train_indices]                  # Training response
  Y_test <- Y[test_indices]                    # Test response
  X_test_new <- cbind(1, X_test)               # Add intercept to test features
  
  p_tr <- ncol(X_train) + 1                    # Number of predictors including intercept
  n_tr <- nrow(X_train)                        # Training sample size
  ###### parameters ###### 
  c_w = 0.5                                    # Constant for smoothing bandwidth
  lr_1 = 2                                     # Learning rate for warm-up (subsample stage)
  T_1 = 200                                    # Iterations for warm-up
  lr_2 = 0.1     # Step size (learning rate) 
  T_2  =  4*ceiling(log(n_tr /(log (p_tr))))   # Iteration count based on sample/feature size
  
  n_init <-  nrow(X_init)                      # initial size  
  
  # ------------------------- initialization ----------------------------------------
  BIC_vec <- c()                               # Store BIC values over s
  beta_mat <- NULL                             # Store initial betas over s
  for(s in s_vec){ 
    w_pre <- c_w * ((s * log(p_tr) + log(n_init)) / n_init)^{0.4}                                                                   # Smoothing bandwidth for subsample 
    QT_pre <- QT_highdim(X_init,Y_init,lr=lr_1,beta0=NULL,tau=tau,w=w_pre,s=s,T=T_1,kernel="Gaussian",intercept=TRUE)   # Warm-up sparse QR
    Q_tau <- sum(tau * QT_pre$residuals * (QT_pre$residuals >= 0) + (tau - 1) * QT_pre$residuals * (QT_pre$residuals < 0))        # Quantile check loss
    BIC_vec <- c(BIC_vec,log(Q_tau) + s * log(n_init) * log(p_tr) / (2 * n_init))                                                   # BIC-type criterion
    beta_mat <- rbind(beta_mat,QT_pre$beta)                                                                                       # Store candidate beta
  }
  s_select <- s_vec[which.min(BIC_vec)]                                 # Selected sparsity
  w_select  <- c_w * ((s_select * log(p_tr) + log(n_tr)) / n_tr)^{0.4}  # Bandwidth on full data  
  beta0_init <- beta_mat[which.min(BIC_vec),]                           # Initial beta from best subsample fit
  
  # ------------------------- CDP ---------------------------------  
  B_QT <- 0.05 * sqrt(log(p_tr) + log(n_tr))
  QT_cdp <- noisyQT_highdim( X_train, Y_train, epsilon = epsilon, delta = delta,
                             lr = lr_2, beta0 = beta0_init, tau = tau, w = w_select,
                             s = s_select, T = T_2, B = B_QT, kernel = "Gaussian", intercept = TRUE )
  beta_hat_cdp <- QT_cdp$beta
  Y_pred_cdp <- as.numeric(X_test_new %*% beta_hat_cdp)
  Y_star_te <- as.numeric(X_test_new %*% beta_true) 
  excess_cdp <- nv_loss(Y_pred_cdp, Y_test, b, h) - nv_loss(Y_star_te, Y_test, b, h)
  
  # ------------------------- HDP ---------------------------------
  excess_hdp <- c()
  for(K in K_vec){
    QT_hdp <- noisyQT_highdim_hdp(X_train,Y_train,epsilon=epsilon, delta= delta, K=K,lr = lr_2, beta0 = beta0_init, tau = tau,
                                  w = w_select, s = s_select, T =  T_2, B = B_QT,  kernel = "Gaussian", intercept = TRUE) 
    beta_hat_hdp <- QT_hdp$beta
    Y_pred_hdp <- as.numeric(X_test_new %*% beta_hat_hdp) 
    excess_hdp <- c(excess_hdp, nv_loss(Y_pred_hdp, Y_test, b, h) - nv_loss(Y_star_te, Y_test, b, h))
  }
  
  # ------------------------- LDP ---------------------------------
  QT_ldp <- noisyQT_highdim_ldp(X_train,Y_train,epsilon=epsilon, delta= delta, lr = lr_2 , beta0 = beta0_init, tau = tau,
                                w = w_select, s = s_select, T =  T_2, B = B_QT,  kernel = "Gaussian", intercept = TRUE)
  beta_hat_ldp <- QT_ldp$beta
  Y_pred_ldp <- as.numeric(X_test_new %*% beta_hat_ldp) 
  excess_ldp <- nv_loss(Y_pred_ldp, Y_test, b, h) - nv_loss(Y_star_te, Y_test, b, h)
  
  return(c(excess_cdp,excess_hdp,excess_ldp))
}


# ----------------------------------------------------------------------------
#                         parallel computing
# ----------------------------------------------------------------------------
### 1.parameter settings
n_vec <- seq(4000,60000,by=4000)        
p <- 100
b <- 2                                                          
h <- 1   
K_vec <- c(2,4)
K_len <- length(K_vec)
epsilon <- 0.5
delta <- 0.0001
repetitions <- 500                             

output_dir <- here::here( "results", "figure_data")
dir.create( output_dir, recursive = TRUE, showWarnings = FALSE )
Ncores = 167 # default: 1
Re  <- NULL  
for(n in n_vec){
  ##  parallel ---------------- 4 mins -----------------
  re_vec <- rep(NA,repetitions) 
  re_vec <- mclapply(1:repetitions, function(r){            # require 'parallel'
    cat('---- Test Time = ', r, '----- \r')
    set.seed(2025*r)
    # True sparsity level
    s_star <- 5
    p_X <- p-1
    itcp <- 4                                                                                    
    beta_star <- matrix(c(sample(c(2,-2),s_star-1,replace=TRUE),rep(0,times=p-s_star)),ncol=1)
    beta_true <- c(itcp,beta_star)
    n_1 = (5*n)/4
    n_all = n_1+30000
    X <-  matrix(rnorm(n_all * p_X), nrow = n_all)                                                      # Generate features X (i.i.d. standard normal)
    sigma <- sqrt(2)                                                                            # Noise std
    tau <- b/(b+h)                                                                              # Target quantile level from cost ratio
    Noise <-  rnorm(n_all,mean = -sigma*qnorm(tau), sd = sigma)                                     # Shifted noise so that quantile at tau is zero 
    Y <- (itcp + X %*% beta_star) +  Noise  
    
    s_vec <- seq(2 ,10, by = 1)
    test_indices <- sample.int(n_all, size = 30000, replace = FALSE) 
    remain <- setdiff(seq_len(n_all), test_indices) 
    init_indices <- sample(remain, size = n/4, replace = FALSE) 
    train_indices <- setdiff(remain, init_indices)
    
    re <- CDP_LDP_HDP_comp(X, Y,init_indices,train_indices,test_indices,  b=b,h=h, K_vec =K_vec,s_vec=s_vec,beta_true =beta_true, epsilon=epsilon,delta=delta )    
    return(re)
    
  }, mc.cores = Ncores)    
  
  Re_1_num <- colMeans(t(matrix(unlist(re_vec), nrow = 2+K_len )))  # Return four costs: ell_0 ERM, cs-ERM, Oracle, ERM 
  Re <- rbind(Re, Re_1_num) 
   
} # epsilon

colnames(Re) <- c("CDP", paste0("HDP_K", K_vec), "LDP")
output_file <- file.path( output_dir, "Result_LDP.csv" )
write.csv(Re,file = output_file,row.names = FALSE)