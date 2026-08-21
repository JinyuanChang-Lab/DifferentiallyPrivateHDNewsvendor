############################################################################### 
###########################  1. Preparation    ################################ 
############################################################################### 

#-------------------------------------------------------------------------------
#                         1.1 Required packages                         
#------------------------------------------------------------------------------- 
library(here)
library(parallel)
library(lbfgs)
library(quantreg)

#-------------------------------------------------------------------------------
#                           1.2 Load functions                         
#------------------------------------------------------------------------------- 
here::i_am("scripts/simulations/main_sec5_example1.R")
source( here::here("src", "main_functions.R") ) 



################################################################################# 
############################  2. Execution algorithm  ###########################  
################################################################################# 

#-------------------------------------------------------------------------------
#                               2.1 Section 5.1                        
#-------------------------------------------------------------------------------

# Mean quantile loss
MQL <- function(res,tau){ return(mean(tau * res * (res >= 0) + (tau - 1) * res * (res < 0))) }


#--------------------- main function 1 ------------------------
test1 <- function(X, X0,Y,init_indices,train_indices,test_indices, b=b,h=h, s_vec=s_vec ){
  ##### Inputs ##### 
  # X                : Numeric matrix of dimension n*(p-1). Represents the feature matrix.
  # X0               : Numeric matrix of dimension n*(p0-1). Represents the low-dimensional feature matrix.
  # Y                : Numeric vector of length n. Response variable.
  # train_indices    : Numeric vector. Integer indices for the training subset. 
  # test_indices     : Numeric vector. Integer indices for the test subset.   
  # b                : Numeric scalar. Underage cost.
  # h                : Numeric scalar. Holding cost.
  # s_vec            : Numeric vector. Integer vector of sparsity candidates for model selection via BIC.
  
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
  ###### methods ###### 
  c_w = 0.5                                    # Constant for smoothing bandwidth
  lr_1 = 2                                     # Learning rate for warm-up (subsample stage)
  T_1 = 200                                    # Iterations for warm-up
  lr_2 = 0.1     # Step size (learning rate) 
  T_21 =  4*ceiling(log(n_tr /(log (p_tr))))   # Iteration count based on sample/feature size
  
  n_init <-  nrow(X_init)                      # initial size  
  
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
  
  
  # 1)  SQR with HT on full training set
  QT_new <- QT_highdim(X_train,Y_train, lr = lr_2, beta0 = beta0_init, tau = tau, w = w_select, s = s_select, T =  T_21,  kernel = "Gaussian", intercept = TRUE)
  SQR_beta <- QT_new$beta  
  SQR_pred <- X_test_new%*% SQR_beta
  cost_SQR <- (b+h)*MQL(Y_test-SQR_pred,tau)            # Newsvendor-style cost
  
  # Prepare baseline (low-dimensional) design
  X0_train <- X0[train_indices,]
  X0_test <- X0[test_indices,]  
  X0_test_new <- cbind(1, X0_test)
  
  
  # 2) cs-ERM_new: smoothed ERM via L-BFGS 
  cs_ERM_new <- lbfgs( loss_fn, grad_fn, vars = rep(0, p_tr), X = X_train, Y = Y_train, tau = tau, w = w_select, invisible = TRUE  )
  cs_ERM_new_beta <- cs_ERM_new$par
  cs_ERM_new_pred <- X_test_new%*% cs_ERM_new_beta
  cost_cs_ERM_new <- (b+h)*MQL(Y_test-cs_ERM_new_pred,tau)
  
  # 3) Quantile regression on baseline features
  fit_rq  <- rq.fit(x = cbind(1, X0_train), y = Y_train, tau = tau)
  quantreg_beta <- fit_rq$coefficients
  quantreg_pred <- X0_test_new%*% quantreg_beta
  cost_quantreg <- (b+h)*MQL(Y_test-quantreg_pred,tau) 
  
  # 4) Quantile regression on high-dimensional features
  fit_rq_high  <- rq.fit(x = cbind(1, X_train), y = Y_train, tau = tau)
  quantreg_high_beta <- fit_rq_high$coefficients
  quantreg_high_pred <- X_test_new%*% quantreg_high_beta
  cost_quantreg_high <- (b+h)*MQL(Y_test-quantreg_high_pred,tau) 
  
  
  return( c( cost_SQR,cost_cs_ERM_new,cost_quantreg,cost_quantreg_high )  )  # Return four costs: ell_0 ERM, cs-ERM, Oracle, ERM
}


# ----------------------------------------------------------------------------
#                         parallel computing
# ----------------------------------------------------------------------------
### 1.parameter settings
p_vec <- c(10,  100,  200, 300, 400,500,600,700,800,900)            
n <- 1000 
b <- 2                                                          
h <- 1 
tau <- b/(b+h)
s_star <- 5 
repetitions <- 500   

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
    p_X <- p-1                                                                                    # Number of covariates used in X (excluding intercept)
    itcp <- 4                                                                                     # True intercept
    beta_star <- matrix(c(sample(c(2,-2),s_star-1,replace=TRUE),rep(0,times=p-s_star)),ncol=1)    # True coefficient vector (sparse)
    n_1 = (5*n)/4
    n_all = n_1+30000
    X <-  matrix(rnorm(n_all * p_X), nrow = n_all)                                                # Generate features X (i.i.d. standard normal)
    sigma <- sqrt(2)                                                                              # Noise std
    tau <- b/(b+h)                                                                                # Target quantile level from cost ratio
    Noise <-  rnorm(n_all,mean = -sigma*qnorm(tau), sd = sigma)                                   # Shifted noise so that quantile at tau is zero
    Y <- as.numeric(itcp + X %*% beta_star + Noise)                                                         # Generate response Y
    X0 <- X[,1:(s_star-1)]                                                                        # Low-dimensional baseline features
    
    s_vec <- seq(2 ,10, by = 1)                                                                   # Candidate sparsity levels
    
    test_indices <- sample.int(n_all, size = 30000, replace = FALSE) 
    remain <- setdiff(seq_len(n_all), test_indices) 
    init_indices <- sample(remain, size = n/4, replace = FALSE) 
    train_indices <- setdiff(remain, init_indices)
    
    re <- test1(X, X0,Y,init_indices,train_indices,test_indices, b=b,h=h, s_vec=s_vec )    
    return(re) 
    
  }, mc.cores = Ncores) # number of parallel cores  
  
  Re_1_num <- colMeans(t(matrix(unlist(re_vec), nrow = 4)))  # Return four costs: ell_0 ERM, cs-ERM, Oracle, ERM
  
  Re_1 <- data.frame(
    p  = p,
    ell_0_ERM = Re_1_num[1],
    cs_ERM  = Re_1_num[2],
    Oracle = Re_1_num[3],
    ERM  = Re_1_num[4] 
  )
  Re <- rbind(Re, Re_1) 
   
} # p


output_file <- file.path( output_dir, "Example1.csv" )
write.csv(Re,file = output_file,row.names = FALSE)










