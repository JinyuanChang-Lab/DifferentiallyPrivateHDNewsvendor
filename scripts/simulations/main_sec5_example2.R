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
here::i_am("scripts/simulations/main_sec5_example2.R")
source( here::here("src", "main_functions.R") ) 
 


################################################################################# 
############################  2. Execution algorithm  ###########################  
################################################################################# 

#-------------------------------------------------------------------------------
#                               2.1 Section 5.2                        
#-------------------------------------------------------------------------------

#--------------------- main function 2 ------------------------
test2 <- function(X_init, Y_init, X_tr,  Y_tr, X_te, Y_te,  b=b,h=h, s_vec=s_vec, epsilon=epsilon,delta=delta){ 
  ##### Inputs ##### 
  # X_tr                : Numeric matrix of dimension n_tr*(p-1). Represents the training feature matrix.
  # Y_tr                : Numeric vector of length n_tr. Training response variable.
  # X_te                : Numeric matrix of dimension n_te*(p-1). Represents the test feature matrix. 
  # Y_te                : Numeric vector of length n_te. Test response variable. 
  # b                   : Numeric scalar. Underage cost.
  # h                   : Numeric scalar. Holding cost.
  # s_vec               : Numeric vector. Integer vector of sparsity candidates for model selection via BIC.
  # epsilon/delta       : Numeric scalar. Differential privacy budget.
  
  
  ##### data #####
  tau = b/(b+h)                                                  # Target quantile determined by cost ratio
  cstar <- b - h                                                 # Unit margin (b - h)
  X_test_new <- cbind(1, X_te)                                   # Add intercept to test features
  p_tr <- ncol(X_tr) + 1                                         # Number of predictors including intercept
  n_tr <- nrow(X_tr)                                             # Training sample size
  
  ##### methods #####  
  c_w = 0.5                                                      # Constant for smoothing bandwidth                         
  lr_1 = 2                                                       # Learning rate for warm-up (subsample stage)
  T_1 = 200                                                      # Iterations for warm-up
  lr_2 = 0.1                                                     # Step size (learning rate) 
  T_22 =  4*ceiling(log(n_tr   /(log (p_tr))))                   # Iteration count based on sample/feature size 
  
  n_init <-  nrow(X_init)
  BIC_vec <- c()                                                 # Store BIC values over s
  beta_mat <- NULL                                               # Store initial betas over s
  for(s in s_vec){ 
    w_pre <- c_w * ((s * log(p_tr) + log(n_init)) / n_init)^{0.4}                                                                       # Smoothing bandwidth for subsample 
    QT_pre <- QT_highdim(X_init, Y_init,lr=lr_1,beta0=NULL,tau=tau,w=w_pre,s=s,T=T_1,kernel="Gaussian",intercept=TRUE)      # Warm-up sparse QR 
    Q_tau <- sum(tau * QT_pre$residuals * (QT_pre$residuals >= 0) + (tau - 1) * QT_pre$residuals * (QT_pre$residuals < 0))            # Quantile check loss  
    BIC_vec <- c(BIC_vec,log(Q_tau) + s * log(n_init) * log(p_tr) / (2 * n_init))                                                       # BIC-type criterion
    beta_mat <- rbind(beta_mat,QT_pre$beta)                                                                                           # Store candidate beta
  }
  s_select <- s_vec[which.min(BIC_vec)]                          # Selected sparsity
  w_select  <- c_w*((s_select*log(p_tr)+log(n_tr))/n_tr)^{0.4}   # Bandwidth on full data  
  beta0_init <- beta_mat[which.min(BIC_vec),]                    # Initial beta from best subsample fit
  
  c_gamma  <- 0.05                       # Candidate clipping scales 
  B_QT <- c_gamma * sqrt(log(p_tr) + log(n_tr))
  QT_cont_priv <- noisyQT_highdim(X_tr,Y_tr,epsilon=epsilon, delta= delta, lr = lr_2, beta0 = beta0_init, tau = tau,
                                  w = w_select, s = s_select, T =  T_22, B = B_QT,  kernel = "Gaussian", intercept = TRUE)
  DPSQR_beta <- QT_cont_priv$beta
  Y_te_pred <- X_test_new%*% DPSQR_beta                          # Predicted response on test set
  
  sold   <- pmin(Y_te, Y_te_pred)                                # Sold units
  over   <- pmax(Y_te_pred - Y_te, 0)                            # Overstock
  under  <- pmax(Y_te - Y_te_pred, 0)                            # Stockout
  revenue <- mean(cstar * sold)                                  # Average revenue
  cost <- mean(h * over + b * under )                            # Average cost
  pi <- revenue   -  cost                                        # Average profit
  
  return( c(revenue,cost,pi) )                                   # Return revenue, cost, profit
}


# ----------------------------------------------------------------------------
#                         parallel computing
# ----------------------------------------------------------------------------
### 1.parameter settings
epsilon_vec <- seq(0.1,0.9,by=0.1)           
n <- 1000 
p <- 200
b <- 2                                                          
h <- 1  
s_star <- 5 
delta <- 0.0001
repetitions <- 500                             

output_dir <- here::here( "results", "figure_data")
dir.create( output_dir, recursive = TRUE, showWarnings = FALSE )
Ncores = 167 # default: 1
Re  <- NULL  
for(epsilon in epsilon_vec){
  ##  parallel ---------------- 7 mins -----------------
  re_vec <- rep(NA,repetitions) 
  re_vec <- mclapply(1:repetitions, function(r){            # require 'parallel'
    cat('---- Test Time = ', r, '----- \r')
    set.seed(2025*r)
    # True sparsity level
    s_star <- 5                                                                                 # True sparsity level                 
    p_X <- p-1                                                                                  # Number of covariates used in X (excluding intercept)
    itcp <- 4                                                                                   # True intercept
    beta_star <- matrix(c(sample(c(2,-2),s_star-1,replace=TRUE),rep(0,times=p-s_star)),ncol=1)  # True coefficient vector (sparse)
    n_1 = (5*n)/4
    n_all = n_1+30000
    X <-  matrix(rnorm(n_all * p_X), nrow = n_all)                                                      # Generate features X (i.i.d. standard normal)
    sigma <- sqrt(2)                                                                            # Noise std
    tau <- b/(b+h)                                                                              # Target quantile level from cost ratio
    Noise <-  rnorm(n_all,mean = -sigma*qnorm(tau), sd = sigma)                                     # Shifted noise so that quantile at tau is zero 
    Y <- (itcp + X %*% beta_star)*(1+ 3/(1+exp(4*(epsilon**2) ))) +  Noise                      # Generate response Y with epsilon-driven scaling  
    s_vec <- seq(2 ,10, by = 1)                                                                 # Candidate sparsity levels
    
    test_indices <- sample.int(n_all, size = 30000, replace = FALSE) 
    remain <- setdiff(seq_len(n_all), test_indices) 
    init_indices <- sample(remain, size = n/4, replace = FALSE) 
    train_indices <- setdiff(remain, init_indices)
    
    X_init <- X[init_indices,] 
    X_tr <- X[train_indices,]                                                                   # Training features
    X_te <- X[test_indices,]                                                                    # Test features
    Y_init <- Y[init_indices]      
    Y_tr <- Y[train_indices]                                                                    # Training response
    Y_te <- Y[test_indices]                                                                     # Test response
    
    
    re <- test2(X_init, Y_init, X_tr,  Y_tr, X_te, Y_te, b=b,h=h, s_vec=s_vec, epsilon=epsilon,delta=delta)     # Run evaluation    
    return(re) 
    
  }, mc.cores = Ncores) # number of parallel cores   
  
  Re_1_num <- colMeans(t(matrix(unlist(re_vec), nrow = 3)))  # Return four costs: ell_0 ERM, cs-ERM, Oracle, ERM
  
  Re_1 <- data.frame(
    eps  = epsilon,
    revenue = Re_1_num[1],
    cost  = Re_1_num[2],
    pi = Re_1_num[3]  
  )
  Re <- rbind(Re, Re_1) 
   
} # epsilon

output_file <- file.path( output_dir, "Example2.csv" )
write.csv(Re,file = output_file,row.names = FALSE)

