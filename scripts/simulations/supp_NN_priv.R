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
source(here::here("src", "supp_functions_NN.R")) 

################################################################################# 
############################  2. Execution algorithm  ###########################  
################################################################################# 
 
#-------------------------------------------------------------------------------
#                               2.1 true f0 models                       
#-------------------------------------------------------------------------------

f0 <- function(X, model_type) { 
  
  x2 <- X[,2]; x3 <- X[,3] 
  x4 <- X[,4]; x5 <- X[,5] 
  z1 <- (x2 + x3) / sqrt(2) 
  z2 <- (x4 - x5) / sqrt(2) 
  if(model_type == 1){
    
    g <- function(u) u * pnorm(  u) 
    s <- function(u) pnorm(  u) - 0.5 
    re <- 0.5 + 0.5 * (g(z1) + g(z2)) +   (s(z1) * s(z2))
    
  }else{
    
    h <- function(u) (u^{3} - 12*u)/8
    re <- 0.5 +  0.5*(h(z1) + h(z2)) 
    
  }
  return(re)
}


#-------------------------------------------------------------------------------
#             2.2 DGP models :y = f0(x) + eps, with Q_tau(eps)=0
#-------------------------------------------------------------------------------

gen_data <- function(n, p, tau, model_type) {
  
  X_org <- matrix(rnorm(n * (p-1)), nrow = n, ncol = p-1) 
  X <- cbind(1,X_org)
  f_0 <- f0(X,model_type) 
  sigma <- 1
  Noise <-  rnorm(n,mean = -sigma*qnorm(tau), sd = sigma)   
  Y <- f_0 + Noise
  list(X = X, Y = Y, f_0 = f_0)
  
}

# newsvendor loss
nv_loss <- function(q, d, b, h) {
  mean(h * pmax(q - d, 0) + b * pmax(d - q, 0))
}

### test
test_NN_priv <- function(X, Y, f_0, train_indices, test_indices,
                          b, h, c_w = 0.1, lr = 0.1, T_t = 100,  # regural parameters
                          m = 200, sd_theta0 = 0.1, theta0 = NULL,  # regural parameters
                          c_B=0.1, c_gamma=0.1, epsilon_1 = 0.5, epsilon_2 = 0.7, epsilon_3 = 0.9, delta = 0.0001   # privacy parameters
) {
  
  # train data and test data
  X_train <- X[train_indices, , drop = FALSE]
  X_test  <- X[test_indices,  , drop = FALSE]
  Y_train <- Y[train_indices]
  Y_test  <- Y[test_indices]
  f0_test <- f_0[test_indices] # true decision
  
  # dimension and training sample size
  p_tr <- ncol(X_train)
  n_tr <- nrow(X_train)
  tau <- b/(b+h)
  
  w_dnn <- c_w * ((p_tr + log(n_tr)) / n_tr)^(0.4) # smooth parameter
  if (is.null(theta0)) { theta0 <- rnorm( (m * (p_tr + 1) +1), sd = sd_theta0) } # initial theta value
  B <- c_B*sqrt(log(n_tr))
  theta0 <- clipping_inf(theta0, B)
  
  ## nn_np
  nn_out_np <- fit_nn_np(X = X_train, Y = Y_train, eta0 = lr, 
                           T_all = T_t, tau = tau, varpi = w_dnn, theta0 = theta0)
  f_hat_np <- predict_nn(nn_out_np, X_test)
  ER_nn_np <-  nv_loss(f_hat_np, Y_test, b, h) - nv_loss(f0_test, Y_test, b, h)
  
  ## nn_priv
  gamma <- c_gamma*sqrt(p_tr+log(n_tr))
  
  
  nn_out_priv_1 <- fit_nn_priv(X = X_train, Y = Y_train, eta0 = lr, 
                                 T_all = T_t, tau = tau, varpi = w_dnn, theta0 = theta0,
                                 gamma = gamma, B = B, epsilon = epsilon_1, delta = delta)
  f_hat_priv_1 <- predict_nn(nn_out_priv_1, X_test)
  ER_nn_priv_1 <- nv_loss(f_hat_priv_1, Y_test, b, h) - nv_loss(f0_test, Y_test, b, h) 
  
  nn_out_priv_2 <- fit_nn_priv(X = X_train, Y = Y_train, eta0 = lr, 
                                 T_all = T_t, tau = tau, varpi = w_dnn, theta0 = theta0,
                                 gamma = gamma, B = B, epsilon = epsilon_2, delta = delta)
  f_hat_priv_2 <- predict_nn(nn_out_priv_2, X_test)
  ER_nn_priv_2 <- nv_loss(f_hat_priv_2, Y_test, b, h) - nv_loss(f0_test, Y_test, b, h)
  
  nn_out_priv_3 <- fit_nn_priv(X = X_train, Y = Y_train, eta0 = lr, 
                                 T_all = T_t, tau = tau, varpi = w_dnn, theta0 = theta0,
                                 gamma = gamma, B = B, epsilon = epsilon_3, delta = delta)
  f_hat_priv_3 <- predict_nn(nn_out_priv_3, X_test)
  ER_nn_priv_3 <- nv_loss(f_hat_priv_3, Y_test, b, h) - nv_loss(f0_test, Y_test, b, h)
  
  
  return(c(ER_nn_np,ER_nn_priv_1,ER_nn_priv_2,ER_nn_priv_3))
}


# ----------------------------------------------------------------------------
#                         parallel computing
# ----------------------------------------------------------------------------




### 1.parameter settings           
n_vec <-  seq(10000,100000 ,by=10000)   
p <- 5 
b <- 2                                                          
h <- 1  
tau <- b/(b+h)   
model_vec <- c(1,2)
repetitions <- 500             


output_dir <- here::here( "results", "figure_data")
dir.create( output_dir, recursive = TRUE, showWarnings = FALSE )
Ncores = 167 # default: 1
Re_mean <- NULL 
for(model in model_vec){
  for(n in n_vec){
    re_vec <- mclapply(1:repetitions, function(r){
      cat(
        "Model =", model,
        "| n =", n,
        "| repetition =", r,
        "\n"
      )
      set.seed(2025 * r)
      
      n_all  <- n + 30000
      data <- gen_data(n_all, p, tau, model_type = model)
      X   <- data$X
      Y   <- data$Y
      f_0 <- data$f_0
      
      test_indices  <- sample.int(n_all, size = 30000, replace = FALSE)
      train_indices <- setdiff(seq_len(n_all), test_indices)
      
      ##  NN
      out  <- test_NN_priv( X, Y, f_0, train_indices, test_indices, b, h)
      
      return(out)
    }, mc.cores = Ncores)
    
    
    mat <- do.call(rbind, re_vec)
    mean_vec <- colMeans(mat, na.rm = TRUE)
    names(mean_vec) <- c("NN_np", "NN_priv_0.5", "NN_priv_0.7", "NN_priv_0.9")
    mean_row <- c(model=model,  n = n,    mean_vec) 
    Re_mean <- rbind(Re_mean, mean_row) 
    
  } # n
}# model


output_file <- file.path( output_dir, "Result_nn_priv.csv" )
write.csv(Re_mean,file = output_file,row.names = FALSE)

