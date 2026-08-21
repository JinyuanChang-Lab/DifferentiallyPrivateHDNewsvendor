############################################################################### 
###########################  1. Preparation    ################################ 
############################################################################### 

#-------------------------------------------------------------------------------
#                         1.1 Required packages                         
#------------------------------------------------------------------------------- 
library(here)
library(parallel) 
library(quantregForest)
library(qrnn)

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
test_NN <- function(X, Y, f_0, train_indices, test_indices,
                     b, h, c_w = 0.1, lr = 0.1, T_t = 100,  # regural parameters
                     m = 200, sd_theta0 = 0.1, theta0 = NULL ) {
  
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
  theta0 <- clipping_inf(theta0, 1000)
  
  ## dnn_np
  dnn_out_np <- fit_nn_np(X = X_train, Y = Y_train, eta0 = lr, 
                           T_all = T_t, tau = tau, varpi = w_dnn, theta0 = theta0)
  f_hat_np <- predict_nn(dnn_out_np, X_test)
  ER_nn_np <- nv_loss(f_hat_np, Y_test, b, h) - nv_loss(f0_test, Y_test, b, h)  
  
  ## SAA
  f_hat_saa <- as.numeric(quantile(Y_train, probs = tau, type = 1))
  ER_saa  <- nv_loss(f_hat_saa, Y_test, b, h) - nv_loss(f0_test, Y_test, b, h)  
  
  ## quantregForest
  fit_qrf <- quantregForest(x = X_train, y = Y_train)
  f_hat_qrf <- predict(fit_qrf, X_test, what = tau)
  ER_qrf  <- nv_loss(f_hat_qrf, Y_test, b, h) - nv_loss(f0_test, Y_test, b, h)  
  
  # qrnn
  Xtrain_qrnn <- X_train[, -1, drop = FALSE]
  Xtest_qrnn  <- X_test[, -1, drop = FALSE]
  fit_qrnn <- qrnn.fit(
    x = Xtrain_qrnn,
    y = as.matrix(Y_train),
    n.hidden = 5,
    tau = tau,
    iter.max = 500,
    n.trials = 3,
    trace = FALSE
  )
  f_hat_qrnn <- qrnn.predict(x = Xtest_qrnn, parms = fit_qrnn)
  ER_qrnn  <- nv_loss(f_hat_qrnn, Y_test, b, h) - nv_loss(f0_test, Y_test, b, h)  
  
  return(c(ER_nn_np,ER_saa,ER_qrf,ER_qrnn))
}



################################################################################# 
############################  3. Parallel computing   ###########################  
################################################################################# 

### 1.parameter settings           
n_vec <-   seq(500,2000,by=500)  
p <- 5 
b <- 2                                                          
h <- 1  
tau <- b/(b+h)    
repetitions <- 500       
model_vec <- c( 1,2 )


output_dir <- here::here("results", "figure_data")
dir.create(output_dir,recursive = TRUE,showWarnings = FALSE)
Ncores <- 167
result_list <- vector( "list", length(model_vec) * length(n_vec))
result_index <- 1
for (model in model_vec) {
  for (n in n_vec) {
    re_vec <- mclapply(1:repetitions, function(r){
        cat(
          "Model =", model,
          "| n =", n,
          "| repetition =", r,
          "\n"
        )
        
        set.seed(2025 * r)
        
        n_all  <- n + 30000
        data <- gen_data( n_all, p, tau, model_type = model)
        X   <- data$X
        Y   <- data$Y
        f_0 <- data$f_0
        
        test_indices  <- sample.int(n_all, size = 30000, replace = FALSE)
        train_indices <- setdiff(seq_len(n_all), test_indices)
        
        out  <- test_NN( X, Y, f_0, train_indices, test_indices, b, h)
        
        return(out)
      },
      mc.cores = Ncores
    )
    
    mat <- do.call(rbind, re_vec)
    
    colnames(mat) <- c(
      "ER_nn_np",
      "ER_saa",
      "ER_qrf",
      "ER_qrnn"
    )
    
    result_list[[result_index]] <- data.frame(
      model = model,
      n = n,
      repetition = seq_len(repetitions),
      mat,
      check.names = FALSE
    )
    
    result_index <- result_index + 1
  }
}

results_all <- do.call( rbind, result_list)
output_file <- file.path( output_dir, "Result_nn_np.csv")
write.csv(results_all,file = output_file,row.names = FALSE)

