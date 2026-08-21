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
here::i_am("scripts/simulations/main_sec5_example4.R")
source( here::here("src", "main_functions.R") )  


################################################################################# 
############################  2. Execution algorithm  ###########################  
################################################################################# 

#-------------------------------------------------------------------------------
#                               2.1 Section 5.4                        
#-------------------------------------------------------------------------------

#--------------------- main function 4 ------------------------
test4 <- function(
    X_init, X_tr, X_te,
    Y_cont_tr_init, Y_cont_tr, Y_cont_te,
    Y_atq_tr_init,  Y_atq_tr,  Y_atq_te,
    Y_away_tr_init, Y_away_tr, Y_away_te,
    b = b, h = h, s_vec = s_vec,
    beta_true = beta_true,
    epsilon = epsilon, delta = delta
) {
  
  # ----------------------------
  # data / constants
  # ----------------------------
  tau <- b / (b + h)
  X_test_new <- cbind(1, X_te)
  
  p_tr <- ncol(X_tr) + 1
  n_tr <- nrow(X_tr)
  n_init <- nrow(X_init)
  
  # oracle decision on test set: q*(x) = (1,x)^T beta_true
  Y_star_te <- as.numeric(X_test_new %*% beta_true)
  
  # newsvendor loss
  nv_loss <- function(q, d, b, h) {
    mean(h * pmax(q - d, 0) + b * pmax(d - q, 0))
  }
  
  # ----------------------------
  # algorithm hyper-parameters
  # ----------------------------
  c_w   <- 0.5
  lr_1  <- 2
  T_1   <- 200
  lr_2  <- 0.1
  T_main <- 4 * ceiling(log(n_tr / log(p_tr)))
  
  # ----------------------------
  # helper: one scenario
  # ----------------------------
  eval_one_dgp <- function(Y_tr_init, Y_tr, Y_te) {
    
    # ---- 1) BIC selection on init sample ----
    BIC_vec <- numeric(0)
    beta_mat <- NULL
    
    for (s in s_vec) {
      w_pre <- c_w * ((s * log(p_tr) + log(n_init)) / n_init)^0.4
      QT_pre <- QT_highdim( X_init, Y_tr_init,lr = lr_1, beta0 = NULL,tau = tau, w = w_pre,s = s, T = T_1, kernel = "Gaussian", intercept = TRUE )
      res <- QT_pre$residuals
      Q_tau <- sum(tau * res * (res >= 0) + (tau - 1) * res * (res < 0))
      BIC_vec <- c(BIC_vec,log(Q_tau + 1e-12) + s * log(n_init) * log(p_tr) / (2 * n_init))
      beta_mat <- rbind(beta_mat, QT_pre$beta)
    }
    
    idx <- which.min(BIC_vec)
    s_select  <- s_vec[idx]
    beta0_init <- beta_mat[idx, ]
    
    # bandwidth for main training
    w_select <- c_w * ((s_select * log(p_tr) + log(n_tr)) / n_tr)^0.4
    
    # ---- 2) Non-private SQR on training sample ----
    QT_fit <- QT_highdim(X_tr, Y_tr,lr = lr_2, beta0 = beta0_init,tau = tau, w = w_select,s = s_select, T = T_main,
                          kernel = "Gaussian", intercept = TRUE)
    
    beta_hat <- QT_fit$beta
    Y_pred <- as.numeric(X_test_new %*% beta_hat)
    
    # sample excess risk (truncate at 0 to avoid negative due to finite-sample noise)
    excess <- nv_loss(Y_pred, Y_te, b, h) - nv_loss(Y_star_te, Y_te, b, h)
    excess <- pmax(excess, 0)
    
    # ---- 3) DP-SQR on training sample ----
    c_gamma <- 0.05
    B_QT <- c_gamma * sqrt(log(p_tr) + log(n_tr))
    
    QT_priv <- noisyQT_highdim(X_tr, Y_tr,epsilon = epsilon, delta = delta, lr = lr_2, beta0 = beta0_init,tau = tau,
                               w = w_select,s = s_select, T = T_main,B = B_QT,kernel = "Gaussian", intercept = TRUE)
    
    beta_hat_dp <- QT_priv$beta
    Y_pred_dp <- as.numeric(X_test_new %*% beta_hat_dp)
    
    excess_dp <- nv_loss(Y_pred_dp, Y_te, b, h) - nv_loss(Y_star_te, Y_te, b, h)
    excess_dp <- pmax(excess_dp, 0)
    
    c(excess, excess_dp)
  }
  
  # ----------------------------
  # three DGPs
  # ----------------------------
  res_cont <- eval_one_dgp(Y_cont_tr_init, Y_cont_tr, Y_cont_te)
  res_atq  <- eval_one_dgp(Y_atq_tr_init,  Y_atq_tr,  Y_atq_te)
  res_away <- eval_one_dgp(Y_away_tr_init, Y_away_tr, Y_away_te)
  
  c(res_cont, res_atq, res_away)
}


# ----------------------------------------------------------------------------
#                         parallel computing
# ----------------------------------------------------------------------------
### settings
n_vec <- seq(3000, 8000, by = 1000)
p <- 200
b <- 2
h <- 1
pi <- 0.5
s_star <- 5
epsilon <- 0.5
delta <- 1e-4
repetitions <- 500
 

output_dir <- here::here( "results", "figure_data")
dir.create( output_dir, recursive = TRUE, showWarnings = FALSE )
Ncores = 167 # default: 1
Re <- NULL
for (n in n_vec) {
  
  re_list <- mclapply(1:repetitions, function(r) {
    cat("---- Test Time =", r, "----- \r")
    set.seed(2025 * r)
    
    tau <- b/(b+h)
    p_X <- p - 1
    itcp <- 4
    
    beta_star <- matrix( c(sample(c(2, -2), s_star - 1, replace = TRUE),  rep(0, times = p - s_star)),  ncol = 1 )
    beta_true <- c(itcp, beta_star)
    
    n_1   <- (5 * n) / 4
    n_all <- n_1 + 30000
    
    X <- matrix(rnorm(n_all * p_X), nrow = n_all)
    q_star <- as.numeric(itcp + X %*% beta_star)
    
    sigma <- sqrt(2)
    Noise <- rnorm(n_all, mean = -sigma * qnorm(tau), sd = sigma)
    
    # split indices: test / init / train
    test_indices <- sample.int(n_all, size = 30000, replace = FALSE)
    remain <- setdiff(seq_len(n_all), test_indices)
    init_indices <- sample(remain, size = n/4, replace = FALSE)
    train_indices <- setdiff(remain, init_indices)
    
    X_init <- X[init_indices, , drop = FALSE]
    X_tr   <- X[train_indices, , drop = FALSE]
    X_te   <- X[test_indices, , drop = FALSE]
    
    # continuous
    Y_cont <- as.numeric(q_star + Noise)
    Y_cont_tr_init <- Y_cont[init_indices]
    Y_cont_tr      <- Y_cont[train_indices]
    Y_cont_te      <- Y_cont[test_indices]
    
    # atq / away
    u <- runif(n_all)
    s_vec <- seq(2, 10, by = 1)
    
    Y_atq <- Y_cont
    flag_atq <- (u < pi)
    Y_atq[flag_atq] <- q_star[flag_atq]
    
    a <- 1
    Y_away <- Y_cont
    left_flag  <- (u < (pi * tau))
    right_flag <- (u >= (pi * tau)) & (u < pi)
    Y_away[left_flag]  <- q_star[left_flag]  - a
    Y_away[right_flag] <- q_star[right_flag] + a
    
    Y_atq_tr_init  <- Y_atq[init_indices]
    Y_atq_tr       <- Y_atq[train_indices]
    Y_atq_te       <- Y_atq[test_indices]
    
    Y_away_tr_init <- Y_away[init_indices]
    Y_away_tr      <- Y_away[train_indices]
    Y_away_te      <- Y_away[test_indices]
    
    # ---- robust wrapper: never crash the whole parallel batch ----
    out <- tryCatch({
      test4(
        X_init, X_tr, X_te,
        Y_cont_tr_init, Y_cont_tr, Y_cont_te,
        Y_atq_tr_init,  Y_atq_tr,  Y_atq_te,
        Y_away_tr_init, Y_away_tr, Y_away_te,
        b = b, h = h, s_vec = s_vec,
        beta_true = beta_true,
        epsilon = epsilon, delta = delta
      )
    }, error = function(e) {
      message(sprintf("[ERROR] n=%d r=%d : %s", n, r, conditionMessage(e)))
      rep(NA_real_, 6)
    })
    
    # enforce numeric length-6 output
    if (!is.numeric(out) || length(out) != 6) {
      message(sprintf("[ERROR] n=%d r=%d : invalid output (type/length).", n, r))
      out <- rep(NA_real_, 6)
    }
    
    out
  }, mc.cores = Ncores)  # number of parallel cores
  
  # re_list: length = repetitions, each element numeric(6) or NA_real_(6)
  mat <- do.call(rbind, re_list)  # repetitions x 6
  
  Re_1_num <- colMeans(mat, na.rm = TRUE)
  
  Re_1 <- data.frame(
    n = n,
    cont_nonp = Re_1_num[1],
    cont_dp   = Re_1_num[2],
    atq_nonp  = Re_1_num[3],
    atq_dp    = Re_1_num[4],
    away_nonp = Re_1_num[5],
    away_dp   = Re_1_num[6]
  )
  
  Re <- rbind(Re, Re_1)
 
}

output_file <- file.path( output_dir, "Example4.csv" )
write.csv(Re,file = output_file,row.names = FALSE)

