############################################################################### 
###########################  1. Preparation    ################################ 
############################################################################### 

#-------------------------------------------------------------------------------
#                         1.1 Required packages                         
#------------------------------------------------------------------------------- 
library(here) 
library(quantreg)
library(dplyr)
library(tidyr)

#-------------------------------------------------------------------------------
#                        1.2 Load functions and data                       
#-------------------------------------------------------------------------------
here::i_am("scripts/simulations/supp_runtime_comparison.R")
source(here::here("src", "main_functions.R"))   


################################################################################# 
############################  2. Execution algorithm  ###########################  
################################################################################# 

prepare_once <- function(p, n, b = 2, h = 1, s_star = 5, epsilon = 0.5, delta = 0.0001) {
  set.seed(2025)
  
  tau <- b / (b + h)
  p_X <- p - 1
  itcp <- 4
  n_1 <- (5 * n) / 4
  n_all <- n_1 + 30000
  sigma <- sqrt(2)
  
  beta_star <- matrix( c(sample(c(2, -2), s_star - 1, replace = TRUE), rep(0, times = p - s_star)), ncol = 1 )
  
  X <- matrix(rnorm(n_all * p_X), nrow = n_all)
  Noise <- rnorm(n_all, mean = -sigma * qnorm(tau), sd = sigma)
  Y <- as.numeric(itcp + X %*% beta_star + Noise) 
  
  test_indices <- sample.int(n_all, size = 30000, replace = FALSE)
  remain <- setdiff(seq_len(n_all), test_indices)
  init_indices <- sample(remain, size = n / 4, replace = FALSE)
  train_indices <- setdiff(remain, init_indices)
  
  X_init  <- X[init_indices, , drop = FALSE]
  X_train <- X[train_indices, , drop = FALSE]
  X_test  <- X[test_indices, , drop = FALSE]
  Y_init  <- Y[init_indices]
  Y_train <- Y[train_indices]
  Y_test  <- Y[test_indices] 
  X_train_new <- cbind(1, X_train)
  X_test_new  <- cbind(1, X_test) 
  
  p_tr <- ncol(X_train) + 1
  n_tr <- nrow(X_train)
  s_vec <- seq(2, 10, by = 1)
  
  c_w <- 0.5
  lr_1 <- 2
  T_1 <- 200
  lr_2 <- 0.1
  T_2 <- 4 * ceiling(log(n_tr / log(p_tr)))
  n_init <- nrow(X_init)
  
  BIC_vec <- c()
  beta_mat <- NULL
  
  for (s in s_vec) {
    w_pre <- c_w * ((s * log(p_tr) + log(n_init)) / n_init)^0.4
    QT_pre <- QT_highdim( X_init, Y_init, lr = lr_1, beta0 = NULL, tau = tau, w = w_pre, s = s, T = T_1 )
    Q_tau <- sum(tau * QT_pre$residuals * (QT_pre$residuals >= 0) +
                   (tau - 1) * QT_pre$residuals * (QT_pre$residuals < 0))
    BIC_vec <- c(BIC_vec, log(Q_tau) + s * log(n_init) * log(p_tr) / (2 * n_init))
    beta_mat <- rbind(beta_mat, QT_pre$beta)
  }
  s_select <- s_vec[which.min(BIC_vec)]
  w_select <- c_w * ((s_select * log(p_tr) + log(n_tr)) / n_tr)^0.4
  beta0_init <- beta_mat[which.min(BIC_vec), ]
  
  c_gamma = 0.05
  B_QT <- c_gamma * sqrt(log(p_tr) + log(n_tr))
  
  list(
    tau = tau, b = b, h = h,
    X_train = X_train, Y_train = Y_train,
    X_train_new = X_train_new, X_test_new = X_test_new, 
    p_tr = p_tr, w_select = w_select,
    s_select = s_select, beta0_init = beta0_init,
    lr_2 = lr_2, T_2 = T_2, B_QT= B_QT,
    epsilon = epsilon, delta = delta
  )
}


run_SQR <- function(dat) {
  QT_new <- QT_highdim(
    dat$X_train, dat$Y_train,
    lr = dat$lr_2,
    beta0 = dat$beta0_init,
    tau = dat$tau,
    w = dat$w_select,
    s = dat$s_select,
    T = dat$T_2 
  ) 
  return(QT_new)
}

run_DPSQR <- function(dat) {
  DPQT_new <- noisyQT_highdim(
    dat$X_train, dat$Y_train,
    dat$epsilon, dat$delta,
    lr = dat$lr_2,
    beta0 = dat$beta0_init,
    tau = dat$tau,
    w = dat$w_select,
    s = dat$s_select,
    T = dat$T_2,
    B = dat$B_QT
  ) 
  return(DPQT_new)
} 

run_ERM <- function(dat) {
  fit <-   rq.fit(
    x = dat$X_train_new,
    y = dat$Y_train,
    tau = dat$tau
  )
  return(fit)
}

n_vec <- c(1000, 2000)
p_vec <- c(25,50,100,200,400 )
n_rep <- 100 

time_res <- NULL

for (n in n_vec) {
  for (p in p_vec) {
    cat("Running n =", n, ", p =", p, "\n")
    
    dat <- prepare_once(p = p, n = n)
    
    t_sqr <- numeric(n_rep)
    t_dpsqr <- numeric(n_rep)
    t_erm <- numeric(n_rep)
    
    for (r in 1:n_rep) {
      t_sqr[r]   <- system.time(run_SQR(dat))[["elapsed"]]
      t_dpsqr[r] <- system.time(run_DPSQR(dat))[["elapsed"]]
      t_erm[r]   <- system.time(run_ERM(dat))[["elapsed"]]
    }
    
    tmp <- data.frame(
      n = n,
      p = p,
      Method = c("l0_ERM", "DP_l0_ERM", "ERM"), 
      mean_s   = round(c(mean(t_sqr), mean(t_dpsqr), mean(t_erm)), 4)
    )
    
    time_res <- rbind(time_res, tmp)
  }
}

time_res_print <- time_res
time_res_print$mean_s   <- sprintf("%.3f", time_res_print$mean_s)


time_mean_wide <- time_res_print |>
  select(n, p, Method, mean_s) |>
  mutate(  Method = factor( Method, levels = c("l0_ERM", "DP_l0_ERM", "ERM") ) ) |>
  pivot_wider(names_from  = Method,values_from = mean_s) |>
  arrange(n, p)
 
time_mean_export <- time_mean_wide |>
  mutate(
    across(
      c(l0_ERM, DP_l0_ERM, ERM),
      \(x) formatC(x, format = "f", digits = 3)
    )
  )
 

# save the table
dir.create(here::here("results", "tables"),recursive = TRUE,showWarnings = FALSE)
table_file <- here::here( "results", "tables",  "runtime_comparison.csv")
write.csv( time_mean_export, file = table_file,  row.names = FALSE,quote = FALSE )



 