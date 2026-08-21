############################################################################### 
###########################  1. Preparation    ################################ 
############################################################################### 

#-------------------------------------------------------------------------------
#                         1.1 Required packages                         
#------------------------------------------------------------------------------- 
library(here)  
library(quantregForest)
library(qrnn)
library(foreach)
library(doParallel)
library(doRNG)
library(doSNOW)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggtext)

#-------------------------------------------------------------------------------
#                        1.2 Load functions and data                       
#-------------------------------------------------------------------------------
here::i_am("scripts/real_data/supp_NN_real_data.R")
source( here::here("src", "supp_functions_NN.R") )
processed_data <- readRDS( here::here("data", "processed", "processed_data.rds") )
RD_X <- processed_data$RD_X
RD_Y <- processed_data$RD_Y


################################################################################## 
############################### 2. Numerical performance #########################  
################################################################################## 
## low dim
RD_X_low <- RD_X[,1:4]
RD_X_low <- cbind(1,RD_X_low)
colnames(RD_X_low)[1] <- "itcp" 

nv_loss <- function(q, d, b, h) {  mean(h * pmax(q - d, 0) + b * pmax(d - q, 0)) }

RD_NN <- function(X, Y,  train_indices, test_indices,
                     b=b, h=h, c_w = 0.1, lr = 0.1, T_t = 100,  # regural parameters
                     m = 200, sd_theta0 = 0.1, theta0 = NULL  ) {
  
  # train data and test data
  X_train <- X[train_indices, , drop = FALSE]
  X_test  <- X[test_indices,  , drop = FALSE]
  Y_train <- Y[train_indices]
  Y_test  <- Y[test_indices] 
  
  # dimension and training sample size
  p_tr <- ncol(X_train)
  n_tr <- nrow(X_train)
  
  w_dnn <- c_w * ((p_tr + log(n_tr)) / n_tr)^(0.4) # smooth parameter
  if (is.null(theta0)) { theta0 <- rnorm( (m * (p_tr + 1) +1), sd = sd_theta0) } # initial theta value
  theta0 <- clipping_inf(theta0, 1000)
  
  ## dnn_np
  tau <- b/(b+h)
  dnn_out_np <- fit_nn_np(X = X_train, Y = Y_train, eta0 = lr, 
                           T_all = T_t, tau = tau, varpi = w_dnn, theta0 = theta0)
  f_hat_np <- predict_nn(dnn_out_np, X_test)
  COST_np <- nv_loss(f_hat_np, Y_test, b, h) 
  
  ## SAA
  f_hat_saa <- as.numeric(quantile(Y_train, probs = tau, type = 1))
  COST_saa  <- nv_loss(f_hat_saa, Y_test, b, h)   
  
  ## quantregForest
  fit_qrf <- quantregForest(x = X_train, y = Y_train)
  f_hat_qrf <- predict(fit_qrf, X_test, what = tau)
  COST_qrf  <- nv_loss(f_hat_qrf, Y_test, b, h)   
  
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
  COST_qrnn  <- nv_loss(f_hat_qrnn, Y_test, b, h)    
  
  return(c(COST_np,COST_saa,COST_qrf,COST_qrnn))
}





n_0 <- nrow(RD_X_low)
repetitions <- 500
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
                  .packages = c("quantregForest","qrnn")) %dorng% {    
                    h     <- 1
                    b_vec <- seq(2 ,18, by = 2)    
                    train_indices <- sample.int(n_0, size = floor(0.7 * n_0))
                    test_indices  <- setdiff(seq_len(n_0), train_indices)
                    
                    sapply(b_vec, function(b)
                      RD_NN(X=RD_X_low, Y=RD_Y,  train_indices=train_indices, test_indices=test_indices,b = b,h=h))
                  }
close(pb)
stopCluster(cl)
end <- Sys.time()
print(end-start)   ## 22 mins

b_vec <- seq(2 ,18, by = 2)    
k <- length(b_vec) 
arr <- array(test_1, dim = c( 4 , k, repetitions))     
test_1_mean  <- apply(arr, c(1, 2), mean, na.rm = TRUE) 
test_1_mean

FNN <-   (test_1_mean[1,])
SAA <-   (test_1_mean[2,])
QRF <-   (test_1_mean[3,]) 
QRNN <-   (test_1_mean[4,]) 

plot_df_1 <- data.frame( b = b_vec, FNN = FNN, SAA = SAA, QRF = QRF,  QRNN = QRNN)
plot_df_long <- plot_df_1 %>% pivot_longer(cols = c("FNN", "SAA","QRF", "QRNN"), names_to = "Method",values_to = "Cost") #library(tidyr) 
method_colors <- c(  "SAA" = "deepskyblue","FNN" = "tomato","QRNN" = "#4B3F72", "QRF" = "#2ca02c")

figure1 <- ggplot(plot_df_long, aes(x = b, y = Cost, color = Method)) +
  geom_line(size = 1, linetype = "solid") +
  geom_point(size = 2) +
  labs(
    x = "b (Underage cost)",
    y = "Average out-of-sample cost",
    color = "Method"
  ) +
  scale_color_manual(
    values = method_colors,
    breaks = c("FNN", "QRF", "QRNN", "SAA"),
    labels = c("FNN", "QRF", "QRNN", "SAA")
  ) +
  scale_x_continuous(breaks = b_vec) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(size = 1, color = "black"),
    axis.title = element_text(face = "bold", size = 14, color = "black"),
    axis.text = element_text(face = "bold", size = 12, color = "black"),
    legend.title = element_text(face = "bold", size = 12, color = "black"),
    legend.text = element_markdown(face = "bold", size = 12, color = "black"),
    legend.position = "right",
    plot.title = element_text(face = "bold", color = "black"),
    plot.subtitle = element_text(face = "bold", color = "black"),
    strip.text = element_text(face = "bold", color = "black")
  ) +
  guides(color = guide_legend(title = "Method"))

#figure1
 
## Save the results  
figure_file_RD_NN <- here::here( "results", "figures", "figure_NN_RD.pdf" )
ggsave(filename = figure_file_RD_NN, plot = figure1, width = 5, height = 3.5, units = "in", dpi = 300)
