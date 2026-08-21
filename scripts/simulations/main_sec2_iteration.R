################################################################################
################################# Required packages ############################ 
################################################################################
library(here) 
library(ggplot2)
library(dplyr)
library(grid)
library(patchwork)


################################################################################
################################# 1. Main functions ############################ 
################################################################################


#-------------------------------------------------------------------------------
#             1.1 Gradient descent under smoothed check loss                          
#-------------------------------------------------------------------------------

GD_smooth <- function(X ,Y, lr_sm = 1, beta0 = NULL, tau = 0.5, varpi = 0.1, T = 1000,  intercept = TRUE) {
  ##### Inputs ####
  # X         : Numeric matrix of dimension n*p (with itcp) or n*(p-1) (without itcp). Represents the feature matrix.
  # Y         : Numeric vector of length n. Response variable.
  # lr_sm     : Numeric scalar. Learning rate (step size) for the iterative update.
  # beta0     : Numeric vector of length p. Initial values for the coefficient vector.
  # tau       : Numeric scalar in (0,1). Quantile level for estimation. 
  # varpi     : Numeric scalar. Smoothing parameter used in the loss function.
  # T         : Integer. Number of iterations to perform.
  # intercept : Logical. If TRUE, an intercept (column of ones) is included in X (n*(p-1)) and the first element of beta0 corresponds to the intercept.
  
  
  n = nrow(X)                                                    # number of observations
  if(intercept){                                                 # whether to include an intercept
    X_new <- cbind(1, X)                                         # prepend a column of ones to X
  }else{                                                         # otherwise
    X_new <- X                                                   # use X as is
  }
  p = ncol(X_new)                                                # number of coefficients to estimate (incl. intercept if any)
  
  if (is.null(beta0)) {                                          # if no initial beta provided
    beta0 <- rep(0, p)                                           # initialize all coefficients to zero
    if (intercept) beta0[1] <- mean(Y)                           # set intercept start at sample mean of Y
  }
  beta0_sm <- beta0                                              # working copy of current iterate for smoothing GD
  
  beta_seq_sm <- matrix(0, nrow = p, ncol = T + 1)               # container to store all iterates (columns)
  beta_seq_sm[, 1] <- beta0_sm                                   # store the initial iterate
  
  res_sm <- Y - X_new %*% beta0_sm                               # initial residuals r = Y - X_new %*% beta
  
  
  # Gradient weights for the (weighted) smoothed check loss.
  kernel_cdf <- function(x, tau) { 
    Ker <-  function(x){pnorm(x)}                                
    return(Ker(x) - tau)
  }                                                              # Gaussian kernel 
  
  # main iteration loop
  for (count in 1:T) {                                                                          # iterate T times
    beta0_sm <- beta0_sm - (lr_sm / n) * (t(X_new) %*% kernel_cdf(-res_sm / varpi, tau))        # smoothed GD update 
    res_sm <- Y - X_new %*% beta0_sm                                                            # update residuals  
    beta_seq_sm[, count + 1] <- beta0_sm                                                        # store current iterate
  }
  
  return(list(                                                                                  # return a list with the trajectory and iteration count
    beta_seq_sm = beta_seq_sm[, 1:(T + 1)],                                                     # matrix of iterates (p × (T+1))
    niter = T                                                                                   # total iterations performed
  ))
}


#-------------------------------------------------------------------------------
#             1.2 Gradient descent under standard check loss                          
#-------------------------------------------------------------------------------

GD_standard <- function(X,Y, lr_ns = 1, beta0 = NULL, tau = 0.5,  T = 1000,  intercept = TRUE) {
  ##### Inputs ####
  # X         : Numeric matrix of dimension n*p (with itcp) or n*(p-1) (without itcp). Represents the feature matrix.
  # Y         : Numeric vector of length n. Response variable.
  # lr_ns     : Numeric scalar. Learning rate (step size) for the iterative update.
  # beta0     : Numeric vector of length p. Initial values for the coefficient vector.
  # tau       : Numeric scalar in (0,1). Quantile level for estimation.  
  # T         : Integer. Number of iterations to perform.1
  # intercept : Logical. If TRUE, an intercept (column of ones) is included in X (n*(p-1)) and the first element of beta0 corresponds to the intercept.
  
  n = nrow(X)                                       # number of observations      
  if(intercept){                                    # whether to include an intercept
    X_new <- cbind(rep(1, n), X)                    # prepend a column of ones to X  
  }else{                                            # otherwise
    X_new <- X                                      # use X as is
  }
  p = ncol(X_new)                                   # number of coefficients to estimate (incl. intercept if any)
  
  if (is.null(beta0)) {                             # if no initial beta provided
    beta0 <- rep(0, p)                              # initialize all coefficients to zero
    if (intercept) beta0[1] <- mean(Y)              # set intercept start at sample mean of Y 
  }
  
  
  beta0_st <- beta0                                 # working copy of current iterate for standard GD
  beta_seq_st <- matrix(0, nrow = p, ncol = T + 1)  # container to store all iterates (columns)
  beta_seq_st[, 1] <- beta0_st                      # store the initial iterate
  res_st <- Y - X_new %*% beta0_st                      # initial residuals r = Y - X_new %*% beta
  
  # quantile_weight function
  quantile_weight <- function(x, tau) {
    z <- numeric(length(x))                                                      # allocate weight vector
    z[x > 0] <- tau                                                              # ψ_tau(r) = τ for r > 0   
    z[x < 0] <- tau - 1                                                          # ψ_tau(r) = τ−1 for r < 0
    zero_idx <- which(x == 0)                                                    # subgradient set at r = 0 is [τ−1, τ]
    if (length(zero_idx) > 0) { 
      z[zero_idx] <- runif(length(zero_idx), min = tau - 1, max = tau)
    }                                                                            # draw uniformly in [τ−1, τ] to pick a valid subgradient
    return(-z)
  }
  
  for (count in 1:T) {                                                               # iterate T times
    beta0_st <- beta0_st - (lr_ns / n) * (t(X_new) %*% quantile_weight(res_st, tau)) # standard GD update
    res_st <- Y - X_new %*% beta0_st                                                 # update residuals with the new β       
    beta_seq_st[, count + 1] <- beta0_st                                             # store current iterate
  }
  
  return(list(                                                                   # return a list with the trajectory and iteration count
    beta_seq_st = beta_seq_st[, 1:(T + 1)],                                      # matrix of iterates (p × (T+1))
    niter = T                                                                    # total iterations performed
  ))
}



################################################################################
################################# 2. comparison results ######################## 
################################################################################

###### Basic settings ######
tau = 0.7                                 # tau =b , b+h=1
n <- 500                                  # Sample size
p <- 19                                   # Ture dimension is 19+1=20 
lr_sm <- 2                                # Learning rate for smooth check loss
lr_ns_1 <- 2                              # Learning rate 1 for standard check loss
lr_ns_2 <- 0.5                            # Learning rate 1 for standard check loss
T <- 90                                   # Number of iterations
varpi <-  ((p+1 + log(n)) / n)^0.4        # Smoothing parameter

#-------------------------------------------------------------------------------
#                           2.1 Gaussian setting                         
#-------------------------------------------------------------------------------

set.seed(3)                                                       # Set seed for reproducibility
itcp = 2                                                          # intercept = 2
beta_0 <-  rep(1 , p) * (2 * rbinom(p, 1, 0.5) - 1)       # beta: {1, -1, -1, 1 ...}  
X <- matrix(rnorm(n * p ), nrow = n, ncol = p)                    # Generate X ~ N(0,1)
sd_noise = sqrt(2)                                                # the standard deviation of Gaussian noise
Y <- itcp + X %*% beta_0 +   rnorm(n,mean = 0, sd = sd_noise)     # Generate Y
beta_true <- c(itcp, beta_0) + c(sd_noise*qnorm(tau) ,rep(0,p))   # true beta^*
beta_norm <- sqrt(sum(beta_true^2))                               # The L2 norm of true beta
# Compute the results 
out_smooth <- GD_smooth(X,Y, lr_sm = lr_sm, beta0 = NULL, tau = tau, varpi = varpi, T = T , intercept = TRUE)     # smoothed GD with learning rate 2
out_standard_1 <- GD_standard(X,Y, lr_ns = lr_ns_1, beta0 = NULL, tau = tau, T = T , intercept = TRUE)            # standard GD with learning rate 2
out_standard_2 <- GD_standard(X,Y, lr_ns = lr_ns_2, beta0 = NULL, tau = tau, T = T , intercept = TRUE)            # standard GD with learning rate 0.5
smooth_G <- sqrt(colSums((out_smooth$beta_seq_sm - beta_true) ^ 2)) / beta_norm                                   # sequence of the relative l2 error
nonsmooth_G1 <- sqrt(colSums((out_standard_1$beta_seq_st - beta_true) ^ 2)) / beta_norm                           # sequence of the relative l2 error
nonsmooth_G2 <- sqrt(colSums((out_standard_2$beta_seq_st - beta_true) ^ 2)) / beta_norm                           # sequence of the relative l2 error

#-------------------------------------------------------------------------------
#                           2.2 Uniform setting                         
#-------------------------------------------------------------------------------

set.seed(3)                                                                  # Set seed for reproducibility
itcp = 2                                                                     # intercept = 2
beta_0 <- rep(1 , p) * (2 * rbinom(p, 1, 0.5) - 1)                  # beta: {1, -1, -1, 1 ...}  
X <- matrix(runif(n * p, min = -2, max = 2), nrow = n, ncol = p)             # Generate X ~ U[-2,2]
Y <- itcp + X %*% beta_0 +  (2.5+X[, 1]) *  rnorm(n)                         # Generate Y
beta_true <- c(itcp, beta_0) + c(c(2.5*qnorm(tau), qnorm(tau)) ,rep(0,p-1))  # true beta^*
beta_norm <- sqrt(sum(beta_true^2))                                          # The L2 norm of true beta  
# Compute the results 
out_smooth_U <- GD_smooth(X,Y, lr_sm = lr_sm, beta0 = NULL, tau = tau, varpi = varpi, T = T , intercept = TRUE)   # smoothed GD with learning rate 2
out_standard_U_1 <- GD_standard(X,Y, lr_ns = lr_ns_1, beta0 = NULL, tau = tau, T = T , intercept = TRUE)          # standard GD with learning rate 2
out_standard_U_2 <- GD_standard(X,Y, lr_ns = lr_ns_2, beta0 = NULL, tau = tau, T = T , intercept = TRUE)          # standard GD with learning rate 0.5
smooth_U <- sqrt(colSums((out_smooth_U$beta_seq_sm - beta_true) ^ 2)) / beta_norm                                 # sequence of the relative l2 error
nonsmooth_U1 <- sqrt(colSums((out_standard_U_1$beta_seq_st - beta_true) ^ 2)) / beta_norm                         # sequence of the relative l2 error
nonsmooth_U2 <- sqrt(colSums((out_standard_U_2$beta_seq_st - beta_true) ^ 2)) / beta_norm                         # sequence of the relative l2 error 


################################################################################
################################# 3. comparison figure  ######################## 
################################################################################



## Color mapping and order 
color_order  <- c('Smoothed (learning rate = 2)', 'Standard (learning rate = 2)', 'Standard (learning rate = 0.5)')
color_values <- c('Smoothed (learning rate = 2)'   = '#4B3F72',  'Standard (learning rate = 2)'   = 'tomato',  'Standard (learning rate = 0.5)' = 'deepskyblue')

## Build data for the two panels  
niter <- 1:(T + 1)
# Panel 1: x ~ N(0,1)
G_combined_data <- rbind(
  data.frame(iteration = niter, error = as.vector(log(smooth_G[niter])), type = 'Smoothed (learning rate = 2)'),
  data.frame(iteration = niter, error = as.vector(log(nonsmooth_G1[niter])), type = 'Standard (learning rate = 2)'),
  data.frame(iteration = niter, error = as.vector(log(nonsmooth_G2[niter])), type = 'Standard (learning rate = 0.5)')
)
G_combined_data_sub <- G_combined_data %>% filter(iteration >= 50 & iteration <= 80) # subplot

U_combined_data <- rbind(
  data.frame(iteration = niter, error = as.vector(log(smooth_U[niter])), type = 'Smoothed (learning rate = 2)'),
  data.frame(iteration = niter, error = as.vector(log(nonsmooth_U1[niter])), type = 'Standard (learning rate = 2)'),
  data.frame(iteration = niter, error = as.vector(log(nonsmooth_U2[niter])), type = 'Standard (learning rate = 0.5)')
)
U_combined_data_sub <- U_combined_data %>% filter(iteration >= 35 & iteration <= 75) # subplot

## Two-level theme function 
theme_comb <- function(a){
  if(a==1){
    theme(
      text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
    )
  }else{
    theme(
      legend.position = "none",
      text = element_text(size = 10),
      axis.text = element_text(size = 10),
      axis.line.x = element_line(color = "black", size = 0.8),
      axis.line.y = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(0.15, "cm")
    )
  }
}

## Legend title expression
legend_title_expr <- expression(frac(bold(b), bold(b) ~ bold("+") ~ bold(h)) ~ bold("=") ~ bold("0.7"))
line_size_main <- 0.8

## Panel 1 main plot + inset
main_plot_1 <- ggplot(G_combined_data, aes(x = iteration, y = error, color = type)) +
  geom_line(size = line_size_main) +
  labs(x = 'Number of iterations', y = 'Estimation error', title = expression(x %~% N(0, 1))) +
  scale_color_manual(name = legend_title_expr, values = color_values,
                     breaks = color_order, limits = color_order) +
  theme_minimal() + theme_comb(1)

G_inset_plot <- ggplot(G_combined_data_sub, aes(x = iteration, y = error, color = type)) +
  geom_line(size = 1) +
  labs(x = 'Number of iterations', y = expression('log relative'~L[2]~'-error')) +
  scale_color_manual(values = color_values, breaks = color_order, limits = color_order) +
  theme_void() + theme_comb(2)

main_plot_1 <- main_plot_1 +
  annotation_custom(ggplotGrob(G_inset_plot), xmin = 40, xmax = 90, ymin = -2.1, ymax = -0.4) +
  annotation_custom(rectGrob(gp = gpar(lwd = 1, col = "black", lty = 2, fill = NA)),
                    xmin = 50, xmax = 80, ymin = -2.9, ymax = -2.6) +
  annotate("segment", x = 60, y = -2.4, xend = 62, yend = -2.2,
           arrow = arrow(type = "closed", length = unit(0.08, "inches")),
           size = 0.5, color = "black")

## Panel 2 main plot + inset
main_plot_2 <- ggplot(U_combined_data, aes(x = iteration, y = error, color = type)) +
  geom_line(size = line_size_main) +
  labs(x = 'Number of iterations', y = 'Estimation error', title = expression(x %~% U * "[-2, 2]")) +
  scale_color_manual(name = legend_title_expr, values = color_values,
                     breaks = color_order, limits = color_order) +
  theme_minimal() + theme_comb(1)

U_inset_plot <- ggplot(U_combined_data_sub, aes(x = iteration, y = error, color = type)) +
  geom_line(size = 1) +
  labs(x = 'Number of iterations', y = expression('log relative'~L[2]~'-error')) +
  scale_color_manual(values = color_values, breaks = color_order, limits = color_order) +
  theme_void() + theme_comb(2)

main_plot_2 <- main_plot_2 +
  annotation_custom(ggplotGrob(U_inset_plot), xmin = 40, xmax = 90, ymin = -2.0, ymax = -0.4) +
  annotation_custom(rectGrob(gp = gpar(lwd = 1, col = "black", lty = 2, fill = NA)),
                    xmin = 35, xmax = 75, ymin = -2.7, ymax = -2.3) +
  annotate("segment", x = 45, y = -2.2, xend = 47, yend = -2.0,
           arrow = arrow(type = "closed", length = unit(0.08, "inches")),
           size = 0.5, color = "black")

## Final combined figure (two panels side by side; shared legend)
final_plot_07 <- (main_plot_1 | main_plot_2) + plot_layout(guides = "collect")
#print(final_plot_07)

## Save the results
here::i_am("scripts/simulations/main_sec2_iteration.R")
figure_file <- here::here( "results", "figures", "iteration.pdf" )
ggsave(filename = figure_file, plot = final_plot_07, width = 14, height = 4, units = "in", dpi = 300)
 