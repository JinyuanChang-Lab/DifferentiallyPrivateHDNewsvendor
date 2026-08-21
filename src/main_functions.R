################################################################################
########################## 1. Useful tool functions ############################ 
################################################################################

#-------------------------------------------------------------------------------
#                           1.1 Kernel functions                         
#------------------------------------------------------------------------------- 

# Return: ``\bar{K}(x) - τ''  
kernel_cdf <- function(x, tau, kernel = "Laplacian", w = NULL) { 
  # """Gradient weights for the (weighted) smoothed check loss."""
  kernels <- list(
    Laplacian = function(x) 0.5 + 0.5 * sign(x) * (1 - exp(-abs(x))),                #  Laplace(0,1) CDF; closed form: 0.5 + 0.5 * sign(x) * (1 - exp(-|x|))
    Gaussian = function(x) pnorm(x),                                                 #  Standard normal CDF Φ(x)
    Logistic = function(x) 1 / (1 + exp(-x)),                                        #  Logistic CDF: 1 / (1 + e^{-x})
    Uniform = function(x) ifelse(abs(x) <= 1, 0.5 * (1 + x), as.numeric(x > 1)),     #  Linear (tent) CDF for Uniform kernel on [-1,1]; right tail is 1{x>1}
    Epanechnikov = function(x) ifelse(
      abs(x) <= sqrt(5),
      0.25 * (2 + 3 * x / sqrt(5) - (x / sqrt(5))^3),
      as.numeric(x > sqrt(5))
    )                                                                                #  Epanechnikov-based smooth CDF (normalized to have compact support on [-√5, √5])
  )
  # Validate kernel name
  if (!kernel %in% names(kernels)) {
    stop(paste("Unsupported kernel type:", kernel))
  }
  
  Ker <- kernels[[kernel]]                                                           #  pick the CDF-like smoothing function K
  result <- Ker(x) - tau                                                             #  gradient weight = K(x) - τ  
  if (!is.null(w)) {
    result <- w * result                           
  }                                                                                  #  optional sample weights
  return(result)                                                                     #  vector of weights, same length as x
}

#-------------------------------------------------------------------------------
#                           1.2  Clip functions                         
#------------------------------------------------------------------------------- 

# Return: ``w_{B}(X)''
clipping_inf_new <- function(X, B) {
  return(matrix(pmax(-B, pmin(X, B)), nrow = nrow(X), ncol = ncol(X)))  # Clip covariates elementwise to [-B, B]
}

#-------------------------------------------------------------------------------
#                         1.3  Hard Thresholding                         
#------------------------------------------------------------------------------- 

# Return: v_{S} (non-DP)
ht <- function(v, s) {
  d <- length(v)                                           # Get the length of vector v
  top_indices <- order(abs(v), decreasing = TRUE)[1:s]     # Indices of the top s elements with the largest absolute values
  v_S <- rep(0, d)                                         # Initialize a zero vector of length d
  v_S[top_indices] <- v[top_indices]                       # Keep the top s elements and set others to zero
  return(v_S)                                              # Return the hard-thresholded vector
}

#-------------------------------------------------------------------------------
#              1.4 Generate laplace r.v.s  with specified scale                        
#------------------------------------------------------------------------------- 

rlaplace <- function(n, scale) {
  u <- runif(n, 0, 1)                                          # Generate n uniform random numbers between 0 and 1
  return(-scale * sign(u-0.5) * log(1 - 2 * abs(u-0.5)))       # Apply the inverse CDF method to generate Laplace-distributed samples
}

#-------------------------------------------------------------------------------
#                    1.5 Noisy Hard Thresholding (NoisyHT)                        
#-------------------------------------------------------------------------------

# Return: v_{S} (DP)
noisyht <- function(v, s, epsilon, delta, lambda_scale) {
  d <- length(v)                                                               # Get the length of vector v                                                    
  S <- c()                                                                     # Initialize the selected index set
  
  for (i in 1:s) {
    w <- rlaplace(d,scale=lambda_scale*2*sqrt(5*s*log(1/delta))/epsilon)       # Generate Laplace noise for selection
    candidates <- setdiff(1:d, S)                                              # Indices not yet selected
    values <- abs(v[candidates]) + w[candidates]                               # Add noise to absolute values of v for differential privacy
    j_max <- candidates[which.max(values)]                                     # Select index with maximum noisy score
    S <- c(S, j_max)                                                           # Add the selected index to S
  }
  
  v_S <- rep(0, d)                                                             # Initialize a zero vector for the output                 
  noise <- rlaplace(d,scale=lambda_scale*2*sqrt(5*s*log(1/delta))/epsilon)     # Generate Laplace noise for output 
  
  for (j in S) {
    v_S[j] <- v[j] + noise[j]                                                  # Add Laplace noise to the selected entries
  }
  
  return(v_S)                                                                  # Return the noisy hard-thresholded vector
}

#-------------------------------------------------------------------------------
#                     1.6 Quantile (check) loss function                       
#-------------------------------------------------------------------------------

# Return: u * (τ - I(u < 0))
check_loss <- function(u, tau) {
  u * (tau - as.numeric(u < 0))  # Quantile (check) loss function: ρ_τ(u) = u * (τ - I(u < 0))
}


#################################################################################
########################## 2. Main functions ####################################
#################################################################################

#-------------------------------------------------------------------------------
#      2.1 non-private sparse smoothed quantile regression  (Algorithm 1)                     
#-------------------------------------------------------------------------------
QT_highdim <- function(X, Y, lr = 0.1, beta0 = NULL, tau = 0.5, w = 0.1, s = 5, T = 1e3, 
                       kernel = 'Gaussian',intercept = TRUE) {
  ##### Inputs ####
  # X         : Numeric matrix of dimension n*p (with itcp) or n*(p-1) (without itcp). Represents the feature matrix.
  # Y         : Numeric vector of length n. Response variable.
  # lr        : Numeric scalar. Learning rate (step size) for the iterative update.
  # beta0     : Numeric vector of length p. Initial values for the coefficient vector.
  # tau       : Numeric scalar in (0,1). Quantile level for estimation. 
  # w         : Numeric scalar. Smoothing parameter used in the loss 
  # s         : Numeric scalar. Sparisty level.
  # T         : Integer. Number of iterations to perform.
  # kernel    : Character. Kernel used for smoothness.
  # intercept : Logical. If TRUE, an intercept (column of ones) is included in X (n*(p-1)) and the first element of beta0 corresponds to the intercept.
  
  
  n = nrow(X)                                      # Number of observations
  
  if(intercept){
    X_new <- cbind(rep(1, n), X)                   # Add intercept column if required
  }else{
    X_new <- X                                     # Use original X if no intercept
  }
  
  p_new = ncol(X_new)                              # Number of predictors (including intercept if any)
  
  if (is.null(beta0)) {       
    beta0 <- rep(0, p_new)                         # Initialize coefficients if not provided 
  }
  
  beta_seq <- matrix(0,nrow = p_new,ncol=T+1)      # Store coefficient path over iterations
  beta_seq[, 1] <- beta0                           # Store initial coefficients
  res <- Y - X_new %*% beta0                       # Compute initial residuals
  
  count <- 0                                                                                 # Initialize iteration counter              
  while (count < T) {                                                                        # Main iterative loop
    beta0  <- beta0  - (lr / n) * t(X_new) %*% kernel_cdf(-res / w, tau, kernel)             # Gradient descent update using smoothed check loss
    beta0  <- ht(beta0, s = s)                                                               # Apply hard thresholding to enforce sparsity
    res <- Y - X_new %*% beta0                                                               # Update residuals
    beta_seq[, count + 1] <- beta0                                                           # Store current coefficients
    count <- count + 1                                                                       # Increase iteration count
  }
  return(list(beta=beta0,                               # Final estimated coefficients
              beta_seq=beta_seq[, 1:(count + 1)],       # Coefficient sequence over iterations
              residuals=res,                            # Final residuals
              niter=count))                             # Number of iterations performed
}


#-------------------------------------------------------------------------------
#           2.2 DP sparse smoothed quantile regression  (Algorithm 3)                  
#-------------------------------------------------------------------------------
noisyQT_highdim <- function(X,Y,epsilon, delta, lr = 1, beta0 = NULL, tau = 0.5,
                            w = 0.1, s = 5, T = 1e3, B = 1, 
                            kernel = "Gaussian", intercept = TRUE) {
  ##### Inputs ####
  # X         : Numeric matrix of dimension n*p (with itcp) or n*(p-1) (without itcp). Represents the feature matrix.
  # Y         : Numeric vector of length n. Response variable.
  # epsilon   : Numeric scalar. Privacy parameter.
  # delta     : Numeric scalar. Privacy parameter.
  # lr        : Numeric scalar. Learning rate (step size) for the iterative update.
  # beta0     : Numeric vector of length p. Initial values for the coefficient vector.
  # tau       : Numeric scalar in (0,1). Quantile level for estimation. 
  # w         : Numeric scalar. Smoothing parameter used in the loss 
  # s         : Numeric scalar. Sparisty level.
  # T         : Integer. Number of iterations to perform.
  # B         : Numeric scalar. Truncation parameter.
  # kernel    : Character. Kernel used for smoothness.
  # intercept : Logical. If TRUE, an intercept (column of ones) is included in X (n*(p-1)) and the first element of beta0 corresponds to the intercept.
  
  
  n = nrow(X)                                                 # Number of observations
  
  if(intercept){
    X_new <- cbind(rep(1, n), X)                              # Add intercept column if required
    trunc_X1 <- cbind(rep(1, n),clipping_inf_new(X, B))       # Clipped version of X for privacy
  }else{
    X_new <- X                                                # Use original X if no intercept
    trunc_X1 <- clipping_inf_new(X, B)                        # Clipped version of X for privacy
  }
  
  p_new = ncol(X_new)                                         # Number of predictors (including intercept if any)
  
  if (is.null(beta0)) {
    beta0 <- rep(0, p_new)                                    # Initialize coefficients if not provided 
  }
  
  beta_seq <- matrix(0, nrow = p_new, ncol = T + 1)           # Store coefficient path over iterations
  beta_seq[, 1] <- beta0                                      # Store initial coefficients
  res <- Y - X_new %*% beta0                                  # Compute initial residuals
  
  lambda_scale <- 2 * (lr / n) * B * (max(tau, 1 - tau))      # Noise scale for DP mechanism
  
  count <- 0                                                                                                 # Initialize iteration counter       
  while (count < T) {                                                                                        # Main iterative loop
    beta0  <- beta0  - (lr / n) * (t(trunc_X1) %*% kernel_cdf(-res / w, tau, kernel = kernel))               # Gradient update
    beta0  <- noisyht(beta0, s = s, epsilon = epsilon / T, delta = delta / T, lambda_scale = lambda_scale)   # Apply NoisyHT to enforce sparsity
    res <- Y - X_new %*% beta0                                                                               # Update residuals
    beta_seq[, count + 1] <- beta0                                                                           # Store current coefficients
    count <- count + 1                                                                                       # Increase iteration count
  }
  
  # Return the results
  return(list(
    beta = beta0,                                # Final estimated DP coefficients                        
    beta_seq = beta_seq[, 1:(count + 1)],        # Coefficient sequence over iterations
    residuals = res,                             # Final residuals
    niter = count                                # Number of iterations performed
  ))
}


#-------------------------------------------------------------------------------
#     2.3 Objective function based on the smoothed check loss for L-BFGS
#-------------------------------------------------------------------------------
loss_fn <- function(beta, X, Y, tau, w) {
  ##### Inputs ####
  # beta      : Numeric vector. Regression coefficients including intercept.
  # X         : Numeric matrix of dimension n*p. Feature matrix without intercept.
  # Y         : Numeric vector of length n. Response variable.
  # tau       : Numeric scalar in (0,1). Quantile level.
  # w         : Numeric scalar. smoothing parameter.
  
  beta <- as.numeric(beta)
  
  if (!is.numeric(w) || length(w) != 1 || w <= 0) {
    stop("w must be a positive scalar.")
  }
  if (!is.numeric(tau) || length(tau) != 1 || tau <= 0 || tau >= 1) {
    stop("tau must be a scalar in (0, 1).")
  }
  
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (length(beta) != p + 1) {
    stop("length(beta) must equal ncol(X) + 1.")
  }
  if (length(Y) != n) {
    stop("length(Y) must equal nrow(X).")
  }
  
  X_new <- cbind(1, X)
  res <- as.vector(Y - X_new %*% beta)
  u <- res / w
  
  mean(res * (pnorm(u) + tau - 1) + w * dnorm(u))
}


#-------------------------------------------------------------------------------
#      2.4 Gradient function based on the smoothed check loss  for L-BFGS
#-------------------------------------------------------------------------------
grad_fn <- function(beta, X, Y, tau, w) {
  ##### Inputs ####
  # beta      : Numeric vector. Regression coefficients including intercept.
  # X         : Numeric matrix of dimension n*p. Feature matrix without intercept.
  # Y         : Numeric vector of length n. Response variable.
  # tau       : Numeric scalar in (0,1). Quantile level.
  # w         : Numeric scalar.  smoothing parameter.
  
  beta <- as.numeric(beta)
  
  if (!is.numeric(w) || length(w) != 1 || w <= 0) {
    stop("w must be a positive scalar.")
  }
  if (!is.numeric(tau) || length(tau) != 1 || tau <= 0 || tau >= 1) {
    stop("tau must be a scalar in (0, 1).")
  }
  
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (length(beta) != p + 1) {
    stop("length(beta) must equal ncol(X) + 1.")
  }
  if (length(Y) != n) {
    stop("length(Y) must equal nrow(X).")
  }
  
  X_new <- cbind(1, X)
  res <- as.vector(Y - X_new %*% beta)
  grad <- t(X_new) %*% (pnorm(-res / w) - tau) / n
  
  as.numeric(grad)
}




