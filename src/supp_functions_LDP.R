################################################################################
########################## 1. Useful tool functions ############################ 
################################################################################

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

# Return: ``w_{B}(X)''
clipping_inf_new <- function(X, B) {
  return(matrix(pmax(-B, pmin(X, B)), nrow = nrow(X), ncol = ncol(X)))  # Clip covariates elementwise to [-B, B]
}

rlaplace <- function(n, scale) {
  u <- runif(n, 0, 1)                                          # Generate n uniform random numbers between 0 and 1
  return(-scale * sign(u-0.5) * log(1 - 2 * abs(u-0.5)))       # Apply the inverse CDF method to generate Laplace-distributed samples
}

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


################################################################################
##############     1. function for LDP (Algorithm A1)   ########################  
################################################################################ 

noisyQT_highdim_ldp <- function(X,Y,epsilon, delta, beta0, lr = 1, tau = 0.5,
                                w = 0.1, s = 5, T = 1e3, B = 1, kernel = "Gaussian", intercept = TRUE){ 
  
  # ----- Construct design matrix and clipped version -----
  n <- nrow(X)
  if (intercept) {
    X_new <- cbind(rep(1, n), X)
    trunc_X1 <- cbind(rep(1, n), clipping_inf_new(X, B))
  } else {
    X_new <- X
    trunc_X1 <- clipping_inf_new(X, B)
  }
  p_new <- ncol(X_new)
   
  # ----- Fixed support set S0 -----
  supp0 <- which(beta0 != 0) 
  S0 <- sort(unique(supp0))
  s0 <- length(S0)
  
  beta_seq <- matrix(0, nrow = p_new, ncol = T + 1)
  beta_seq[, 1] <- beta0
  
  # ----- Per-iteration LDP budget -----
  eps_it <- epsilon / T
  del_it <- delta / T
  bar_tau <- max(tau, 1 - tau)
  B_ldp <- 2 * B * bar_tau * sqrt(s0)
  sigma_it <- B_ldp * sqrt(2 * log(1.25 / del_it)) / eps_it
  res <- as.numeric(Y - X_new %*% beta0)
  
  # ----- Main iterative loop -----
  count <- 0                                                                                                 # Initialize iteration counter       
  while (count < T)  { 
    # psi_i = K(-(res_i)/w) - tau
    psi <- as.numeric(kernel_cdf(-res / w, tau, kernel = kernel))  # length n
    
    # Individual gradients restricted to S0: g_i,S0 = trunc_X1[i,S0] * psi_i
    XiS0 <- trunc_X1[, S0, drop = FALSE]          # n x s0
    G_S0 <- XiS0 * psi                             # n x s0 (row-wise scaling)
    
    # Local Gaussian perturbation (each user adds iid N(0, sigma_it^2) entries)
    W <- matrix(stats::rnorm(n * s0, mean = 0, sd = sigma_it), nrow = n, ncol = s0)
    G_tilde_S0 <- G_S0 + W
    
    # Aggregate perturbed gradients and update beta on S0 only
    grad_hat_S0 <- colMeans(G_tilde_S0)            # s0 vector = (1/n) sum_i tilde g_i,S0
    beta0[S0] <- beta0[S0] - lr * grad_hat_S0      # beta^{t+1} = beta^t - eta0 * (1/n) sum_i tilde g_i
    
    # Keep fixed support explicitly
    beta0[-S0] <- 0
    
    beta_seq[, count + 1] <- beta0
    res <- as.numeric(Y - X_new %*% beta0)
    count <- count + 1
  }
  
  return(list(
    beta = beta0,                                # Final estimated DP coefficients                        
    beta_seq = beta_seq[, 1:(count + 1)],        # Coefficient sequence over iterations
    residuals = res,                             # Final residuals
    niter = count                                # Number of iterations performed
  ))
}


################################################################################
#################   1. function for HDP (Algorithm A2)    ###################### 
################################################################################ 

noisyQT_highdim_hdp <- function(X, Y, epsilon, delta, K, lr = 1, beta0 = NULL, tau = 0.5,
                                w = 0.1, s = 5, T = 1000, B = 1, kernel = "Gaussian", intercept = TRUE) { 
  
  n <- nrow(X)  
  if (intercept) {
    X_new <- cbind(rep(1, n), X)
    trunc_X1 <- cbind(rep(1, n), clipping_inf_new(X, B))
  } else {
    X_new <- X
    trunc_X1 <- clipping_inf_new(X, B)
  }
  p_new <- ncol(X_new)   
  
  # ----- Initialize beta -----
  if (is.null(beta0)) {
    beta0 <- rep(0, p_new) 
  }
  
  # ----- Build disjoint groups -----
  perm <- sample.int(n)
  gid <- rep(1:K, length.out = n)
  gid <- gid[order(order(perm))]  # shuffle labels according to perm
  G_list <- split(seq_len(n), gid)
  
  # Basic sanity checks: disjointness and coverage
  all_idx <- sort(unlist(G_list, use.names = FALSE)) 
  
  # ----- Storage -----
  beta_seq <- matrix(0, nrow = p_new, ncol = T + 1)  # Store beta at the end of each outer iteration
  beta_seq[, 1] <- beta0
  
  # Per-outer-iteration privacy budget (parallel composition across disjoint groups within each sweep)
  eps_it <- epsilon / T
  del_it <- delta / T
  
  # ----- Main loop: outer iterations -----
  count <- 0 
  while (count < T) {
    # Within-iteration sweep over groups (k = 1,...,K)
    for (k in seq_len(K)) {
      idx <- G_list[[k]]
      nk <- length(idx)
      
      # Current group data
      Xg <- X_new[idx, , drop = FALSE]
      trunc_Xg <- trunc_X1[idx, , drop = FALSE]
      Yg <- Y[idx]
      
      # Residuals on the group
      res_g <- Yg - Xg %*% beta0
      
      # Group-wise gradient update
      beta0 <- beta0 - (lr / nk) * (t(trunc_Xg) %*% kernel_cdf(-res_g / w, tau, kernel = kernel))
      
      # Group-specific DP scale (matches your CDP scaling but with nk in place of n)
      lambda_scale_k <- 2 * (lr / nk) * B * (max(tau, 1 - tau))
      
      # Apply NoisyHT at each group update (fresh noise each time)
      beta0 <- noisyht(beta0, s = s, epsilon = eps_it, delta = del_it, lambda_scale = lambda_scale_k)
      
    }
    
    # End of sweep: record outer-iteration iterate
    count <- count + 1
    beta_seq[, count + 1] <- beta0 
  }
  
  # Final residuals on full data
  res <- Y - X_new %*% beta0
  
  return(list(
    beta = beta0,                                  # Final estimated DP coefficients
    beta_seq = beta_seq,                           # Coefficient sequence over outer iterations 
    residuals = res,                               # Final residuals
    niter = T,                                     # Number of outer iterations
    K = K,                                         # Number of groups
    groups = G_list                                # Group partition used
  ))
}

 
 