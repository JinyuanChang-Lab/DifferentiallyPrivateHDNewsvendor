####################
### ---------------  
### 1) Utilities  
### --------------- 
#################### 

# kernel used
barK_gaussian <- function(u) pnorm(u)

# l2-norm of a given vector
l2norm <- function(v) sqrt(sum(v^2))

# vector to list
unpack_theta <- function(theta, p, m) {
  stopifnot(length(theta) == m * (p + 1)+1)
  a <- theta[1]
  b <- theta[2:(m+1)]
  Xi_vec <- theta[(m + 2):(m * (p + 1)+1)] 
  Xi <- matrix(Xi_vec,nrow=p,ncol=m)
  list(a = a, b=b, Xi = Xi)
}

# list to vector
pack_theta <- function(a, b, Xi) c(a, b, as.vector(Xi))

# clipping function under l_infinity norm
clipping_inf  <- function(X, B) {
  return( pmax(-B, pmin(X, B)))  # Clip covariates elementwise to [-B, B]
} 

# clipping function under l_2 norm
clipping_2 <- function(X, gamma) {
  row_norms <- sqrt(rowSums(X^2))
  scales <- pmin(gamma / pmax(row_norms, 1e-12), 1)
  X_tilde <- X * scales
  return(X_tilde)
}

# fit non-linear model under given theta
predict_nn <- function(theta, X) { 
  p <- ncol(X)
  m <- as.integer( (length(theta)-1) / (p + 1))
  
  pars <- unpack_theta(theta, p = p, m = m)
  a  <- pars$a #  1
  b  <- pars$b # m*1
  Xi <- pars$Xi # p*m
  
  Z <- X %*% Xi # n*m
  H <- tanh(Z) # n*m
  re <- as.numeric(H %*% b + a) # n*1
  return(re)
}


########################################  
### ----------------------------------- 
### 2)  Non-private dnn procedure   
### ----------------------------------- 
######################################## 

# tilde theta  
compute_tilde_theta_np <- function(theta, X, Y,  eta0, tau, varpi, barK = barK_gaussian ) { 
  
  # common parameters
  n <- nrow(X)
  p <- ncol(X)
  m <- as.integer( (length(theta) -1) / (p + 1))
  
  # release (a,b,Xi) from long-vector theta
  pars <- unpack_theta(theta, p = p, m = m)
  a  <- pars$a              # 1
  b  <- pars$b              # m
  Xi <- pars$Xi             # p x m
  
  # ---- forward (vectorized) ----
  Z <- X %*% Xi             # (n x m) 
  H <- tanh(Z)              # (n x m)
  f <- as.numeric(H %*% b + a)  # (n)
  
  s <- barK((f - Y) / varpi) - tau      # (n)
  
  Hp <- 1 - H * H                       # tanh'(z) = 1 - tanh(z)^2  (n x m)
  AHp <- sweep(Hp, 2, b, "*")           # (n x m)
  V   <- sweep(AHp, 1, s, "*")          # (n x m)
  
  # ---- gradients (vectorized) ----
  grad_a  <- mean(s) 
  grad_b  <- colMeans(H * s)            # (m) 
  grad_Xi <- crossprod(X, V) / n        # (p x m)
  
  # ---- gradient step ----
  a_tilde  <- a  - eta0 * grad_a
  b_tilde  <- b  - eta0 * grad_b
  Xi_tilde <- Xi - eta0 * grad_Xi
  
  tilde_theta <- pack_theta(a_tilde, b_tilde, Xi_tilde)
  
  return(tilde_theta)
}

# fit T iterations
fit_nn_np <- function(X, Y, eta0, T_all, tau, varpi,  theta0) {
  theta <- theta0
  for (t in seq_len(T_all)) {
    tilde_theta <- compute_tilde_theta_np(theta, X, Y,   eta0, tau, varpi)
    theta <- clipping_inf(tilde_theta, 1000)
  }
  return(theta)
}


########################################### 
### --------------------------------------  
### 3)  (epsilon,delta)-DP dnn procedure  
### -------------------------------------- 
########################################### 

# tilde theta  
compute_tilde_theta_priv <- function(theta, X, Y,  eta0, tau, varpi, gamma, barK = barK_gaussian) { 
  
  # common parameters
  n <- nrow(X)
  p <- ncol(X)
  m <- as.integer( (length(theta) -1) / (p + 1))
  
  # release (a,b,Xi) from long-vector theta
  pars <- unpack_theta(theta, p = p, m = m)
  a  <- pars$a              # 1
  b  <- pars$b              # m
  Xi <- pars$Xi             # p x m
  
  # ---- forward (vectorized) ----
  Z <- X %*% Xi             # (n x m) 
  H <- tanh(Z)              # (n x m)
  f <- as.numeric(H %*% b + a)  # (n)
  
  s <- barK((f - Y) / varpi) - tau      # (n)
  
  Hp <- 1 - H * H                       # tanh'(z) = 1 - tanh(z)^2  (n x m)
  AHp <- sweep(Hp, 2, b, "*")           # (n x m)
  V   <- sweep(AHp, 1, s, "*")          # (n x m)
  
  # ---- gradients (vectorized) ----
  grad_a  <- mean(s) 
  grad_b  <- colMeans(H * s)            # (m) 
  tildeX  <- clipping_2(X, gamma)
  grad_Xi <- crossprod(tildeX, V) / n        # (p x m)
  
  # ---- gradient step ----
  a_tilde  <- a  - eta0 * grad_a
  b_tilde  <- b  - eta0 * grad_b
  Xi_tilde <- Xi - eta0 * grad_Xi
  
  tilde_theta <- pack_theta(a_tilde, b_tilde, Xi_tilde)
  
  return(tilde_theta)
}


fit_nn_priv <- function(X, Y, eta0, T_all, tau, varpi, theta0, 
                         gamma, B, epsilon, delta # privacy parameter
){
  theta <- theta0
  n <- nrow(X)
  p <- ncol(X)
  m <- as.integer( (length(theta) -1) / (p + 1))
  for (t in seq_len(T_all)) {
    tilde_theta <- compute_tilde_theta_priv(theta, X, Y, eta0, tau, varpi, gamma) # intermediate iterate
    Delta <- max(tau,1-tau)*sqrt(1+m+m*B^{2}*gamma^{2})  
    sigma <- 2 * eta0 * Delta * T_all / (n * epsilon) * sqrt(2 * log(1.25 * T_all / delta))
    check_theta <- tilde_theta + rnorm(length(tilde_theta), mean = 0, sd = sigma)
    theta <- clipping_inf(check_theta, B)
  }
  return(theta)
}


