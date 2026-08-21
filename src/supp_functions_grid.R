nv_cost_vec <- function(q, y, b, h) { b * pmax(y - q, 0) + h * pmax(q - y, 0) }

fit_binning_policy <- function(X, Y, n_bins = 3) {
  X <- as.data.frame(X)
   
  keep <- vapply(X, function(z) {
    z <- z[is.finite(z)]
    length(unique(z)) > 1
  }, logical(1))
  
  X <- X[, keep, drop = FALSE]
  p <- ncol(X)
  
  if (p == 0) stop("No usable columns left for binning.")
  
  if (length(n_bins) == 1) n_bins <- rep(n_bins, p)
   
  breaks_list <- vector("list", p)
  
  for (j in seq_len(p)) {
    qs <- quantile(
      X[[j]],
      probs = seq(0, 1, length.out = n_bins[j] + 1),
      na.rm = TRUE,
      type = 7
    )
    
    brks <- unique(as.numeric(qs))
     
    if (length(brks) < 2) {
      stop(sprintf("Column %s has too few distinct values for binning.", names(X)[j]))
    }
    
    brks[1] <- -Inf
    brks[length(brks)] <- Inf
    
    breaks_list[[j]] <- brks
  }
   
  train_bins <- lapply(seq_len(p), function(j) {
    cut(X[[j]], breaks = breaks_list[[j]], include.lowest = TRUE, labels = FALSE)
  })
  
  train_cell <- do.call(paste, c(train_bins, sep = "_"))
  
  list(
    breaks_list = breaks_list,
    train_cell = train_cell,
    Y = Y,
    keep = keep
  )
}


predict_binning_policy <- function(fit, X, q_grid, b, h) {
  X <- as.data.frame(X)
  X <- X[, fit$keep, drop = FALSE]
  
  p <- ncol(X)
  
  test_bins <- lapply(seq_len(p), function(j) {
    cut(X[[j]], breaks = fit$breaks_list[[j]], include.lowest = TRUE, labels = FALSE)
  })
  
  test_cell <- do.call(paste, c(test_bins, sep = "_"))
  
  yhat <- numeric(nrow(X))
  tau <- b / (b + h)
  global_q <- as.numeric(quantile(fit$Y, probs = tau, type = 1, na.rm = TRUE))
  
  for (i in seq_len(nrow(X))) {
    idx <- which(fit$train_cell == test_cell[i])
    
    if (length(idx) == 0) {
      yhat[i] <- global_q
      next
    }
    
    y_local <- fit$Y[idx]
    y_local <- y_local[is.finite(y_local)]
    
    if (length(y_local) == 0) {
      yhat[i] <- global_q
      next
    }
    
    costs <- sapply(q_grid, function(q) {
      mean(b * pmax(y_local - q, 0) + h * pmax(q - y_local, 0), na.rm = TRUE)
    })
    
    if (!any(is.finite(costs))) {
      yhat[i] <- global_q
    } else {
      yhat[i] <- q_grid[which.min(costs)]
    }
  }
  
  yhat
}



fit_knn_policy <- function(X_train, Y_train) {
  X_train <- as.matrix(X_train)
  center <- colMeans(X_train)
  scale_ <- apply(X_train, 2, sd)
  scale_[scale_ == 0] <- 1
  
  X_train_sc <- scale(X_train, center = center, scale = scale_)
  
  list(
    X_train_sc = X_train_sc,
    Y_train = Y_train,
    center = center,
    scale_ = scale_
  )
}


predict_knn_policy <- function(fit, X_test, k, q_grid, b, h) {
  X_test <- as.matrix(X_test)
  X_test_sc <- scale(X_test, center = fit$center, scale = fit$scale_)
  
  yhat <- numeric(nrow(X_test_sc))
  
  for (i in seq_len(nrow(X_test_sc))) {
    dists <- rowSums((fit$X_train_sc - matrix(X_test_sc[i, ], nrow = nrow(fit$X_train_sc), ncol = ncol(fit$X_train_sc), byrow = TRUE))^2)
    nn_idx <- order(dists)[1:k]
    y_local <- fit$Y_train[nn_idx]
    
    costs <- sapply(q_grid, function(q) mean(nv_cost_vec(q, y_local, b, h)))
    yhat[i] <- q_grid[which.min(costs)]
  }
  
  yhat
}