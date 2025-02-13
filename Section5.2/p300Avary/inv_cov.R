source("RcppExports.R")
penalized_quad_r <- function(A, b, lambda_list, max_iter=200, tol=1e-3) {
  penalized_quad(A, b, lambda_list, max_iter, tol)
}


# x=X.target
# K=3
# beta=beta.target
# nfolds = 5
# lambda.choice = "lambda.min"
# nlambda = 100
# max_iter = 200
# tol = 1e-3
# ncores = 1



inv_cov_calc<- function(x, K, beta, nfolds = 5, nlambda = 100, max_iter = 200, tol = 1e-3, ncores = 1,
                        lambda.choice = "lambda.min") {
  pb_calc <- function(X, beta) {
    expxb <- cbind(exp(cbind(1, X) %*% beta), 1)
    pb <- expxb/rowSums(expxb)
    pb
  }
  
  p <- ncol(x)
  n <- nrow(x)
  p.big <- (K-1)*p
  fold <- sample(1:n %% nfolds)+1
  
  pb <- pb_calc(X=x, beta = beta)
  Sigma.full <- matrix(nrow = (K-1)*p, ncol=(K-1)*p)
  for (k1 in 1:(K-1)) {
    for (k2 in 1:(K-1)) {
      if (k1 == k2) {
        Sigma.full[(1+(k1-1)*p):(k1*p), (1+(k2-1)*p):(k2*p)] <- (t(x) %*% diag(pb[,k1]*(1-pb[,k1])) %*% x)/n
      } else {
        Sigma.full[(1+(k1-1)*p):(k1*p), (1+(k2-1)*p):(k2*p)] <- (t(x) %*% diag(-pb[,k1]*pb[,k2]) %*% x)/n
      }
    }
  }
  lambda.max <- sapply(1:p.big, function(j){
    max(Sigma.full[j, -j])
  })
  
  
  loss <- sapply(1:nfolds, function(r){ 
    # message(paste("Fold ", r, " of ", nfolds ,"......", sep = ""))
    
    pb <- pb_calc(X=x[fold != r, ], beta = beta)
    Sigma <- matrix(nrow = (K-1)*p, ncol=(K-1)*p)
    for (k1 in 1:(K-1)) {
      for (k2 in 1:(K-1)) {
        if (k1 == k2) {
          Sigma[(1+(k1-1)*p):(k1*p), (1+(k2-1)*p):(k2*p)] <- (t(x[fold != r, ]) %*% diag(pb[,k1]*(1-pb[,k1])) %*% x[fold != r, ])/sum(fold != r)
        } else {
          Sigma[(1+(k1-1)*p):(k1*p), (1+(k2-1)*p):(k2*p)] <- (t(x[fold != r, ]) %*% diag(-pb[,k1]*pb[,k2]) %*% x[fold != r, ])/sum(fold != r)
        }
      }
    }
    Sigma.valid <- matrix(nrow = (K-1)*p, ncol=(K-1)*p)
    pb <- pb_calc(X=x[fold== r, ], beta = beta)
    for (k1 in 1:(K-1)) {
      for (k2 in 1:(K-1)) {
        if (k1 == k2) {
          Sigma.valid[(1+(k1-1)*p):(k1*p), (1+(k2-1)*p):(k2*p)] <- (t(x[fold == r, ]) %*% diag(pb[,k1]*(1-pb[,k1])) %*% x[fold == r, ])/sum(fold == r)
        } else {
          Sigma.valid[(1+(k1-1)*p):(k1*p), (1+(k2-1)*p):(k2*p)] <- (t(x[fold == r, ]) %*% diag(-pb[,k1]*pb[,k2]) %*% x[fold == r, ])/sum(fold == r)
        }
      }
    }
    
    loss <- foreach(j = 1:p.big, .combine = "rbind",.packages = 'pemultinom') %dopar% {
      if (lambda.max[j] == 0) {
        rep(NA, nlambda)
      } else {
        lambda.min <- 0.01*lambda.max[j]
        lambda_list <- exp(seq(log(lambda.max[j]), log(lambda.min), len = nlambda))
        gamma_j_matrix <- penalized_quad_r(A = Sigma[-j, -j], b = Sigma[j, -j], lambda_list = lambda_list, max_iter = max_iter, tol = tol)
        sapply(1:nlambda, function(l){
          Sigma.valid[j,j] - Sigma.valid[j, -j] %*% gamma_j_matrix[l, ] + gamma_j_matrix[l, ] %*% Sigma.valid[-j, -j] %*% gamma_j_matrix[l, ]/2
        })
      }
    }
    loss
  }, simplify = FALSE)
  message("Determining the optimal lambda and finalizing the inference results......")
  # calculate gamma for each j with the full data and best lambda
  gamma_matrix <- foreach(j = 1:p.big, .combine = "rbind") %dopar% {
    if (lambda.max[j] == 0) {
      rep(0, p.big - 1)
    } else {
      lambda.min <- 0.01*lambda.max[j]
      lambda_list <- exp(seq(log(lambda.max[j]), log(lambda.min), len = nlambda))
      gamma_j_matrix <- penalized_quad_r(A = Sigma.full[-j, -j], b = Sigma.full[j, -j], lambda_list = lambda_list, max_iter = max_iter, tol = tol)
      loss_j <- sapply(1:nfolds, function(r){
        loss[[r]][j,]
      })
      cvm <- rowMeans(loss_j, na.rm = TRUE)
      cvsd <- apply(loss_j, 1, function(x){sd(x, na.rm = TRUE)})
      ind.min <- which.min(cvm)
      lambda.min <- lambda_list[ind.min]
      cvsd.min <- cvsd[ind.min]
      cvm.min <- cvm[ind.min]
      ind.1se <- min(which(cvm <= cvm.min + cvsd.min))
      lambda.1se <- lambda_list[ind.1se]
      if(lambda.choice == "lambda.1se") {
        gamma_j <- gamma_j_matrix[ind.1se, ]
      } else if(lambda.choice == "lambda.min") {
        gamma_j <- gamma_j_matrix[ind.min, ]
      }
      gamma_j
    }
  }
  
  tau2 <- sapply(1:p.big, function(j){
    Sigma.full[j, j] - Sigma.full[j,-j] %*% gamma_matrix[j, ]
  })
  
  Theta <- matrix(nrow = p.big, ncol = p.big)
  for (j in 1:p.big) {
    if (tau2[j] == 0){
      Theta[j, ] <- 0
    } else {
      Theta[j, j] <- 1/tau2[j]
      Theta[j, -j] <- -gamma_matrix[j, ]/tau2[j]
    }
    
  }
  
  return(list(Theta = Theta))
}
