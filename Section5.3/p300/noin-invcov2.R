source("RcppExports.R")

penalized_quad_r <- function(A, b, lambda_list, max_iter=200, tol=1e-3) {
  penalized_quad(A, b, lambda_list, max_iter, tol)
}

inv_cov_calc2_1<- function(x, y, beta, nfolds = 5, nlambda = 100, max_iter = 200, tol = 1e-3, ncores = 1,
                         lambda.choice = c("lambda.1se", "lambda.min")) {
  
  pb_calc <- function(X, beta) {
    expxb <- cbind(exp(cbind(1, X) %*% beta), 1)
    pb <- expxb/rowSums(expxb)
    pb
  }
  
  p <- ncol(x)
  n <- nrow(x)
  K <- length(unique(y))
  p.big <- (K-1)*p
  
  n_list <- as.numeric(table(y))
  fold_by_class <- sapply(1:K, function(k){
    lb <- sample(1:n_list[k] %% nfolds)
    lb[lb == 0] <- nfolds
    lb
  }, simplify = F)
  fold <- numeric(n)
  for (k in 1:K) {
    fold[y == k] <- fold_by_class[[k]]
  }
  
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
    message(paste("Fold ", r, sep = ""))
    
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
    pb <- pb_calc(X=x[fold == r, ], beta = beta)
    for (k1 in 1:(K-1)) {
      for (k2 in 1:(K-1)) {
        if (k1 == k2) {
          Sigma.valid[(1+(k1-1)*p):(k1*p), (1+(k2-1)*p):(k2*p)] <- (t(x[fold == r, ]) %*% diag(pb[,k1]*(1-pb[,k1])) %*% x[fold == r, ])/sum(fold == r)
        } else {
          Sigma.valid[(1+(k1-1)*p):(k1*p), (1+(k2-1)*p):(k2*p)] <- (t(x[fold == r, ]) %*% diag(-pb[,k1]*pb[,k2]) %*% x[fold == r, ])/sum(fold == r)
        }
      }
    }
    
    CVLoss<-matrix(NA,nrow =p.big ,ncol = nlambda)
    for (j in 1:p.big) {
      if (lambda.max[j] == 0) {
        rep(NA, nlambda)
      } else {
        lambda.min <- 0.01*lambda.max[j]
        lambda_list <- exp(seq(log(lambda.max[j]), log(lambda.min), len = nlambda))
        gamma_j_matrix <- penalized_quad_r(A = Sigma[-j, -j], b = Sigma[j, -j], lambda_list = lambda_list, max_iter = max_iter, tol = tol)
        CVLoss[j,]<-sapply(1:nlambda, function(l){
          Sigma.valid[j,j] - Sigma.valid[j, -j] %*% gamma_j_matrix[l, ] + gamma_j_matrix[l, ] %*% Sigma.valid[-j, -j] %*% gamma_j_matrix[l, ]/2
        })
      }
      
    }
    CVLoss
  }, simplify = FALSE)
  
  message("Determining the optimal lambda and finalizing the inference results......")
  # calculate gamma for each j with the full data and best lambda
  gamma_matrix<-matrix(NA,nrow =p.big ,ncol = p.big-1)
  for (j in 1:p.big) {
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
      gamma_matrix[j,]<-gamma_j
    }
  }
  
  Theta <- matrix(nrow = p.big, ncol = p.big)
  for (j in 1:p.big) {
      Theta[j, j] <-1
      Theta[j, -j] <-gamma_matrix[j, ]
    
  }
  
  return(list(Theta = Theta))
}

inv_cov_calc2_2<-function(x, beta,ThetaW ,nfolds = 4, nlambda = 100, max_iter = 200, tol = 1e-3,
                          lambda.choice = "lambda.min") {
  pb_calc <- function(X, beta) {
    expxb <- cbind(exp(cbind(1, X) %*% beta), 1)
    pb <- expxb/rowSums(expxb)
    pb
  }
  
  p <- ncol(x)
  n <- nrow(x)
  K <- length(unique(y))
  p.big <- (K-1)*p
  
  n_list <- as.numeric(table(y))
  fold_by_class <- sapply(1:K, function(k){
    lb <- sample(1:n_list[k] %% nfolds)
    lb[lb == 0] <- nfolds
    lb
  }, simplify = F)
  fold <- numeric(n)
  for (k in 1:K) {
    fold[y == k] <- fold_by_class[[k]]
  }
  
  
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
    message(paste("Fold ", r, " of ", nfolds ,"......", sep = ""))
    
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
    
    CVLoss<-matrix(NA,nrow =p.big ,ncol = nlambda)
    for (j in 1:p.big) {
      if (lambda.max[j] == 0) {
        rep(NA, nlambda)
      } else {
        lambda.min <- 0.01*lambda.max[j]
        lambda_list <- exp(seq(log(lambda.max[j]), log(lambda.min), len = nlambda))
        gamma_j_matrix <- penalized_quad_r(A = Sigma[-j, -j], b = Sigma[j, -j]-t(ThetaW[j,-j])%*%Sigma[-j, -j], lambda_list = lambda_list, max_iter = max_iter, tol = tol)
        CVLoss[j,]<-sapply(1:nlambda, function(l){
          Sigma.valid[j,j]-( Sigma.valid[j, -j]-t(ThetaW[j,-j])%*%Sigma.valid[-j, -j]) %*% gamma_j_matrix[l, ] + gamma_j_matrix[l, ] %*% Sigma.valid[-j, -j] %*% gamma_j_matrix[l, ]/2
        })
      }
      
    }
    
    CVLoss
  }, simplify = FALSE)
  
  
  message("Determining the optimal lambda and finalizing the inference results......")
  # calculate gamma for each j with the full data and best lambda
  gamma_matrix<-matrix(NA,nrow =p.big ,ncol = p.big-1)
  for (j in 1:p.big) {
    if (lambda.max[j] == 0) {
      rep(0, p.big - 1)
    } else {
      lambda.min <- 0.01*lambda.max[j]
      lambda_list <- exp(seq(log(lambda.max[j]), log(lambda.min), len = nlambda))
      gamma_j_matrix <- penalized_quad_r(A = Sigma.full[-j, -j], b = Sigma.full[j, -j]-t(ThetaW[j,-j])%*%Sigma.full[-j, -j], lambda_list = lambda_list, max_iter = max_iter, tol = tol)
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
      gamma_matrix[j,]<-gamma_j
    }
  }

  Theta <- matrix(nrow = p.big, ncol = p.big)
  for (j in 1:p.big) {
      Theta[j, j] <- 1
      Theta[j, -j] <-gamma_matrix[j, ]
    }
    
  return(list(Theta = Theta))
}



