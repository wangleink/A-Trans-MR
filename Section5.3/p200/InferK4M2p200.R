rm(list=ls())
# library(rstudioapi)
# setwd(dirname(getActiveDocumentContext()$path))
source("noin-invcov2.R")
source("noin-invcov.R")
source("NDestep.R")
library(pemultinom)
library(mvtnorm)
library(MASS)
library(glmnet)
library(foreach)
library(doParallel)
library(caret)

REP=200
ncore=20

ntarget<-nsource<- 200
M<-8
M_h<-2


h=1

exact=0


p <- 200
K <- 4
s <- 3

detect<-matrix(0,nrow = REP,ncol = M)
Tarscore<-detectsigma<-rep(NA,REP)
Sourcescore<-matrix(NA,nrow =REP ,ncol = M)


cl <- makeCluster(ncore)
registerDoParallel(cl)
packages <- c("pemultinom", "glmnet","mvtnorm","MASS","foreach","doParallel","caret")
clusterExport(cl, varlist="packages")
clusterEvalQ(cl, lapply(packages, require, character.only=TRUE))


ALLresult <- foreach(run=1:REP, .combine=abind::abind, .multicombine=TRUE, .init=array(0,c(nrow=p*(K-1),ncol=2,ndim=0))) %dopar% {
  set.seed(run)
  print(run)
  
  beta_coef <- matrix(0, nrow = p+1, ncol = K-1)
  index1<-c(1:s)
  index2<-c(1:s)+s
  index3<-c(1:s)+2*s
  beta_coef[1+index1, 1] <- rep(1,s)
  beta_coef[1+index2, 2] <- rep(1,s)
  beta_coef[1+index3, 3] <- rep(1,s)
  coef<-list()
  if(exact==1){
    coef[[1]]<-beta_coef
    for (m in 1:M) {
      w<-beta_coef
      for (k in 1:(K-1)) {
        sign0<-sample(c(-1,1),s,replace = T)
        sign<-c(0,sign0,rep(0,p-s))
        w[,k]<-beta_coef[,k]+h/p*sign
      }
      coef[[m+1]]<-w
      
    }
  }
  if(exact==0){
    if(M_h==0){
      coef[[1]]<-beta_coef
      for (m in 1:M) {
        w<-beta_coef
        for (k in 1:(K-1)) {
          sign0<-sample(c(-1,1),p,replace = T)
          sign<-c(0,sign0)
          #de<-sample(((K-1)*s):p,s,replace = F)+1
          de<-c(((K-1)*s+1):(K*s))+1
          w[,k]<-beta_coef[,k]+2*h/p*sign
          w[de,k]=w[de,k]=(-1)^{k}*4
        }
        coef[[m+1]]<-w
        
      }
    }else{
      coef[[1]]<-beta_coef
      for (m in 1:(M_h)) {
        w<-beta_coef
        for (k in 1:(K-1)) {
          sign0<-sample(c(-1,1),p,replace = T)
          sign<-c(0,sign0)
          w[,k]<-beta_coef[,k]+h/p*sign
        }
        coef[[m+1]]<-w
        
      }
      for (m in (M_h+1):M) {
        w<-beta_coef
        for (k in 1:(K-1)) {
          sign0<-sample(c(-1,1),p,replace = T)
          sign<-c(0,sign0)
          #de<-sample(((K-1)*s):p,s,replace = F)+1
          de<-c(((K-1)*s+1):(K*s))+1
          w[,k]<-beta_coef[,k]+2*h/p*sign
          w[de,k]=w[de,k]=(-1)^{k}*4
        }
        coef[[m+1]]<-w
        
      }
    }
    
    
    
  }
  
  
  sigmat<-diag(p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigmat[i,j]=0.5^(abs(i-j))
    }
  }
  sigmas<-diag(p)
  eps<-rnorm(p,sd=0.3)
  sigmas<-sigmat+eps%*%t(eps)
  
  
  
  
  ##
  # sum(c(y.target==3))
  X.target<-mvrnorm(n=ntarget,mu=rep(0,p),Sigma=sigmat)
  y.target <- sapply(1:ntarget, function(j){
    prob_i <- c(sapply(1:(K-1), function(k){
      exp(sum(c(1,X.target[j, ])*beta_coef[,k]))
    }), 1)
    prob_i <- prob_i/sum(prob_i)
    sample(1:K, size = 1, replace = TRUE, prob = prob_i)
  })
  
  X.source<-list()
  y.source<-matrix(NA,nrow=nsource,ncol =M)
  
  for (m in 1:M) {
    X<-mvrnorm(n=nsource,mu=rep(0,p),Sigma=sigmas)
    y<- sapply(1:nsource, function(j){
      prob_i <- c(sapply(1:(K-1), function(k){
        exp(sum(c(1,X[j, ])*coef[[1+m]][,k]))
      }), 1)
      prob_i <- prob_i/sum(prob_i)
      sample(1:K, size = 1, replace = TRUE, prob = prob_i)
    })
    X.source[[m]]<-X
    y.source[,m]<-y
  }
  
  
  ######################################point estimation
  
  ###########fit target
  fit <- cv.pemultinom(X.target, y.target, ncores = 1)
  Beta.target<-fit$beta.min
  
  #############fit TL
  #transfer step
  if(M_h==0){
    pX<-X.target
    py<-y.target
  }else{
    pX<-X.target
    py<-y.target
    for (m in 1:M_h) {
      pX<-rbind(pX,X.source[[m]])
      py<-c(py,y.source[,m])
      
    }
    
  }
  fitw <- cv.pemultinom(pX, py, ncores = 1)
  tl.w <- fitw$beta.min
  Beta.tlw<- fitw$beta.min
  
  
  ##########debias step
  #lambdaseq<-seq(0.05*sqrt(log(p)/ntarget),0.5*sqrt(log(p)/ntarget),length.out=10)
  lambdaseq<-seq(0.005,0.5*sqrt(log(p)/ntarget),length.out=10)
  cv.fit<-cv.Destep2(X.target,y.target,betaw = tl.w,kfold = 5,lambdaseq=lambdaseq ,T1=50,T2=50)
  
  deresult<-DeStep_feng2(X.target,y.target,betaw = tl.w,lambda = cv.fit[[2]],T1=50,T2=50)
  tl.delta<-deresult[[1]]
  Beta.tldelta<-tl.delta
  Beta.tl<-Beta.tlw+Beta.tldelta
  
  
  ###############################
  
  ###############################  Inference
  pb_calc <- function(X, beta) {
    expxb <- cbind(exp(cbind(1, X) %*% beta), 1)
    pb <- expxb/rowSums(expxb)
  }
  p.big=(K-1)*p
  
  ##########
  #Inv_TL
  Inv_tl1<-inv_cov_calc2_1(x=pX,y=py,beta=Beta.tl, nfolds = 4,lambda.choice = "lambda.min") 
  Theta.1<-Inv_tl1[[1]]
  Inv_tl2<-inv_cov_calc2_2(x=X.target,ThetaW=Theta.1,beta=Beta.tl, nfolds = 4,lambda.choice = "lambda.min") 
  Theta.2<-Inv_tl2[[1]]
  
  p.big=(K-1)*p
  Theta.int<-matrix(1,nrow = p.big, ncol = p.big)
  for (j in 1:p.big) {
    Theta.int[j,-j]<-Theta.2[j,-j]+Theta.1[j,-j]
  }
  ###
  pb.T <- pb_calc(X=pX, beta =Beta.tl)
  Sigma.T <- matrix(nrow = (K-1)*p, ncol=(K-1)*p)
  for (k1 in 1:(K-1)) {
    for (k2 in 1:(K-1)) {
      if (k1 == k2) {
        Sigma.T[(1+(k1-1)*p):(k1*p), (1+(k2-1)*p):(k2*p)] <- (t(pX) %*% diag(pb.T[,k1]*(1-pb.T[,k1])) %*% pX)/(length(py))
      } else {
        Sigma.T[(1+(k1-1)*p):(k1*p), (1+(k2-1)*p):(k2*p)] <- (t(pX) %*% diag(-pb.T[,k1]*pb.T[,k2]) %*% pX)/(length(py))
      }
    }
  }
  Theta <- matrix(nrow = p.big, ncol = p.big)
  tau2 <- sapply(1:p.big, function(j){
    Sigma.T[j, j] - Sigma.T[j,-j] %*% Theta.int[j,-j ]
  })
  for (j in 1:p.big) {
    if (tau2[j] == 0){
      Theta[j, ] <- 0
    } else {
      Theta[j, j] <- 1/tau2[j]
      Theta[j, -j] <- -Theta.int[j,-j]/tau2[j]
    }
  }
  Theta.tl<-Theta
  
  beta.ini<-Beta.tl[-1,1]
  for (int in 2:(K-1)) {
    beta.ini<-c(beta.ini,Beta.tl[-1,int])
  }
  pb<-pb_calc(X=X.target, beta = Beta.tl)
  id=c(y.target==1)
  yk<-rep(0,ntarget)
  yk[id]<-1
  D<-t(X.target)%*%(yk-pb[,1])
  for (int in 2:(K-1)) {
    id=c(y.target==int)
    yk<-rep(0,ntarget)
    yk[id]<-1
    D<-c(D,t(X.target)%*%(yk-pb[,int]))
  }
  dterm<-Theta.tl%*%D/(ntarget)
  beta.debias<-beta.ini+dterm
  
  beta<-matrix(NA,nrow = p+1,ncol = K-1)
  for (int in 1:(K-1)) {
    beta[,int]<-c(0,beta.debias[((int-1)*p+1):(int*p)])
  }
  
  pb<-pb_calc(X=pX, beta = Beta.tl)
  Sigma.full <- matrix(NA,nrow = p.big, ncol=p.big)
  for (k1 in 1:(K-1)) {
    for (k2 in 1:(K-1)) {
      if (k1 == k2) {
        Sigma.full[(1+(k1-1)*p):(k1*p), (1+(k2-1)*p):(k2*p)] <- (t(pX) %*% diag(pb[,k1]*(1-pb[,k1])) %*% pX)/(length(py))
      } else {
        Sigma.full[(1+(k1-1)*p):(k1*p), (1+(k2-1)*p):(k2*p)] <- (t(pX) %*% diag(-pb[,k1]*pb[,k2]) %*% pX)/(length(py))
      }
    }
  }
  sedebias.tl<-sqrt(diag(Theta.tl%*%Sigma.full%*%Theta.tl/ntarget))
  DBeta.tl<-beta.debias
  
  
  
  
  
  
  R<-matrix(NA,nrow=p*(K-1),ncol=2)
  R[,1]<-sedebias.tl
  R[,2]<-DBeta.tl
  
  return(R)
  
}##end run
stopCluster(cl)


###############################    Inference result     #####################

SEdebias.tar<-Dbeta.tar<-SEdebias.tlw<-Dbeta.tlw<-SEdebias.tl<-Dbeta.tl<-matrix(NA,nrow = REP,ncol = (K-1)*p)


for (i in 1:REP) {
  SEdebias.tl[i,]<-ALLresult[,1,i]
  Dbeta.tl[i,]<-ALLresult[,2,i]
}


beta_coef <- matrix(0, nrow = p+1, ncol = K-1)
index1<-c(1:s)
index2<-c(1:s)+s
index3<-c(1:s)+2*s
beta_coef[1+index1, 1] <- rep(1,s)
beta_coef[1+index2, 2] <- rep(1,s)
beta_coef[1+index3, 3] <- rep(1,s)

beta.true<-beta_coef[-1,1]
for (int in 2:(K-1)) {
  beta.true<-c(beta.true,beta_coef[-1,int])
}
########################################CP CI

C.tar<-C.tlw<-C.tl<-matrix(NA,nrow = REP,ncol = (K-1)*p)
for(i in 1:REP)
{
  C.tl[i,]<-as.numeric(c(Dbeta.tl[i,]<beta.true+qnorm(0.975)*SEdebias.tl[i,] & Dbeta.tl[i,]>beta.true-qnorm(0.975)*SEdebias.tl[i,]))
}

CP.tl<-apply(C.tl,2,mean)


####SD SE

SD.tl<-apply(Dbeta.tl,2,sd)
SE.tl<-apply(SEdebias.tl,2,mean)

result=matrix(c(mean(CP.tl),
                mean(SE.tl),
                mean(SD.tl)),3,1,byrow=T)
rownames(result)=c("CP","SE","SD")
colnames(result)=c("tl")

result

save.image("InferK4M2p200.RData")

