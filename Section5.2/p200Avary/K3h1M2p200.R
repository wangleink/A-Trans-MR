rm(list=ls())
# library(rstudioapi)
# setwd(dirname(getActiveDocumentContext()$path))
source("inv_cov.R")
source("NDestep2.R")
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
K <- 3
s <- 3

detect<-matrix(0,nrow = REP,ncol = M)
Tarscore<-detectsigma<-rep(NA,REP)
Sourcescore<-matrix(NA,nrow =REP ,ncol = M)

beta.target<-beta.tl<-beta.tlw<-beta.tldelta<-
beta.pool<-beta.poolw<-beta.pooldelta<-array(NA,dim = c(REP,p+1,K-1))

beta.pool2<-beta.tl2<-array(NA,dim = c(REP,p+1,K-1))

SEdebias<-beta.tldebias<-matrix(NA,nrow = REP,ncol = (K-1)*p)

cl <- makeCluster(ncore)
registerDoParallel(cl)
packages <- c("pemultinom", "glmnet","mvtnorm","MASS","foreach","doParallel","caret")
clusterExport(cl, varlist="packages")
clusterEvalQ(cl, lapply(packages, require, character.only=TRUE))



ALLresult <- foreach(run=1:REP, .combine=abind::abind, .multicombine=TRUE, .init=array(0,c(nrow=p+1,ncol=(K-1)*6+1,ndim=0))) %dopar% {
  set.seed(run)
  print(run)
  
  beta_coef <- matrix(0, nrow = p+1, ncol = K-1)
  index1<-c(1:s)
  index2<-c(1:s)+s
  index3<-c(1:s)+2*s
  beta_coef[1+index1, 1] <- rep(1,s)
  beta_coef[1+index2, 2] <- rep(1,s)
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
          w[de,k]=(-1)^{k}*4
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
          w[de,k]=(-1)^{k}*4
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
  
  
  
  ##########################TL2
  detectid<-sdet_feng(X0=X.target,Y0=y.target,XK=X.source,YK=y.source,M=M,c=0.5,fold=3)
  Tarscore<- detectid[[2]]
  Sourcescore<-detectid[[3]]
  detectsigma<- detectid[[4]]
  ddetectid<-detectid[[1]]
  id<-rep(0,M)
  tmp=1
  for(int in 1:M)
  {
    if(length(ddetectid)==0)
    {
      break
    }
    if(length(ddetectid)==1)
    {
      id[ddetectid]=1
      break
    }
    if(tmp>length(ddetectid)){break}
    if(ddetectid[tmp]==int){id[int]=1
    tmp=tmp+1}
  }
  detectid<-id
  
  pX<-X.target
  py<-y.target
  for (int in 1:M) {
    if(detectid[int]==1){
    pX<-rbind(pX,X.source[[int]])
    py<-c(py,y.source[,int])}
    
  }
  
  ##step1
  fitw2 <- cv.pemultinom(pX, py, ncores = 1)
  tl.w2 <- fitw2$beta.min
  Beta.tlw2<- fitw2$beta.min
  
  
  ##########debias step
  #lambdaseq<-seq(0.05*sqrt(log(p)/ntarget),0.5*sqrt(log(p)/ntarget),length.out=10)
  lambdaseq<-seq(0.005,0.5*sqrt(log(p)/ntarget),length.out=10)
  cv.fit2<-cv.Destep2(X.target,y.target,betaw = tl.w2,kfold = 5,lambdaseq=lambdaseq ,T1=50,T2=50)
  deresult2<-DeStep_feng2(X.target,y.target,betaw = tl.w2,lambda = cv.fit2[[2]],T1=50,T2=50)
  tl.delta2<-deresult2[[1]]
  Beta.tldelta2<-tl.delta2
  Beta.tl2<-Beta.tlw2+Beta.tldelta2
  
  DD<-rep(10,p+1)
  DD[1:M]<-detectid
  
  
  
  
  
  
  
  
  
  
  ###########################pooled
  pX<-X.target
  py<-y.target
  for (m in 1:M) {
    pX<-rbind(pX,X.source[[m]])
    py<-c(py,y.source[,m])
    
  }
  fitw <- cv.pemultinom(pX, py, ncores = 1)
  tl.w <- fitw$beta.min
  Beta.poolw<- fitw$beta.min
  
  lambdaseq<-seq(0.005,0.5*sqrt(log(p)/ntarget),length.out=10)
  cv.fit<-cv.Destep2(X.target,y.target,betaw = tl.w,kfold = 5,lambdaseq=lambdaseq ,T1=50,T2=50)
  
  deresult<-DeStep_feng2(X.target,y.target,betaw = tl.w,lambda = cv.fit[[2]],T1=50,T2=50)
  tl.delta<-deresult[[1]]
  Beta.pooldelta<-tl.delta
  Beta.pool<-Beta.poolw+Beta.pooldelta
  
  
  betahat<-matrix(NA,nrow=p+1,ncol=(K-1)*6+1)
  betahat[,1:2]<-Beta.target
  betahat[,3:4]<-Beta.tl
  betahat[,5:6]<-Beta.pool
  betahat[,7:8]<-Beta.tlw
  betahat[,9:10]<-Beta.poolw
  betahat[,11:12]<-Beta.tl2
  betahat[,13]<-DD
  return(betahat)
  
}##end run
stopCluster(cl)

dd<-matrix(NA,nrow = REP,ncol = M)

for (i in 1:REP) {
  beta.target[i,,]<-ALLresult[,1:2,i]
  beta.tl[i,,]<-ALLresult[,3:4,i]
  beta.pool[i,,]<-ALLresult[,5:6,i]
  beta.tlw[i,,]<-ALLresult[,7:8,i]
  beta.poolw[i,,]<-ALLresult[,9:10,i]
  beta.tl2[i,,]<-ALLresult[,11:12,i]
  dd[i,]<-ALLresult[1:M,13,i]
}


beta_coef <- matrix(0, nrow = p+1, ncol = K-1)
index1<-c(1:s)
index2<-c(1:s)+s
index3<-c(1:s)+2*s
beta_coef[1+index1, 1] <- rep(1,s)
beta_coef[1+index2, 2] <- rep(1,s)


####MSE and Bias
Bias.target<-Bias.tl<-Bias.tl2<-Bias.pool<-Bias.pool2<-array(NA,c(REP,p+1,K-1))
L2.target<-L2.tl<-L2.tl2<-L2.pool<-L2.pool2<-matrix(NA,nrow = REP,ncol = K-1)
L1.target<-L1.tl<-L1.tl2<-L1.pool<-L1.pool2<-matrix(NA,nrow = REP,ncol = K-1)
Lmax.target<-Lmax.tl<-Lmax.tl2<-Lmax.pool<-Lmax.pool2<-matrix(NA,nrow = REP,ncol = K-1)
###Bias and MSE
for (int in 1:REP) {
  Bias.target[int,,]<-beta.target[int,,]-beta_coef
  Bias.tl[int,,]<-beta.tl[int,,]-beta_coef
  Bias.tl2[int,,]<-beta.tl2[int,,]-beta_coef
  Bias.pool[int,,]<-beta.pool[int,,]-beta_coef
  Bias.pool2[int,,]<-beta.poolw[int,,]-beta_coef
  
  Lmax.target[int,]<-apply(abs(Bias.target[int,,]),2,max)
  Lmax.tl[int,]<-apply(abs(Bias.tl[int,,]),2,max)
  Lmax.tl2[int,]<-apply(abs(Bias.tl2[int,,]),2,max)
  Lmax.pool[int,]<-apply(abs(Bias.pool[int,,]),2,max)
  Lmax.pool2[int,]<-apply(abs(Bias.pool2[int,,]),2,max)
  
  L2.target[int,]<-colSums(Bias.target[int,,]^2)
  L2.tl[int,]<-colSums(Bias.tl[int,,]^2)
  L2.tl2[int,]<-colSums(Bias.tl2[int,,]^2)
  L2.pool[int,]<-colSums(Bias.pool[int,,]^2)
  L2.pool2[int,]<-colSums(Bias.pool2[int,,]^2)
}


table<-rbind(c(apply(Lmax.target,2,mean),apply(Lmax.tl,2,mean),apply(Lmax.tl2,2,mean),apply(Lmax.pool,2,mean),apply(Lmax.pool2,2,mean)),
             c(apply(L2.target,2,mean),apply(L2.tl,2,mean),apply(L2.tl2,2,mean),apply(L2.pool,2,mean),apply(L2.pool2,2,mean)) )
row.names(table)<-c("Lmaxerror","L2error")
colnames(table)<-c("targetK1","targetK2","tlK1","tlK2","tl2K1","tl2K2","poolK1","poolK2","pool2K1","pool2K2")
print(table)




save.image("K3h1M2p200.RData")




