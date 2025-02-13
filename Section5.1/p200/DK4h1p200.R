rm(list=ls())
source("inv_cov.R")
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
h=1

exact=0


p <- 200
K <- 4
s <- 3

detect<-matrix(0,nrow = REP,ncol = M)
Tarscore<-detectsigma<-rep(NA,REP)
Sourcescore<-matrix(NA,nrow =REP ,ncol = M)

beta.target<-beta.tl<-beta.tlw<-beta.tldelta<-
  beta.pool<-beta.poolw<-beta.pooldelta<-array(NA,dim = c(REP,p+1,K-1))
SEdebias<-beta.tldebias<-matrix(NA,nrow = REP,ncol = (K-1)*p)

cl <- makeCluster(ncore)
registerDoParallel(cl)
packages <- c("pemultinom", "glmnet","mvtnorm","MASS","foreach","doParallel","caret")
clusterExport(cl, varlist="packages")
clusterEvalQ(cl, lapply(packages, require, character.only=TRUE))



ALLresult <- foreach(run=1:REP, .combine=abind::abind, .multicombine=TRUE, .init=array(0,c(nrow=4,ncol=M,ndim=0))) %dopar% {
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
    coef[[1]]<-beta_coef
    for (m in 1:(ceiling(M/2))) {
      w<-beta_coef
      for (k in 1:(K-1)) {
        sign0<-sample(c(-1,1),p,replace = T)
        sign<-c(0,sign0)
        w[,k]<-beta_coef[,k]+h/p*sign
      }
      coef[[m+1]]<-w
      
    }
    for (m in (ceiling(M/2)+1):M) {
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
  
  
  # MM<-matrix(NA,nrow = M,ncol = K) 
  # 
  # for (i in 1:M) {
  #   MM[i,]<-c(sum(y.source[,i]==1),sum(y.source[,i]==2),sum(y.source[,i]==3))
  # }
  # 
  # print(MM)
  # c(sum(y.target==1),sum(y.target==2),sum(y.target==3))
  
  
  
  ######################################detection
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
  
  
  
  DD<-matrix(NA,nrow=4,ncol=M)
  DD[1,]<-detectid
  DD[2,]<-Sourcescore
  DD[3,]<-c(Tarscore,rep(0,M-1))
  DD[4,]<-c(detectsigma,rep(0,M-1))
  return(DD)
  
  
}##end run
stopCluster(cl)

for (i in 1:REP) {
  detect[i,]<-ALLresult[1,,i]
  Sourcescore[i,]<-ALLresult[2,,i]
  Tarscore[i]<-ALLresult[3,1,i]
  detectsigma[i]<-ALLresult[4,1,i]
}


Sourcescore
Tarscore
detectsigma


AVEdetect<-apply(detect,2,mean)
print(AVEdetect)

save.image("DK4h1p200.RData")






